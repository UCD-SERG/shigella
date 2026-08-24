#!/usr/bin/env Rscript
# =====================================================================
# make_banks.R -- turn the chapter 2 fits into the two parameter banks that
# chapter 3's likelihood consumes.
#
#   arm T  cross-covariances estimated   -> correlated bank
#   arm S  cross-covariances fixed at 0  -> independent bank (= chapter 1)
#
# Both are `par_new` draws, not `mu_par`. par_new is the predictive
# distribution for a participant who has not been seen; mu_par is the
# population mean. Substituting the mean would assert that everyone has the
# average trajectory, which makes one titre far more informative about the time
# since infection than it is.
#
# THREE THINGS THIS SCRIPT REFUSES TO GUESS
#
#   1  the layout of par_new. Stan flattens `vector[10]` as par_new[1..10] but
#      `matrix[2,5]` column-major as par_new[1,1], par_new[2,1], ... which
#      interleaves the isotypes. The script detects which and then checks the
#      result against the population medians; a scrambled mapping fails loudly
#      instead of producing a plausible-looking bank.
#   2  which chains to keep. The lp rule from chapter 2 is applied, and the
#      number kept is printed. arm S has no minor mode -- its lp spread is 71,
#      inside the threshold -- so all its chains are retained, and that is
#      reported rather than assumed.
#   3  whether the correlation survived. cor(log y1) and the rest are printed
#      for both banks. arm T should show the chapter 2 values; arm S near zero.
#      If arm T comes out near zero the `iter` pairing is broken and nothing
#      downstream is meaningful.
#
#   Rscript R/02_make_banks.R
#   N_DRAWS=1000 Rscript R/02_make_banks.R
#   CH2_DIR=~/"chapter 2 work"/sees/ch2 Rscript R/02_make_banks.R
# =====================================================================

suppressPackageStartupMessages(library(posterior))

LP_DROP <- as.numeric(Sys.getenv("LP_DROP", "100"))
N_DRAWS <- as.integer(Sys.getenv("N_DRAWS", "500"))
OUT <- Sys.getenv("BANK_DIR", "data")
CH2 <- Sys.getenv("CH2_DIR", "~/chapter 2 work/sees/ch2")
ARMS <- list(T = "correlated", S = "independent")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(..., "\n", sep = "")
hr <- function() cat(strrep("-", 74), "\n")

# Population medians from the chapter 2 arm T fit, used only as a sanity check
# on the par_new LAYOUT -- if the isotypes were swapped or interleaved the
# values would be wildly wrong. They are not a target: arm S has its own
# population parameters and is expected to differ by tens of percent. The
# tolerance below is deliberately loose for that reason.
EXPECT <- list(
  HlyE_IgG = c(y0 = 22.6, y1 = 148.3, t1 = 22.5, r = 2.77, kdec = 0.0492),
  HlyE_IgA = c(y0 = 15.0, y1 = 51.7, t1 = 16.2, r = 2.89, kdec = 0.2289)
)

# =====================================================================
extract_par_new <- function(fit) {
  v <- fit$metadata()$stan_variables
  if (!"par_new" %in% v) stop("par_new not in this fit", call. = FALSE)
  pn <- as.matrix(fit$draws("par_new", format = "draws_matrix"))
  cn <- colnames(pn)

  # layout A: par_new[1] ... par_new[10]        -> IgG 1:5, IgA 6:10
  # layout B: par_new[1,1] ... par_new[2,5]     -> [isotype, parameter]
  if (all(grepl("^par_new\\[[0-9]+\\]$", cn))) {
    idx <- list(
      HlyE_IgG = sprintf("par_new[%d]", 1:5),
      HlyE_IgA = sprintf("par_new[%d]", 6:10)
    )
    layout <- "vector[10]"
  } else if (any(grepl("^par_new\\[[0-9]+,[0-9]+\\]$", cn))) {
    idx <- list(
      HlyE_IgG = sprintf("par_new[1,%d]", 1:5),
      HlyE_IgA = sprintf("par_new[2,%d]", 1:5)
    )
    layout <- "matrix[2,5]"
  } else {
    say("unrecognised par_new naming:")
    print(head(cn, 12))
    stop(call. = FALSE)
  }
  miss <- setdiff(unlist(idx), cn)
  if (length(miss)) {
    say("missing columns: ", paste(miss, collapse = " "))
    stop(call. = FALSE)
  }
  list(draws = pn, idx = idx, layout = layout)
}

to_curve_params <- function(fit, n_draws, lp_drop) {
  lp <- apply(fit$draws("lp__", format = "draws_array"), 2, median)
  K <- length(lp)
  keep <- which(lp > max(lp) - lp_drop)
  say(sprintf(
    "  chains %d kept of %d   (lp spread %.0f, threshold %g)",
    length(keep), K, max(lp) - min(lp), lp_drop
  ))
  if (length(keep) < K) {
    d <- setdiff(seq_len(K), keep)
    say(sprintf(
      "    dropped: %s",
      paste(sprintf("chain %d (%.0f below)", d, lp[d] - max(lp)), collapse = "; ")
    ))
  } else {
    say("    no chain falls outside the threshold")
  }

  P <- extract_par_new(fit)
  say("  par_new layout: ", P$layout)
  nit <- nrow(P$draws) / K
  rows <- unlist(lapply(keep, function(k) ((k - 1) * nit + 1):(k * nit)))
  rows <- rows[round(seq(1, length(rows), length.out = min(n_draws, length(rows))))]

  # log scale -> natural scale. The fourth coordinate is log k_dec, so alpha
  # has to be recovered: k_dec = alpha * y1^(r-1).
  out <- do.call(rbind, lapply(names(P$idx), function(nm) {
    p <- P$draws[rows, P$idx[[nm]], drop = FALSE]
    storage.mode(p) <- "double"
    dimnames(p) <- NULL
    y0 <- exp(p[, 1])
    y1 <- y0 + exp(p[, 2])
    t1 <- exp(p[, 3])
    kd <- exp(p[, 4])
    r <- exp(p[, 5]) + 1
    data.frame(
      antigen_iso = nm, iter = seq_along(rows),
      y0 = y0, y1 = y1, t1 = t1,
      alpha = kd / y1^(r - 1), r = r,
      row.names = NULL
    )
  }))
  attr(out, "chains_kept") <- keep
  attr(out, "layout") <- P$layout
  out
}

check_layout <- function(cp) {
  ok <- TRUE
  for (nm in names(EXPECT)) {
    d <- cp[cp$antigen_iso == nm, ]
    got <- c(
      y0 = median(d$y0), y1 = median(d$y1), t1 = median(d$t1),
      r = median(d$r), kdec = median(d$alpha * d$y1^(d$r - 1))
    )
    exp_ <- EXPECT[[nm]]
    ratio <- got / exp_
    say(sprintf(
      "  %-9s y0 %6.1f  y1 %7.1f  t1 %5.1f d  r %5.2f  k %7.4f/d",
      nm, got["y0"], got["y1"], got["t1"], got["r"], got["kdec"]
    ))
    say(sprintf(
      "  %-9s expected %5.1f   %7.1f    %5.1f      %5.2f    %7.4f",
      "", exp_["y0"], exp_["y1"], exp_["t1"], exp_["r"], exp_["kdec"]
    ))
    bad <- names(ratio)[ratio < 0.4 | ratio > 2.5]
    if (length(bad)) {
      ok <- FALSE
      say(sprintf(
        "  *** %s off by more than a factor of 2.5: %s ***",
        nm, paste(bad, collapse = ", ")
      ))
    }
  }
  ok
}

report_correlation <- function(cp, label) {
  g <- cp[cp$antigen_iso == "HlyE_IgG", ]
  g <- g[order(g$iter), ]
  a <- cp[cp$antigen_iso == "HlyE_IgA", ]
  a <- a[order(a$iter), ]
  r <- c(
    y0 = cor(log(g$y0), log(a$y0)),
    y1 = cor(log(g$y1), log(a$y1)),
    t1 = cor(log(g$t1), log(a$t1)),
    alpha = cor(log(g$alpha), log(a$alpha)),
    r = cor(log(g$r - 1), log(a$r - 1))
  )
  say(sprintf("  cross-isotype correlation in the bank:"))
  say(sprintf(
    "    y0 %+.3f   y1 %+.3f   t1 %+.3f   alpha %+.3f   r %+.3f",
    r["y0"], r["y1"], r["t1"], r["alpha"], r["r"]
  ))
  invisible(r)
}

# =====================================================================
res <- list()
for (arm in names(ARMS)) {
  p <- path.expand(file.path(CH2, sprintf("fit_ch2_461_%s.rds", arm)))
  hr()
  say(sprintf("arm %s  (%s)   %s", arm, ARMS[[arm]], p))
  if (!file.exists(p)) {
    say("  not found; skipping")
    next
  }

  fit <- readRDS(p)
  cp <- to_curve_params(fit, N_DRAWS, LP_DROP)

  say("")
  if (!check_layout(cp)) {
    say("")
    say("  The par_new mapping does not reproduce the population medians.")
    say("  Do not use this bank. Check the layout reported above against the")
    say("  parameters block of the Stan file before going further.")
    quit(status = 1)
  }
  say("")
  r <- report_correlation(cp, arm)

  fn <- file.path(OUT, sprintf("curve_params_arm_%s.rds", arm))
  attr(cp, "arm") <- arm
  attr(cp, "role") <- ARMS[[arm]]
  attr(cp, "lp_drop") <- LP_DROP
  attr(cp, "source") <- sprintf(
    "%s, seed %s, %d chains, %d kept",
    p, unique(fit$metadata()$seed),
    length(fit$metadata()$id),
    length(attr(cp, "chains_kept"))
  )
  saveRDS(cp, fn)
  say(sprintf(
    "  -> %s   %d draws x 2 isotypes   %.2f MB",
    fn, nrow(cp) / 2, file.size(fn) / 1e6
  ))
  res[[arm]] <- r
  rm(fit)
  invisible(gc(verbose = FALSE))
}

hr()
if (length(res) == 2) {
  say("The contrast the two banks provide:")
  say(sprintf("  %-8s %8s %8s %8s %8s %8s", "arm", "y0", "y1", "t1", "alpha", "r"))
  for (a in names(res)) {
    say(sprintf(
      "  %-8s %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f",
      a, res[[a]]["y0"], res[[a]]["y1"], res[[a]]["t1"],
      res[[a]]["alpha"], res[[a]]["r"]
    ))
  }
  say("")
  # alpha is a derived quantity and the least well identified, so it is the
  # one where between-chain variation in the thinned bank shows up. Judge the
  # bank on the four sampled parameters and flag alpha separately.
  s_main <- res$S[c("y0", "y1", "t1", "r")]
  if (max(abs(s_main)) > 0.10) {
    say("  *** arm S should be near zero on y0, y1, t1 and r; it is not. ***")
  } else if (abs(res$S["alpha"]) > 0.10) {
    say(sprintf(
      "  arm S is near zero on the four sampled parameters (max %.3f);\n  alpha is %+.3f. Run R/06_check_bank_correlation.R before using this.",
      max(abs(s_main)), res$S["alpha"]
    ))
  } else if (min(res$T[c("y1", "alpha", "r")]) < 0.3) {
    say("  *** arm T should show the chapter 2 correlations; it does not. ***")
  } else {
    say("  arm T carries the association, arm S does not. The banks are usable.")
  }
}
