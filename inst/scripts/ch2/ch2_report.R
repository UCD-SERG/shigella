#!/usr/bin/env Rscript
# =====================================================================
# ch2_report.R -- score every simulation seed and the real-data fits under all
# three schemes, and print the table that goes in front of Ezra.
#
#   all chains   no exclusion
#   lp filter    drop chains whose median log posterior is more than LP_DROP
#                below the best  (chain-picking, Yao/Vehtari/Gelman 2022)
#   stacking     weight chains by leave-one-participant-out predictive
#                performance     (the same paper's recommended alternative)
#
# FIXES IN THIS VERSION, in the order they were found:
#   1  the score cache (fit_seedNNN_scored.rds) is excluded from the fit glob.
#      Without this the cache is treated as a fit and scored again, giving
#      fit_seedNNN_scored_scored.rds. Applies to BOTH globs -- the scoring loop
#      and the threshold-sensitivity loop.
#   2  log_lik is read once per fit, not once per chain. Reading it inside the
#      chain loop pulled the whole iter x chain x obs array on every pass --
#      1.06 GB for a 40-chain fit, forty times over.
#   3  r_eff coming back NA for a single chain is caught and replaced with 1
#      rather than left to loo's warning path.
#   4  the summary frame names its columns explicitly; cv/se/wd/bi all carry
#      the names all/lp/stacking and collided when bound with t().
#   5  a seed can finish between the scoring loop's glob and the sensitivity
#      loop's, leaving a fit with no cache. That fit is skipped and the count
#      of seeds actually used is printed.
#
# WHY CACHING. Stacking needs PSIS-LOO for every chain of every fit, which is
# the only slow part -- roughly 20 s for a 10-chain simulation fit and several
# minutes for a 40-chain real-data fit. Everything downstream is arithmetic on
# a few hundred kilobytes, so each fit is scored once and re-runs touch only
# seeds that appeared since.
#
# USAGE
#   Rscript ch2_report.R                    simulations only
#   Rscript ch2_report.R --real             simulations plus arm T, T2, T3
#   Rscript ch2_report.R --refresh          ignore the cache and recompute
#   LP_DROP=50 Rscript ch2_report.R         a different threshold
#   Rscript ch2_report.R --quiet            summary table only
#
#   nice -n 19 Rscript ch2_report.R --real  when the simulation grid is running
# =====================================================================

suppressPackageStartupMessages({library(posterior); library(loo)})

args    <- commandArgs(trailingOnly = TRUE)
REAL    <- "--real"    %in% args
REFRESH <- "--refresh" %in% args
QUIET   <- "--quiet"   %in% args
LP_DROP <- as.numeric(Sys.getenv("LP_DROP", "100"))
SIM_DIRS <- strsplit(Sys.getenv("SIM_DIRS",
              "sim_out_null_prod,sim_out_alt_prod"), ",")[[1]]
REAL_FITS <- strsplit(Sys.getenv("REAL_FITS",
              "fit_ch2_461_T.rds,fit_ch2_461_T2.rds,fit_ch2_461_T3.rds"), ",")[[1]]
SUBSET   <- Sys.getenv("SUBSET", "sees_subset_v2.rds")
NM5 <- c("baseline", "peak", "peak-time", "decay-rate", "decay-shape")
S_RESAMPLE <- 20000L

say <- function(...) if (!QUIET) cat(..., "\n", sep = "")
hr  <- function(ch = "-") if (!QUIET) cat(strrep(ch, 78), "\n")

# =====================================================================
# scoring one fit  ->  a small cached object
# =====================================================================
score_fit <- function(path, id_map = NULL) {
  cache <- sub("\\.rds$", "_scored.rds", path)
  if (!REFRESH && file.exists(cache) &&
      file.mtime(cache) > file.mtime(path)) return(readRDS(cache))

  say("  scoring ", basename(path), " ...")
  f  <- readRDS(path)
  lp <- apply(f$draws("lp__", format = "draws_array"), 2, median)
  K  <- length(lp)
  cc <- as.matrix(f$draws("cross_corr", format = "draws_matrix"))
  nit <- nrow(cc) / K
  # iter x chain x 5, so every scheme below is a column subset or a resample
  ccarr <- array(NA_real_, c(nit, K, 5))
  for (j in 1:5) ccarr[, , j] <- matrix(cc[, j], nrow = nit, ncol = K)

  # ---- stacking weights ----------------------------------------------------
  w <- rep(1/K, K); kbad <- rep(NA_real_, K); loo_unit <- "none"
  if ("log_lik" %in% f$metadata()$stan_variables) {
    # Read log_lik ONCE. Calling f$draws() inside the chain loop would pull the
    # whole iter x chain x obs array on every pass -- 1.06 GB for a 40-chain fit
    # with 3,318 observations, forty times over. Reading once and subsetting a
    # matrix costs 26 MB per chain.
    LL <- f$draws("log_lik", format = "draws_array")
    nobs <- dim(LL)[3]
    LLm <- matrix(as.numeric(LL), nrow = nit * K)   # (iter within chain) x obs
    rm(LL); invisible(gc(verbose = FALSE))

    # Participant-level, not observation-level. With 3.6 observations per
    # participant per isotype, leaving out one observation leaves the others to
    # pin down that participant's curve, so the chains look alike and the
    # weights carry no information. Leaving out a participant is the unit the
    # model comparison uses and the one that separates the modes.
    grp <- if (!is.null(id_map) && length(id_map) == nobs) id_map else NULL
    loo_unit <- if (is.null(grp)) "observation" else "participant"

    loos <- vector("list", K)
    for (k in seq_len(K)) {
      m <- LLm[((k - 1) * nit + 1):(k * nit), , drop = FALSE]
      if (!is.null(grp)) m <- t(rowsum(t(m), group = grp))  # sum within participant
      a <- array(m, c(nit, 1L, ncol(m)))
      # r_eff often comes back NA for a single chain; loo then substitutes 1,
      # which is the conservative choice and what we want here.
      re <- suppressWarnings(try(relative_eff(exp(a)), silent = TRUE))
      if (inherits(re, "try-error") || all(is.na(re))) re <- rep(1, ncol(m))
      loos[[k]] <- suppressWarnings(loo(a, r_eff = re))
      kbad[k] <- mean(loos[[k]]$diagnostics$pareto_k > 0.7)
      if (!QUIET && K > 20 && k %% 10 == 0) cat("    chain", k, "\n")
    }
    rm(LLm); invisible(gc(verbose = FALSE))
    w <- as.numeric(suppressWarnings(
      loo_model_weights(loos, method = "stacking",
                        optim_control = list(reltol = 1e-10))))
  } else {
    say("    no log_lik; stacking unavailable, weights left uniform")
    loo_unit <- "unavailable"
  }

  dg <- try(f$diagnostic_summary(), silent = TRUE)
  out <- list(
    path = path, K = K, nit = nit, lp = lp, cc = ccarr, w = w,
    pareto_bad = kbad, loo_unit = loo_unit,
    baseline_by_chain = apply(ccarr[, , 1, drop = FALSE], 2, median),
    n_div = if (inherits(dg, "try-error")) NA else sum(dg$num_divergent),
    ebfmi_low = if (inherits(dg, "try-error")) NA else sum(dg$ebfmi < 0.3))
  saveRDS(out, cache)
  out
}

# =====================================================================
# three sets of draws from one scored fit
# =====================================================================
draws_by_scheme <- function(s, lp_drop = LP_DROP, seed = 1) {
  keep <- which(s$lp > max(s$lp) - lp_drop)
  set.seed(seed)
  pick <- sample.int(s$K, S_RESAMPLE, replace = TRUE, prob = s$w)
  iter <- sample.int(s$nit, S_RESAMPLE, replace = TRUE)
  lapply(1:5, function(j) list(
    all      = as.vector(s$cc[, , j]),
    lp       = as.vector(s$cc[, keep, j, drop = FALSE]),
    stacking = s$cc[cbind(iter, pick, j)]))
}

summ <- function(x) c(med = median(x), q05 = unname(quantile(x, .05)),
                      q95 = unname(quantile(x, .95)))

# =====================================================================
# simulations
# =====================================================================
grids <- list()
for (d in SIM_DIRS) {
  fits <- sort(Sys.glob(file.path(d, "fit_seed*.rds")))
  # The score cache lives beside the fit as fit_seedNNN_scored.rds and the glob
  # above matches it too. Without this line the cache gets treated as a fit and
  # scored again, producing fit_seedNNN_scored_scored.rds.
  fits <- fits[!grepl("_scored\\.rds$", fits)]
  if (!length(fits)) { say("\n", d, ": no saved fits"); next }
  res <- sort(Sys.glob(file.path(d, "res_seed*.rds")))
  if (!length(res)) { say("\n", d, ": no res_seed to read truth from"); next }
  r1 <- readRDS(res[1])
  truth <- r1$truth[grepl("cross_corr", r1$variable)]

  # participant map: the simulation reuses the real visit schedules, so the
  # number of observations varies by seed and the real subset cannot supply it.
  id_map <- NULL

  hr("=")
  say(d, "   ", length(fits), " fits   truth ",
      paste(sprintf("%.3f", truth), collapse = " "))
  hr("=")

  rows <- list()
  for (p in fits) {
    s  <- score_fit(p, id_map)
    dw <- draws_by_scheme(s)
    sid <- as.integer(gsub("\\D", "", basename(p)))
    keep_n <- sum(s$lp > max(s$lp) - LP_DROP)
    for (j in 1:5) {
      for (sc in c("all", "lp", "stacking")) {
        q <- summ(dw[[j]][[sc]])
        rows[[length(rows) + 1]] <- data.frame(
          dir = d, seed = sid, j = j, param = NM5[j], scheme = sc,
          truth = truth[j], med = q["med"], q05 = q["q05"], q95 = q["q95"],
          covered = truth[j] >= q["q05"] && truth[j] <= q["q95"],
          width = q["q95"] - q["q05"], keep_n = keep_n, K = s$K,
          w_eff = 1/sum(s$w^2), loo_unit = s$loo_unit,
          n_div = s$n_div, row.names = NULL)
      }
    }
    if (!QUIET) {
      cv <- sapply(c("all","lp","stacking"), function(sc)
        sum(sapply(1:5, function(j) {
          q <- summ(dw[[j]][[sc]]); truth[j] >= q["q05"] && truth[j] <= q["q95"] })))
      cat(sprintf("    seed %03d   keep %2d/%2d   w_eff %4.1f   포함 %d / %d / %d\n",
                  sid, keep_n, s$K, 1/sum(s$w^2), cv[1], cv[2], cv[3]))
    }
  }
  grids[[d]] <- do.call(rbind, rows)
}

# =====================================================================
# summary
# =====================================================================
hr("=")
cat(sprintf("\nCOVERAGE   nominal 0.90   LP_DROP = %g\n\n", LP_DROP))
cat(sprintf("  %-22s %5s %5s %10s %10s %10s\n",
            "grid", "seeds", "int", "all", "lp filter", "stacking"))
final <- list()
for (d in names(grids)) {
  g <- grids[[d]]
  ns <- length(unique(g$seed)); ni <- ns * 5
  cv <- sapply(c("all","lp","stacking"), function(sc) mean(g$covered[g$scheme == sc]))
  se <- sqrt(cv * (1 - cv) / ni)
  wd <- sapply(c("all","lp","stacking"), function(sc) mean(g$width[g$scheme == sc]))
  bi <- sapply(c("all","lp","stacking"), function(sc)
               mean(g$med[g$scheme == sc] - g$truth[g$scheme == sc]))
  cat(sprintf("  %-22s %5d %5d %10.3f %10.3f %10.3f\n", d, ns, ni,
              cv["all"], cv["lp"], cv["stacking"]))
  cat(sprintf("  %-22s %5s %5s %10s %10s %10s\n", "  MCSE", "", "",
              sprintf("(%.3f)", se["all"]), sprintf("(%.3f)", se["lp"]),
              sprintf("(%.3f)", se["stacking"])))
  cat(sprintf("  %-22s %5s %5s %10.3f %10.3f %10.3f\n", "  mean width", "", "",
              wd["all"], wd["lp"], wd["stacking"]))
  cat(sprintf("  %-22s %5s %5s %+10.3f %+10.3f %+10.3f\n\n", "  mean bias", "", "",
              bi["all"], bi["lp"], bi["stacking"]))
  final[[d]] <- data.frame(dir = d, n_seed = ns, n_int = ni,
    cov_all = cv["all"], cov_lp = cv["lp"], cov_stack = cv["stacking"],
    se_all  = se["all"], se_lp  = se["lp"], se_stack  = se["stacking"],
    wid_all = wd["all"], wid_lp = wd["lp"], wid_stack = wd["stacking"],
    bias_all = bi["all"], bias_lp = bi["lp"], bias_stack = bi["stacking"],
    row.names = NULL)
}

# the pre-specified rule, applied
cat("  pre-specified rule: within 2 MCSE of 0.90 counts as nominal\n\n")
for (d in names(grids)) {
  g <- grids[[d]]; ni <- length(unique(g$seed)) * 5
  for (sc in c("all","lp","stacking")) {
    p <- mean(g$covered[g$scheme == sc]); se <- sqrt(p*(1-p)/ni)
    cat(sprintf("    %-22s %-9s %.3f   |diff| %.3f   2MCSE %.3f   %s\n",
                d, sc, p, abs(p - .90), 2*se,
                if (abs(p - .90) <= 2*se) "pass" else "*** outside ***"))
  }
}

# =====================================================================
# threshold sensitivity -- free, the lp values are cached
# =====================================================================
cat(sprintf("\n\nTHRESHOLD SENSITIVITY   (recomputed from cache, no refitting)\n\n"))
cat(sprintf("  %-22s %10s %10s %10s %10s\n", "grid", "LP_DROP 50", "100", "200", "stacking"))
for (d in SIM_DIRS) {
  fits <- sort(Sys.glob(file.path(d, "fit_seed*.rds")))
  # The score cache lives beside the fit as fit_seedNNN_scored.rds and the glob
  # above matches it too. Without this line the cache gets treated as a fit and
  # scored again, producing fit_seedNNN_scored_scored.rds.
  fits <- fits[!grepl("_scored\\.rds$", fits)]
  if (!length(fits)) next
  res <- sort(Sys.glob(file.path(d, "res_seed*.rds"))); if (!length(res)) next
  truth <- readRDS(res[1]); truth <- truth$truth[grepl("cross_corr", truth$variable)]
  cov_at <- function(L) {
    h <- 0; n <- 0
    for (p in fits) {
      # A seed can finish between the scoring loop's glob and this one, in which
      # case the fit exists but its cache does not. Skip it here rather than
      # abort; the count of seeds used is printed alongside.
      cache <- sub("\\.rds$", "_scored.rds", p)
      if (!file.exists(cache)) next
      s <- readRDS(cache)
      keep <- which(s$lp > max(s$lp) - L)
      for (j in 1:5) {
        q <- summ(as.vector(s$cc[, keep, j, drop = FALSE]))
        h <- h + (truth[j] >= q["q05"] && truth[j] <= q["q95"]); n <- n + 1
      }
    }
    if (n == 0) return(NA_real_)
    h/n
  }
  st <- mean(grids[[d]]$covered[grids[[d]]$scheme == "stacking"])
  n_used <- sum(file.exists(sub("\\.rds$", "_scored.rds", fits)))
  cat(sprintf("  %-22s %10.3f %10.3f %10.3f %10.3f   (%d of %d scored)\n", d,
              cov_at(50), cov_at(100), cov_at(200), st, n_used, length(fits)))
}

# =====================================================================
# real-data fits
# =====================================================================
if (REAL) {
  id_map <- NULL
  if (file.exists(SUBSET)) {
    d0 <- readRDS(SUBSET)
    id_map <- as.integer(factor(d0[[1]]))
    cat(sprintf("\n\nREAL DATA   participant map from %s: %d participants, %d observations\n",
                SUBSET, length(unique(id_map)), length(id_map)))
  } else {
    cat("\n\nREAL DATA   subset not found; stacking will use observation-level LOO\n")
  }
  for (p in REAL_FITS) {
    if (!file.exists(p)) { cat(sprintf("\n  %s not found\n", p)); next }
    s  <- score_fit(p, id_map)
    dw <- draws_by_scheme(s)
    keep <- which(s$lp > max(s$lp) - LP_DROP)
    cat(sprintf("\n  === %s ===  keep %d/%d   w_eff %.1f   LOO unit %s\n",
                basename(p), length(keep), s$K, 1/sum(s$w^2), s$loo_unit))
    minor <- which(s$baseline_by_chain < 0)
    if (length(minor))
      cat(sprintf("      minor-mode chains: %s   lp gap %s   stacking weight %.4f\n",
                  paste(minor, collapse = ","),
                  paste(sprintf("%.0f", s$lp[minor] - max(s$lp)), collapse = ","),
                  sum(s$w[minor])))
    cat(sprintf("      %-12s %24s %24s %24s\n", "", "all chains", "lp filter", "stacking"))
    for (j in 1:5) {
      f3 <- sapply(c("all","lp","stacking"), function(sc) {
        q <- summ(dw[[j]][[sc]])
        sprintf("%+.3f [%+.3f,%+.3f]", q["med"], q["q05"], q["q95"]) })
      cat(sprintf("      %-12s %24s %24s %24s%s\n", NM5[j], f3[1], f3[2], f3[3],
                  if (all(sapply(c("all","lp","stacking"), function(sc) {
                        q <- summ(dw[[j]][[sc]]); q["q05"] > 0 }))) "  excl 0" else ""))
    }
    meds <- sapply(c("all","lp","stacking"), function(sc)
      sapply(1:5, function(j) summ(dw[[j]][[sc]])["med"]))
    cat(sprintf("      largest median difference across schemes: %.3f\n",
                max(apply(meds, 1, function(r) diff(range(r))))))
  }
}

saveRDS(list(grids = grids, summary = do.call(rbind, final),
             lp_drop = LP_DROP, when = Sys.time()), "ch2_report.rds")
cat("\n\nsaved -> ch2_report.rds\n")
cat("re-run after new seeds finish; cached fits are not rescored.\n")
