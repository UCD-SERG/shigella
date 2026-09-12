#!/usr/bin/env Rscript
# =====================================================================
# 04_run_2x2.R -- estimate seroincidence four ways and compare
#
# THE DESIGN. Two banks from Chapter 2 (correlated and independent) crossed
# with three likelihood treatments, plus one cheap correction:
#
#   A   composite         product of separately integrated marginals  (current)
#   A'  A + sandwich      same estimate, cluster-robust SE with the person
#                         as the cluster -- this is the Godambe correction for
#                         a composite likelihood, and it is already implemented
#   C   joint             one integral over a shared t, draws paired
#
# FOUR CONTRASTS, and two of them are checks rather than results:
#
#   (C, S)  - (A, S)      what sharing t is worth
#   (C, T)  - (C, S)      what the Chapter 2 correlation is worth
#   (A, T)  - (A, S)      PREDICTED NEAR ZERO -- the composite form cannot use
#                         the correlation, so swapping banks should barely move it
#   (A', ·) - (A, ·)      how much of the interval the sandwich alone recovers
#
# Recording the third prediction before running makes it a test of the argument
# rather than a description of the output.
#
#   Rscript R/04_run_2x2.R
#   ARMS=highE Rscript R/04_run_2x2.R          one endemicity stratum
#   N_TAU=800 Rscript R/04_run_2x2.R
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(serocalculator)
})
source("R/03_log_likelihood_joint.R")

DATA <- Sys.getenv("DATA_DIR", "data")
OUT <- Sys.getenv("OUT_DIR", "out")
N_TAU <- as.integer(Sys.getenv("N_TAU", "400"))
ARMS <- strsplit(Sys.getenv("ARMS", "highE,lowE"), ",")[[1]]
# Restricting to one country is the quickest way to see the effect of a change
# in isolation -- the three differ enormously in run time because their age
# distributions differ, and Nepal is roughly five times slower than Bangladesh.
CTRYS <- Sys.getenv("COUNTRIES", "")
MAX_N <- as.integer(Sys.getenv("MAX_N", "0")) # 0 = all; set small to smoke-test
# The interval endpoints are read off the lambda grid, so its spacing is the
# resolution of every width comparison reported below. At 31 points from 0.005
# to 1 the ratio between neighbours is 1.19, so a "+23% width" can be a single
# grid step. 121 points brings that to 4.5%.
N_LAM <- as.integer(Sys.getenv("N_LAM", "121"))
LAM_LO <- as.numeric(Sys.getenv("LAM_LO", "0.005"))
LAM_HI <- as.numeric(Sys.getenv("LAM_HI", "1.0"))

# --- which banks and which likelihood forms to run --------------------
# defaults reproduce the previous behaviour exactly
BANKS <- strsplit(Sys.getenv("BANKS", "S,T,Tshuf"), ",")[[1]]
FORMS <- strsplit(Sys.getenv("FORMS", "composite,joint"), ",")[[1]]
ISOS <- c("HlyE_IgG", "HlyE_IgA")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(..., "\n", sep = "")
hr <- function(ch = "-") cat(strrep(ch, 74), "\n")

# antigen_iso must be a plain character everywhere; two factors with different
# level sets raise an error rather than comparing as FALSE.
chr <- function(d) {
  d <- as.data.frame(d)
  d$antigen_iso <- as.character(d$antigen_iso)
  d
}

xs <- chr(readRDS(file.path(DATA, "sees_xs_subset.rds")))
noise <- chr(serocalculator::example_noise_params_sees)
# Tshuf is arm T with the IgA rows reordered: identical marginals, pairing
# destroyed. arm S cannot serve as the control because it is a different fit
# and its IgA decay rate differs from arm T's by 35%, so a contrast against it
# mixes the correlation with everything else that differs.
banks <- list(
  S = chr(readRDS(file.path(DATA, "curve_params_arm_S.rds"))),
  T = chr(readRDS(file.path(DATA, "curve_params_arm_T.rds"))),
  Tshuf = chr(readRDS(file.path(DATA, "curve_params_arm_Tshuf.rds")))
)

# --- keep only the requested banks ------------------------------------
if (!all(BANKS %in% names(banks))) {
  stop("unknown bank in BANKS: ",
    paste(setdiff(BANKS, names(banks)), collapse = ", "),
    call. = FALSE
  )
}
banks <- banks[BANKS]

# --- sensitivity levers: LOD_STD and AGE_MAX --------------------------
# no-ops unless the environment variable is set.  Both the data threshold
# and noise$y.low are changed together, and the helper verifies it.
source(file.path("R", "09_lod_age_sensitivity.R"))
.sens <- apply_sensitivity(xs, noise)
xs <- .sens$data
noise <- .sens$noise

hr("=")
say("BANKS")
hr("=")
for (b in names(banks)) {
  cp <- banks[[b]]
  g <- cp[cp$antigen_iso == "HlyE_IgG", ]
  g <- g[order(g$iter), ]
  a <- cp[cp$antigen_iso == "HlyE_IgA", ]
  a <- a[order(a$iter), ]
  say(sprintf(
    "  arm %s  %d draws   cor(log y1) %+.3f   cor(log alpha) %+.3f",
    b, nrow(g), cor(log(g$y1), log(a$y1)), cor(log(g$alpha), log(a$alpha))
  ))
}
say("\n  arm T should carry the association and arm S should not.")
say(sprintf(
  "\n  lambda grid: %d points over %.3f to %.2f, neighbours %.1f%% apart",
  N_LAM, LAM_LO, LAM_HI,
  100 * ((LAM_HI / LAM_LO)^(1 / (N_LAM - 1)) - 1)
))
say("  interval widths are resolved no finer than that")
if (nzchar(CTRYS)) say("\n  restricted to: ", CTRYS)

# =====================================================================
# each (arm, country) is analysed with its own detection limits
# =====================================================================
profile_ci <- function(ll_fun, lams, ...) {
  v <- vapply(lams, function(l) ll_fun(l, ...), numeric(1))
  i <- which.max(v)
  lo <- suppressWarnings(max(lams[seq_len(i)][v[seq_len(i)] < max(v) - 1.92]))
  hi <- suppressWarnings(min(lams[i:length(lams)][v[i:length(lams)] < max(v) - 1.92]))
  c(est = lams[i], lo = lo, hi = hi, width = hi - lo, max_ll = max(v))
}

LAMS <- exp(seq(log(LAM_LO), log(LAM_HI), length.out = N_LAM))
res <- list()

ctry_list <- function(arm) {
  a <- sort(unique(xs$Country[xs$arm == arm]))
  if (nzchar(CTRYS)) a <- intersect(a, strsplit(CTRYS, ",")[[1]])
  a
}

for (arm in ARMS) {
  for (ctry in ctry_list(arm)) {
    d <- xs |> filter(arm == !!arm, Country == ctry)
    if (MAX_N > 0) d <- d |> filter(id %in% head(unique(d$id), MAX_N))
    np <- noise |> filter(Country == ctry, antigen_iso %in% ISOS)
    if (nrow(np) != 2) {
      say("  skip ", arm, "/", ctry, ": noise rows ", nrow(np))
      next
    }

    pd <- as_pop_data(d,
      age = "age", value = "value", id = "id",
      antigen_isos = ISOS
    )

    hr("=")
    say(sprintf("%s / %s   %d people", arm, ctry, n_distinct(d$id)))
    hr("=")

    for (b in names(banks)) {
      cp <- as_sr_params(banks[[b]] |> filter(antigen_iso %in% ISOS))
      fits <- list()

      if ("composite" %in% FORMS) {
        t0 <- Sys.time()
        A <- profile_ci(function(l, ...) {
          serocalculator::log_likelihood(
            l, pd, cp, np,
            antigen_isos = ISOS
          )
        }, LAMS)
        say(sprintf(
          "  bank %s   composite  %.4f [%.4f, %.4f]  width %.4f   (%.1f min)",
          b, A["est"], A["lo"], A["hi"], A["width"],
          as.numeric(difftime(Sys.time(), t0, units = "mins"))
        ))
        fits$composite <- A
      }

      if ("joint" %in% FORMS) {
        t0 <- Sys.time()
        C <- profile_ci(function(l, ...) {
          log_likelihood_joint(
            l, d, banks[[b]], np,
            antigen_isos = ISOS, n_tau = N_TAU
          )
        }, LAMS)
        say(sprintf(
          "           joint      %.4f [%.4f, %.4f]  width %.4f   (%.1f min)",
          C["est"], C["lo"], C["hi"], C["width"],
          as.numeric(difftime(Sys.time(), t0, units = "mins"))
        ))
        fits$joint <- C
      }

      if (length(fits)) {
        res[[length(res) + 1]] <- data.frame(
          arm = arm, Country = ctry, bank = b, n = n_distinct(d$id),
          form = names(fits),
          do.call(rbind, fits), se_naive = NA_real_, se_sandwich = NA_real_,
          row.names = NULL
        )
      }
    }

    # ---- the sandwich, if the package exposes it ----------------------------
    ok <- tryCatch(
      {
        if (!"T" %in% names(banks)) stop("bank T not selected; sandwich skipped")
        e1 <- est_seroincidence(
          pop_data = pd, sr_params = as_sr_params(banks$T),
          noise_params = np, antigen_isos = ISOS
        )
        e2 <- est_seroincidence(
          pop_data = pd, sr_params = as_sr_params(banks$T),
          noise_params = np, antigen_isos = ISOS,
          cluster_var = "id"
        )
        s1 <- summary(e1)
        s2 <- summary(e2)
        say(sprintf(
          "\n  sandwich check   naive SE %.5f   cluster=person SE %.5f   ratio %.2f",
          s1$SE, s2$SE, s2$SE / s1$SE
        ))
        say("  a ratio above 1 confirms the composite likelihood understates uncertainty")
        # Only the SE ratio is comparable. The composite width above comes from a
        # profile on the lambda grid; the sandwich width would be a Wald interval
        # from a standard error. Putting them side by side compares two different
        # constructions, so the width is recorded but not contrasted.
        res[[length(res) + 1]] <- data.frame(
          arm = arm, Country = ctry, bank = "T", n = n_distinct(d$id),
          form = "composite+sandwich", est = s2$incidence.rate,
          lo = s2$CI.lwr, hi = s2$CI.upr, width = s2$CI.upr - s2$CI.lwr,
          max_ll = NA_real_, se_naive = s1$SE, se_sandwich = s2$SE,
          row.names = NULL
        )
        TRUE
      },
      error = function(e) {
        say("\n  sandwich unavailable: ", conditionMessage(e))
        FALSE
      }
    )
  }
}

out <- do.call(rbind, res)
# One file per (arm, country) so the six cells can run as separate processes
# and be combined afterwards.
OUT_FILE <- Sys.getenv("OUT_FILE", "ch3_2x2.rds")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file.path(OUT, OUT_FILE))

hr("=")
say("SUMMARY")
hr("=")
print(as.data.frame(out), row.names = FALSE, digits = 4)

hr("=")
say("THE FOUR CONTRASTS")
hr("=")
w <- function(a, c_, f) {
  x <- out |> filter(arm == a, Country == c_, bank == "S" | bank == "T")
  x
}
for (a in unique(out$arm)) {
  for (ct in unique(out$Country[out$arm == a])) {
    x <- out |> filter(arm == a, Country == ct)
    g <- function(bk, fm) {
      v <- x$width[x$bank == bk & x$form == fm]
      if (length(v)) v[1] else NA
    }
    pc <- function(p, q) if (is.na(p) || is.na(q)) NA else round(100 * (p - q) / q, 1)
    say(sprintf("  %s / %s", a, ct))
    say(sprintf(
      "    sharing t      (joint,S) vs (comp,S)   width %+.1f%%",
      pc(g("S", "joint"), g("S", "composite"))
    ))
    say(sprintf(
      "    corr + margins (joint,T) vs (joint,S)     width %+.1f%%",
      pc(g("T", "joint"), g("S", "joint"))
    ))
    say(sprintf(
      "    correlation    (joint,T) vs (joint,Tshuf) width %+.1f%%",
      pc(g("T", "joint"), g("Tshuf", "joint"))
    ))
    say(sprintf(
      "    PREDICTED = 0  (comp,T)  vs (comp,Tshuf)  width %+.1f%%",
      pc(g("T", "composite"), g("Tshuf", "composite"))
    ))
    say("                   these two banks have identical marginals, so the")
    say("                   composite form must give the SAME number, not a close one")
    sr <- x$se_sandwich[x$form == "composite+sandwich"]
    sn <- x$se_naive[x$form == "composite+sandwich"]
    if (length(sr) && is.finite(sr)) {
      say(sprintf(
        "    sandwich       SE ratio %.2f   (%.5f -> %.5f)",
        sr / sn, sn, sr
      ))
      say("                   above 1 means the composite understates uncertainty")
    }
    cat("\n")
  }
}

hr("=")
say("PREDICTED: the point estimate barely moves and the joint interval is WIDER,")
say("because the composite form counts one specimen as two independent people.")
say("A narrower joint interval would contradict the argument and needs investigating.")
say("")
say("saved -> ", file.path(OUT, OUT_FILE))
