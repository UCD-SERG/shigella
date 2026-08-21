#!/usr/bin/env Rscript
# =====================================================================
# run_ch2_461.R -- one Chapter 2 arm on sees_subset_v2.rds (461 people).
#
# Everything is arm L's setting except the chain count, which goes up because
# there are cores to spend. Arm L was chosen over K because the only thing L
# changed was `init`, and an initial value cannot move a posterior; K narrowed
# c_prior_sd from 1 to 0.5, which can. K stays as the sensitivity arm.
#
# Settings and why each one is what it is:
#   warmup 3000        arm M used 5000 and came out WORSE (R-hat 1.170 vs L's
#                      1.086). Warmup is not the binding constraint here.
#   sampling 1000      arm L reached min ESS 109 on 14 chains. At 40 chains the
#                      same 1000 draws give roughly 40/14 x 109 = 310, which is
#                      far more than the elpd comparison needs. Raising it would
#                      buy ESS nobody is short of and cost hours and disk.
#   max_treedepth 12   arm H used 13: saturation fell from 100% to 75%, R-hat
#                      got worse (1.372 vs E's 1.281), and it took twice as long
#                      (27.9 h vs 14.3 h). Saturation is not what limits R-hat.
#   adapt_delta 0.9    divergences across every v4 arm were 0.01-0.02%. Nothing
#                      to buy by raising it, and lowering it risks the one
#                      diagnostic that has been clean throughout.
#   init 0.1           arm L. Cannot change the posterior.
#   c_prior_sd 1       untouched prior in the main arm, so nobody can ask whether
#                      the prior was narrowed until the answer came out right.
#
# ESTIMATE_C is the ONLY thing that differs between the two arms. That is what
# makes the elpd difference attributable to C.
#
# Environment variables (all have defaults):
#   ARM          tag used in filenames                     default R
#   ESTIMATE_C   1 = Chapter 2, 0 = null, 2 = baseline off default 1
#   DATA         rds holding the case data                 default sees_subset_v2.rds
#   MODEL        stan file                                 default model_ch2_v5.stan
#                arm T (the correlated bank) was fitted with v5, which removes the
#                (-1,1) cap on rho_c; arm S was fitted with v4, where estimate_c = 0
#                sets c = 0 without reference to rho_c, so the cap is inoperative.
#   NCHAIN NADAPT NITER ADAPT_DELTA MAX_TD INIT C_PRIOR_SD SEED
#   FUNCS        stan_ch2_functions.R
#   EXEC_DIR     compile dir; MUST differ between concurrently launched arms
# =====================================================================

ge <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
gn <- function(k, d) as.numeric(ge(k, as.character(d)))

ARM        <- ge("ARM", "R")
ESTIMATE_C <- as.integer(gn("ESTIMATE_C", 1))
DATA       <- ge("DATA",  "sees_subset_v2.rds")
MODEL      <- ge("MODEL", "model_ch2_v5.stan")
FUNCS      <- ge("FUNCS", "stan_ch2_functions.R")
NCHAIN     <- as.integer(gn("NCHAIN", 40))
NADAPT     <- as.integer(gn("NADAPT", 3000))
NITER      <- as.integer(gn("NITER",  1000))
ADAPT_D    <- gn("ADAPT_DELTA", 0.9)
MAX_TD     <- as.integer(gn("MAX_TD", 12))
INIT       <- gn("INIT", 0.1)
C_SD       <- gn("C_PRIOR_SD", 1)
SEED       <- as.integer(gn("SEED", 1))

say <- function(...) { cat(sprintf("[ch2_%s] ", ARM), ..., "\n", sep = ""); flush.console() }

for (f in c(DATA, MODEL, FUNCS))
  if (!file.exists(f)) stop("missing: ", f, call. = FALSE)

# Each concurrently running arm needs its own compile directory. Without this
# the two arms race to copy and build the same executable and one of them picks
# up a half-written binary.
if (!nzchar(Sys.getenv("CMDSTAN_EXEC_DIR")))
  Sys.setenv(CMDSTAN_EXEC_DIR = file.path(tempdir(), paste0("ch2_exe_", ARM)))
say("exec dir ", Sys.getenv("CMDSTAN_EXEC_DIR"))

source(FUNCS)
dat <- readRDS(DATA)

say("data ", DATA, ": ", nrow(dat), " rows")

# ---- guard: resolve columns and normalise biomarker levels ------------------
# prep_ch2_standata() keeps only rows whose biomarker is in
# c("HlyE_IgG", "HlyE_IgA"), and it resolves four column names from attributes
# with nepal_sees fallbacks. If any of that misses, EVERY row is dropped and the
# failure surfaces far downstream as max() on an empty vector:
#     Warning: no non-missing arguments to max; returning -Inf
#     Error in matrix(0, nsubj, max_nsmpl): invalid 'ncol' value
# The subset stores the biomarkers with a space; the old runner renamed them and
# this one has to as well. Check the rest here rather than find out in an hour.
.cols <- .ch2_cols(dat)
say("columns resolved: id=", .cols$id, " time=", .cols$time,
    " value=", .cols$value, " biomarker=", .cols$bio)
.absent <- setdiff(unlist(.cols), names(dat))
if (length(.absent))
  stop("resolved column name(s) not present in the data: ",
       paste(.absent, collapse = ", "),
       "\n  available columns: ", paste(names(dat), collapse = ", "), call. = FALSE)

say("biomarkers as stored: ",
    paste(sort(unique(as.character(dat[[.cols$bio]]))), collapse = ", "))
dat[[.cols$bio]] <- gsub(" ", "_", as.character(dat[[.cols$bio]]))
say("biomarkers passed on: ",
    paste(sort(unique(dat[[.cols$bio]])), collapse = ", "))
if (!all(c("HlyE_IgG", "HlyE_IgA") %in% unique(dat[[.cols$bio]])))
  stop("after normalisation the data still lacks HlyE_IgG and/or HlyE_IgA",
       call. = FALSE)

# Dry-run the filter so an empty result stops here with a readable message.
.chk <- prep_ch2_standata(dat, min_visits = 2L)
say("subjects reaching the model: ", .chk$nsubj,
    " | max visits per subject: ", .chk$max_nsmpl)
# NOTE the slot order. prep_ch2_standata defaults to isos = c("HlyE_IgG",
# "HlyE_IgA"), which is what model_ch2_v4.stan requires, so slot 1 is IgG here.
# Chapter 1 goes through prep_data_stan(), which sorts the isotypes and gives
# slot 1 = IgA. mu_par[1, ] therefore means DIFFERENT isotypes in the two
# chapters -- do not line them up by index when comparing.
say("slot 1 = ", .chk$.isos[1], " | slot 2 = ", .chk$.isos[2],
    "   (Chapter 1 is the other way round: slot 1 = IgA)")
if (.chk$nsubj < 100)
  stop("only ", .chk$nsubj, " subjects survived the paired-visit filter",
       call. = FALSE)
rm(.chk)

say("model ", MODEL, " | estimate_c = ", ESTIMATE_C,
    if (ESTIMATE_C == 1) "  (Chapter 2, all five free)"
    else if (ESTIMATE_C == 0) "  (null: C = 0)" else "  (baseline held at zero)")
say(NCHAIN, " chains | warmup ", NADAPT, " | sampling ", NITER,
    " | adapt_delta ", ADAPT_D, " | max_td ", MAX_TD)
say("init ", INIT, " | c_prior_sd ", C_SD, " | seed ", SEED)

t0 <- Sys.time()
fit <- run_mod_stan_ch2(
  data            = dat,
  file_mod        = MODEL,
  estimate_c      = ESTIMATE_C,
  nchain          = NCHAIN,
  parallel_chains = NCHAIN,
  nadapt          = NADAPT,
  niter           = NITER,
  adapt_delta     = ADAPT_D,
  max_treedepth   = MAX_TD,
  c_prior_sd      = C_SD,
  init            = INIT,
  seed            = SEED,
  min_visits      = 2L,
  # Exposes log_prob/grad_log_prob so loo() can moment-match the observations
  # with Pareto k above 0.7 -- an earlier comparison had 17-19% of them over the
  # threshold, which is enough to make a raw elpd difference unreadable. The
  # call is fail-safe: if the methods build cannot link, it compiles normally
  # and LOO falls back to plain PSIS.
  compile_methods = TRUE
)
say("elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 2), " hours")
say("subjects kept: ", attr(fit, "n_subjects_kept"))
say("loo moment-matching available: ", isTRUE(attr(fit, "methods_ok")))

fit$save_object(sprintf("fit_ch2_461_%s.rds", ARM))
say("[saved] fit_ch2_461_", ARM, ".rds")

# ---- diagnostics ------------------------------------------------------------
dg  <- fit$diagnostic_summary()
tot <- NCHAIN * NITER
say("")
say("=== diagnostics ===")
say(sprintf("  divergent %d (%.3f%%) | treedepth %d (%.1f%%)",
            sum(dg$num_divergent), 100 * sum(dg$num_divergent) / tot,
            sum(dg$num_max_treedepth), 100 * sum(dg$num_max_treedepth) / tot))
say(sprintf("  E-BFMI: %d of %d chains below 0.3 (min %.3f)",
            sum(dg$ebfmi < 0.3), length(dg$ebfmi), min(dg$ebfmi)))

vars <- c("mu_par", "prec_logy", "sd_G", "sd_A", "cross_corr", "cross_cov")
ss   <- fit$summary(vars)
say(sprintf("  max R-hat %.3f | min bulk ESS %.0f",
            max(ss$rhat, na.rm = TRUE), min(ss$ess_bulk, na.rm = TRUE)))
saveRDS(as.data.frame(ss), sprintf("summ_ch2_461_%s.rds", ARM))

# ---- cross-biomarker correlations -------------------------------------------
if (ESTIMATE_C != 0) {
  cc <- fit$summary("cross_corr", "median",
                    q = ~quantile(.x, c(0.05, 0.95)), "rhat", "ess_bulk")
  lab <- c("baseline", "peak", "peak-time", "decay-rate", "decay-shape")
  say("")
  say("=== cross-biomarker correlations ===")
  cat(sprintf("\n  %-13s %8s   %-18s %8s %8s %7s\n",
              "parameter", "median", "90% interval", "R-hat", "ESS", "excl 0"))
  for (i in seq_len(nrow(cc))) {
    lo <- cc[[3]][i]; hi <- cc[[4]][i]
    cat(sprintf("  %-13s %+8.3f   [%+.3f, %+.3f] %8.3f %8.0f %7s\n",
                lab[i], cc$median[i], lo, hi, cc$rhat[i], cc$ess_bulk[i],
                if (lo > 0 || hi < 0) "yes" else "no"))
  }
  # Per-chain baseline correlation: this is what exposed the two-mode problem in
  # arms E, F and H (8 of 72 chains sat at about -0.24 while the rest sat at
  # +0.68). Arms K, L and M showed none in 42 chains. Watch it here too.
  d1 <- fit$draws("cross_corr[1]", format = "draws_matrix")
  nd <- nrow(d1) / NCHAIN
  bych <- tapply(as.numeric(d1), rep(seq_len(NCHAIN), each = nd), median)
  say("")
  say("=== baseline correlation by chain ===")
  cat("  ", paste(sprintf("%+.2f", bych), collapse = " "), "\n", sep = "")
  say(sprintf("  spread %.3f -- %s", diff(range(bych)),
              if (diff(range(bych)) < 0.1) "the chains agree"
              else "SOME CHAINS SIT IN THE OTHER MODE; do not read the interval yet"))
}

pn <- try(fit$summary("par_new"), silent = TRUE)
say("")
say("par_new present: ", if (inherits(pn, "try-error")) "NO" else
      paste0("yes, ", nrow(pn), " components"))
say("done")
