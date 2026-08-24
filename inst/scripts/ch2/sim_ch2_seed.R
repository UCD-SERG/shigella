#!/usr/bin/env Rscript
# =====================================================================
# sim_ch2_seed.R -- one simulation replicate for the Chapter 2 model.
#
# WHY THE TRUTH IS NOT THE FITTED CROSS-CORRELATION
# arm R's five cross-correlations sit against the ceiling that rho_c in (-1,1)
# imposes: every q95 is at 99-100% of 1/sqrt(1+g_j). Setting the simulation
# truth there would put the sampler on a boundary, which is what produced 45 of
# 46 global failures in the earlier 100-seed run. So everything except the
# cross-correlations is taken from the fit, and the cross-correlations are set
# to a common interior value (RHO_TRUE, default 0.4). Two things follow:
#   - the answer is not contaminated by the ceiling question that arm T is
#     currently deciding, so this pilot keeps its value whichever way arm T goes
#   - a simulation study wants a grid of true values anyway, not one point;
#     0.4 is the first grid point, not a substitute for the fitted value
#
# WHY THE VISIT TIMES ARE RESAMPLED RATHER THAN INVENTED
# The SEES prospective design is why this model works at all: median first draw
# at 6 days, 76% of first specimens before the population peak, up to 7 visits.
# Shigella failed precisely because its first draw came after the peak. A
# simulation on a synthetic schedule would answer a question about a study that
# does not exist. Whole subjects' visit vectors are resampled with replacement,
# so the joint structure (how many visits, how spaced) is preserved.
#
# WHY MEASUREMENT ERROR IS INCLUDED
# prec_logy is carried from the fit. Simulating without it makes recovery look
# far better than it is -- the precision goes from about 4 to enormous and every
# parameter appears identified.
#
# ENVIRONMENT VARIABLES
#   SEED         replicate number, also the RNG seed          default 1
#   RHO_TRUE     true cross-correlation, all five             default 0.4
#   TRUTH_FIT    fit to take the remaining truth from         default fit_ch2_461_R.rds
#   TRUTH_CHAINS "positive" keeps only the chains in the dominant baseline mode
#                (arm R has 4 of 40 stuck 110 lp units below the rest)
#   MODEL        stan file to FIT with                        default model_ch2_v5.stan
#   DATA         real data, used only for the visit design    default sees_subset_v2.rds
#   NCHAIN NADAPT NITER ADAPT_DELTA MAX_TD INIT C_PRIOR_SD
#   OUTDIR                                                    default sim_out
#
# USAGE
#   SEED=1 Rscript sim_ch2_seed.R
# =====================================================================

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
})

ge <- function(k, d) {
  v <- Sys.getenv(k)
  if (nzchar(v)) v else d
}
gn <- function(k, d) as.numeric(ge(k, as.character(d)))

SEED <- as.integer(gn("SEED", 1))
RHO_TRUE <- gn("RHO_TRUE", 0.4)
TRUTH_FIT <- ge("TRUTH_FIT", "fit_ch2_461_R.rds")
TRUTH_CH <- ge("TRUTH_CHAINS", "positive")
MODEL <- ge("MODEL", "model_ch2_v5.stan")
DATA <- ge("DATA", "sees_subset_v2.rds")
FUNCS <- ge("FUNCS", "stan_ch2_functions.R")
OUTDIR <- ge("OUTDIR", "sim_out")
NCHAIN <- as.integer(gn("NCHAIN", 8))
NADAPT <- as.integer(gn("NADAPT", 3000))
NITER <- as.integer(gn("NITER", 1000))
ADAPT_D <- gn("ADAPT_DELTA", 0.9)
MAX_TD <- as.integer(gn("MAX_TD", 12))
INIT <- gn("INIT", 0.1)
C_SD <- gn("C_PRIOR_SD", 1)
# PREC_MULT scales the observation precision carried from the fit. The default
# of 1 changes nothing. Larger values make each subject's theta better
# determined, which is how we ask whether the sd bias is a prior problem or an
# information problem.
PREC_MULT <- gn("PREC_MULT", 1)
# WISHDF is passed through to prep_ch2_priors. It sets BOTH the width of the
# lognormal prior on the sds and the LKJ eta, so a change here moves two things.
WISHDF <- gn("WISHDF", 20)
# TRUTH_C = "fit" takes the five cross-covariances straight from the fit
# instead of building them from a single scalar. That is what makes the truth
# equal to what we report, rather than an arbitrary grid point.
TRUTH_C <- ge("TRUTH_C", "scalar")
# LP_DROP is the log-density gap that separates "trapped" from "sampling", used
# when TRUTH_CHAINS = "lp". lp__ itself fluctuates with SD about sqrt(d/2), and
# d is 4,657 here, so chains legitimately in the same mode can sit 50 apart. A
# threshold of 20 would cut into the dominant mode; 100 is about two such SDs,
# which keeps it and still drops the trapped chains. Read the printed table and
# gap column before trusting the default.
LP_DROP <- gn("LP_DROP", 100)
# DRY_RUN=1 stops after the data are simulated and checked, before any sampling.
# Thirty seconds of work that can save a fourteen-hour round.
DRY_RUN <- ge("DRY_RUN", "0") == "1"

say <- function(...) {
  cat(sprintf("[sim %03d] ", SEED), ..., "\n", sep = "")
  flush.console()
}
dir.create(OUTDIR, showWarnings = FALSE)
if (!nzchar(Sys.getenv("CMDSTAN_EXEC_DIR"))) {
  Sys.setenv(CMDSTAN_EXEC_DIR = file.path(tempdir(), sprintf("sim_exe_%03d", SEED)))
}

for (f in c(TRUTH_FIT, MODEL, DATA, FUNCS)) {
  if (!file.exists(f)) stop("missing: ", f, call. = FALSE)
}

P <- 5L
K <- 2L # 5 kinetic parameters, 2 biomarkers
say(
  "seed ", SEED, " | rho_true ", RHO_TRUE, " | truth_c ", TRUTH_C,
  " | prec_mult ", PREC_MULT, " | wishdf ", WISHDF, " | init ", INIT,
  " | fit with ", MODEL
)

# ═══════════════════════════════════════════════════ 1. truth from the fit ══
fit <- readRDS(TRUTH_FIT)

keep_chains <- NULL
if (identical(TRUTH_CH, "lp")) {
  # Filter on log posterior density rather than on the sign of the correlation.
  # A chain can report a positive baseline correlation and still be trapped:
  # arm T2 had chains at lp -5015 with cross_corr[1] = +0.18, which the sign
  # rule would have kept. lp is the quantity that says whether a chain is
  # sampling the posterior or sitting somewhere the posterior barely visits.
  lpc <- apply(fit$draws("lp__", format = "draws_array"), 2, median)
  cc1 <- apply(fit$draws("cross_corr[1]", format = "draws_array"), 2, median)
  o <- order(-lpc)
  say("chain lp and baseline, best first:")
  for (i in o) {
    say(sprintf(
      "    chain %2d   lp %9.0f   %+.0f from best   baseline %+.3f",
      i, lpc[i], lpc[i] - max(lpc), cc1[i]
    ))
  }
  gaps <- -diff(lpc[o])
  big <- which(gaps > 3 * sqrt(0.5 * (10 + 5 + 5 + 10 + 10 + 5 + 2))) # rough scale
  if (length(big)) {
    say("  largest gaps after chain ", paste(o[big[order(-gaps[big])][1:min(3, length(big))]],
      collapse = ", "
    ))
  }
  keep_chains <- which(lpc > max(lpc) - LP_DROP)
  say(
    "truth chains: ", length(keep_chains), " of ", length(lpc),
    " (lp within ", LP_DROP, " of the best)"
  )
  drop <- setdiff(seq_along(lpc), keep_chains)
  if (length(drop)) {
    say("  dropped: ", paste(sprintf(
      "chain %d (lp %.0f, base %+.2f)",
      drop, lpc[drop], cc1[drop]
    ), collapse = "; "))
  }
} else if (identical(TRUTH_CH, "positive")) {
  cc1 <- as_draws_df(fit$draws("cross_corr[1]"))
  bych <- cc1 |>
    group_by(.chain) |>
    summarise(m = median(`cross_corr[1]`), .groups = "drop")
  keep_chains <- bych$.chain[bych$m > 0]
  say(
    "truth chains: ", length(keep_chains), " of ", nrow(bych),
    " (dropping the minor baseline mode)"
  )
}

med <- function(v) {
  d <- fit$draws(v, format = "draws_array")
  if (!is.null(keep_chains)) d <- d[, keep_chains, , drop = FALSE]
  out <- apply(d, 3, median)
  names(out) <- dimnames(d)[[3]]
  out
}

# CmdStanR flattens `array[K] vector[P] mu_par` with the FIRST index varying
# FASTEST: mu_par[1,1], mu_par[2,1], mu_par[1,2], mu_par[2,2], ...
# Positional slicing therefore interleaves the two biomarkers and produces a
# scrambled truth -- IgG would get IgA's baseline as its peak, and the decay
# slot would receive log t1. Index by name instead, and never by position.
mp <- med("mu_par")
mu_G <- unname(mp[sprintf("mu_par[1,%d]", 1:P)]) # slot 1 = IgG
mu_A <- unname(mp[sprintf("mu_par[2,%d]", 1:P)]) # slot 2 = IgA
if (anyNA(mu_G) || anyNA(mu_A)) {
  stop("mu_par names did not resolve; found: ", paste(names(mp), collapse = ", "),
    call. = FALSE
  )
}
sd_G <- unname(med("sd_G"))
sd_A <- unname(med("sd_A")) # CONDITIONAL on IgG -- see the Schur note
prec <- unname(med("prec_logy")) * PREC_MULT
if (PREC_MULT != 1) {
  say(
    "prec_logy scaled by ", PREC_MULT, " -> ", paste(sprintf("%.2f", prec), collapse = " "),
    "  (obs sd ", paste(sprintf("%.3f", 1 / sqrt(prec)), collapse = " "), ")"
  )
}
# Lcorr is a matrix, and matrix() fills column-major, which matches CmdStanR's
# Lcorr_G[1,1], Lcorr_G[2,1], ... ordering. This one is safe as written.
Lg_raw <- matrix(unname(med("Lcorr_G")), P, P)
La_raw <- matrix(unname(med("Lcorr_A")), P, P)

# An element-wise median of a Cholesky factor need not be one. Rebuild the
# correlation matrix, renormalise, and re-factor; stop if that fails.
fix_corr <- function(L, nm) {
  R <- L %*% t(L)
  d <- sqrt(diag(R))
  R <- R / outer(d, d)
  diag(R) <- 1
  ev <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 1e-8) {
    say("  ", nm, ": min eigenvalue ", signif(ev, 3), " -- projecting to PD")
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  }
  chol(R) |> t() # lower triangular
}
Lg <- fix_corr(Lg_raw, "Lcorr_G")
La <- fix_corr(La_raw, "Lcorr_A")

# Joint Cholesky, exactly as model_ch2_v4/v5 builds it:
#   L_full = [[L_G, 0], [diag(c) inv(L_G)', L_A|G]],  c_j = rho_j sd_G_j sd_A_j
L_G <- diag(sd_G) %*% Lg
L_Ac <- diag(sd_A) %*% La
if (identical(TRUTH_C, "fit")) {
  c_vec <- unname(med("cross_cov"))
  say("truth c taken from the fit: ", paste(sprintf("%.4f", c_vec), collapse = " "))
} else {
  c_vec <- RHO_TRUE * sd_G * sd_A
}
B <- diag(c_vec) %*% t(solve(L_G))
L_full <- rbind(cbind(L_G, matrix(0, P, P)), cbind(B, L_Ac))
Sig <- L_full %*% t(L_full)
stopifnot(min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values) > 0)

# The cross-correlations the fitted model will be scored against are the
# MARGINAL ones, not rho. cc_j = rho / sqrt(1 + g_j rho^2).
S_G <- L_G %*% t(L_G)
g <- diag(S_G) * diag(solve(S_G))
# Read the marginal cross-correlations straight off Sigma. The scalar formula
# cc = rho/sqrt(1+g rho^2) only holds when every c_j came from the same rho, so
# it cannot be used once the c's are taken from the fit.
cc_true <- diag(Sig[1:P, (P + 1):(2 * P)]) /
  sqrt(diag(Sig)[1:P] * diag(Sig)[(P + 1):(2 * P)])
ceiling_j <- 1 / sqrt(1 + g)

# Print the truth on the natural scale so a scrambled extraction is visible
# before anything expensive happens.
for (kk in 1:K) {
  p <- if (kk == 1) mu_G else mu_A
  ly1 <- log(exp(p[1]) + exp(p[2]))
  a <- exp(p[5])
  say(sprintf(
    "truth %s: y0 %.1f | y1 %.1f | t1 %.1f d | shape %.2f | k %.4f/d",
    c("IgG", "IgA")[kk], exp(p[1]), exp(ly1), exp(p[3]), a + 1, exp(p[4])
  ))
}
say("  (Aiemjoy 2022 for scale: IgG y0 14.5, y1 226, t1 15.6 d, r 2.18)")

say("true cross_corr (marginal): ", paste(sprintf("%.3f", cc_true), collapse = " "))
say("  ceiling at rho = 1      : ", paste(sprintf("%.3f", ceiling_j), collapse = " "))
say("  headroom                : ", paste(sprintf("%.0f%%", 100 * cc_true / ceiling_j),
  collapse = " "
))
if (any(cc_true > 0.9 * ceiling_j)) {
  say("  WARNING: a true value is within 10% of the ceiling; recovery will be poor there")
}

mu_flat <- c(mu_G, mu_A)
truth <- list(
  seed = SEED, rho_true = RHO_TRUE, truth_c = TRUTH_C,
  prec_mult = PREC_MULT, wishdf = WISHDF,
  c_vec = c_vec, mu = mu_flat, sd_G = sd_G,
  sd_A = sd_A, prec_logy = prec, cc_true = cc_true, g = g,
  ceiling = ceiling_j, Sigma = Sig,
  source = TRUTH_FIT, chains_kept = keep_chains
)
saveRDS(truth, file.path(OUTDIR, sprintf("truth_seed%03d.rds", SEED)))

# ═════════════════════════════════════════════ 2. design from the real data ══
real <- readRDS(DATA)
cols <- list(
  id = "id", t = "dayssincefeveronset", y = "result",
  bio = "antigen_iso"
)
sched <- real |>
  group_by(.data[[cols$id]]) |>
  summarise(t = list(sort(unique(.data[[cols$t]]))), .groups = "drop")
sched <- sched$t[lengths(sched$t) >= 2]
N <- length(sched)
say(
  "design: ", N, " subjects, visits ", min(lengths(sched)), "-", max(lengths(sched)),
  ", first draw median ", round(median(sapply(sched, min)), 1), " d"
)

set.seed(SEED)
pick <- sched[sample.int(N, N, replace = TRUE)]

# ═══════════════════════════════════════════════════════ 3. simulate data ══
# The model's slot 4 is log k, not log alpha. log_alpha = par[4] - a * log_y1.
log1pexp <- function(x) ifelse(x > 30, x, log1p(exp(x)))
curve_lp <- function(t, p) {
  ly0 <- p[1]
  ly1 <- log(exp(p[1]) + exp(p[2]))
  t1 <- exp(p[3])
  a <- exp(p[5])
  la <- p[4] - a * ly1
  ifelse(t <= t1,
    ly0 + (ly1 - ly0) / t1 * t,
    ly1 - log1pexp(log(a) + la + log(pmax(t - t1, 1e-12)) + a * ly1) / a
  )
}

Z <- matrix(rnorm(N * 2 * P), N, 2 * P)
PAR <- sweep(Z %*% t(L_full), 2, mu_flat, "+") # N x 10
sd_obs <- 1 / sqrt(prec)
iso <- c("HlyE_IgG", "HlyE_IgA") # slot 1 = IgG

rows <- vector("list", N)
for (i in seq_len(N)) {
  tt <- pick[[i]]
  out <- vector("list", K)
  for (k in seq_len(K)) {
    p <- PAR[i, ((k - 1) * P + 1):(k * P)]
    mu <- curve_lp(tt, p)
    out[[k]] <- data.frame(
      id = i, dayssincefeveronset = tt,
      antigen_iso = iso[k],
      result = exp(rnorm(length(tt), mu, sd_obs[k]))
    )
  }
  rows[[i]] <- do.call(rbind, out)
}
sim <- do.call(rbind, rows)
saveRDS(PAR, file.path(OUTDIR, sprintf("par_true_seed%03d.rds", SEED)))
say("simulated ", nrow(sim), " observations over ", N, " subjects")

# A silent failure mode: a mis-extracted truth produces data on the wrong
# scale, and nothing downstream complains -- the fit just answers a different
# question for fourteen hours. Compare per biomarker and STOP if the median is
# off by more than a factor of about two.
#
# The 5th percentile is expected to differ: the real assay is left-censored at
# its limit of detection and the simulation is not. The median is not censored,
# so that is what the gate uses.
real_bio <- gsub(" ", "_", as.character(real[[cols$bio]]))
gate_ok <- TRUE
for (b in iso) {
  qs <- quantile(log10(sim$result[sim$antigen_iso == b]), c(.05, .5, .95))
  qr <- quantile(log10(real[[cols$y]][real_bio == b]), c(.05, .5, .95), na.rm = TRUE)
  say(sprintf(
    "  log10 %s  simulated %.2f / %.2f / %.2f   real %.2f / %.2f / %.2f   median gap %+.2f",
    b, qs[1], qs[2], qs[3], qr[1], qr[2], qr[3], qs[2] - qr[2]
  ))
  if (abs(qs[2] - qr[2]) > 0.35) gate_ok <- FALSE
}
if (!gate_ok) {
  stop("simulated levels do not match the real data (median off by more than ",
    "0.35 on log10, i.e. a factor of 2.2).\n",
    "  The truth extraction is the first thing to check -- print mu_G and mu_A ",
    "and compare them with fit$summary(\"mu_par\").",
    call. = FALSE
  )
}

# ══════════════════════════════════════════════════════════════ 4. fit ══════
if (DRY_RUN) {
  # The dry run stops before sampling, so without this it would never touch the
  # fitting call and a wishdf mismatch would only surface hours later.
  source(FUNCS)
  has_w <- "wishdf_param" %in% names(formals(run_mod_stan_ch2))
  say("run_mod_stan_ch2() accepts wishdf_param: ", has_w)
  if (!has_w && WISHDF != 20) {
    stop("WISHDF=", WISHDF, " requested but run_mod_stan_ch2() cannot receive it.",
      call. = FALSE
    )
  }
  say("")
  say("DRY RUN: truth extracted, data simulated, level gate passed.")
  say("Re-run without DRY_RUN to fit.")
  quit(save = "no", status = 0)
}

source(FUNCS)

# run_mod_stan_ch2() may or may not accept wishdf_param depending on which
# version of stan_ch2_functions.R is on disk. Passing an argument a function
# does not have is an immediate error; NOT passing one it does have means the
# setting is silently ignored, which is how a whole night was spent last week
# fitting wishdf = 20 while believing it was 8. Check, then act, and refuse to
# run rather than run the wrong thing.
fit_args <- list(
  data = sim, file_mod = MODEL, estimate_c = 1,
  nchain = NCHAIN, parallel_chains = NCHAIN,
  nadapt = NADAPT, niter = NITER,
  adapt_delta = ADAPT_D, max_treedepth = MAX_TD,
  c_prior_sd = C_SD, init = INIT, seed = SEED, min_visits = 2L
)

if ("wishdf_param" %in% names(formals(run_mod_stan_ch2))) {
  fit_args$wishdf_param <- WISHDF
  say("wishdf ", WISHDF, " passed to run_mod_stan_ch2()")
} else if (WISHDF != 20) {
  stop("run_mod_stan_ch2() has no wishdf_param argument, so WISHDF=", WISHDF,
    " would be silently ignored and this run would fit wishdf = 20.\n",
    "  Add the argument to stan_ch2_functions.R and pass it through to\n",
    "  prep_ch2_priors(), or drop WISHDF from this call.",
    call. = FALSE
  )
} else {
  say("run_mod_stan_ch2() has no wishdf_param; using its own default (20)")
}

t0 <- Sys.time()
f <- do.call(run_mod_stan_ch2, fit_args)
hrs <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
say("fitted in ", round(hrs, 2), " hours")

# ══════════════════ SAVE FIRST. SCORE SECOND. ═══════════════════════════════
# The scoring block below is a dozen lines of tibble manipulation and it has
# already destroyed nine fits once: posterior::quantile2() names its column q5,
# not q05, so ss$q05 was NULL, ss$covered got logical(0), and the whole script
# exited -- taking CmdStan's CSVs in tempdir() with it. Seven hours per seed.
# Nothing downstream of the sampler is allowed to be the reason a fit is lost.
fit_file <- file.path(OUTDIR, sprintf("fit_seed%03d.rds", SEED))
f$save_object(fit_file)
say("[saved] ", fit_file)

# ---- scoring, in a box it cannot escape from ------------------------------
scored <- try(
  {
    vars <- c("mu_par", "sd_G", "sd_A", "prec_logy", "cross_corr")

    # Quantiles are computed straight from the draws. Passing them through
    # summarise_draws() means guessing what the column will be called, and the
    # returned vector's own names win over the argument name -- q5, or 95%, or
    # whatever the helper happens to produce. Take the draws and use quantile().
    dr <- f$draws(vars, format = "draws_matrix")
    ss <- f$summary(vars,
      mean = mean, median = stats::median, sd = stats::sd,
      rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
    )
    stopifnot(identical(as.character(ss$variable), colnames(dr)))
    ss$q05 <- apply(dr, 2, stats::quantile, 0.05)
    ss$q95 <- apply(dr, 2, stats::quantile, 0.95)

    truth_named <- c(
      setNames(mu_G, sprintf("mu_par[1,%d]", 1:P)),
      setNames(mu_A, sprintf("mu_par[2,%d]", 1:P)),
      setNames(sd_G, sprintf("sd_G[%d]", 1:P)),
      setNames(sd_A, sprintf("sd_A[%d]", 1:P)),
      setNames(prec, sprintf("prec_logy[%d]", 1:K)),
      setNames(cc_true, sprintf("cross_corr[%d]", 1:P))
    )
    ss$truth <- unname(truth_named[as.character(ss$variable)])
    if (anyNA(ss$truth)) {
      stop(
        "truth did not match: ",
        paste(ss$variable[is.na(ss$truth)], collapse = ", ")
      )
    }

    ss$bias <- ss$mean - ss$truth
    ss$sq_err <- ss$bias^2
    ss$post_var <- ss$sd^2
    ss$covered <- ss$truth >= ss$q05 & ss$truth <= ss$q95
    ss$seed <- SEED
    ss$rho_true <- RHO_TRUE

    dg <- f$diagnostic_summary()
    tot <- NCHAIN * NITER
    ss$n_div <- sum(dg$num_divergent)
    ss$pct_td <- 100 * sum(dg$num_max_treedepth) / tot
    ss$n_lowebfmi <- sum(dg$ebfmi < 0.3)
    ss$hours <- hrs

    saveRDS(as.data.frame(ss), file.path(OUTDIR, sprintf("res_seed%03d.rds", SEED)))

    d1 <- f$draws("cross_corr[1]", format = "draws_matrix")
    bych <- tapply(as.numeric(d1), rep(seq_len(NCHAIN), each = nrow(d1) / NCHAIN), median)
    saveRDS(bych, file.path(OUTDIR, sprintf("chainmode_seed%03d.rds", SEED)))

    list(ss = ss, bych = bych, dg = dg)
  },
  silent = FALSE
)

if (inherits(scored, "try-error")) {
  say("")
  say("SCORING FAILED -- but the fit is safe at ", fit_file)
  say("Re-score with score_from_fit.R once the cause is fixed.")
  quit(save = "no", status = 0) # 0: the expensive part succeeded
}

ss <- scored$ss
bych <- scored$bych
say("")
say("=== seed ", SEED, " ===")
say(sprintf(
  "  divergent %d | treedepth %.1f%% | E-BFMI<0.3 %d/%d | max R-hat %.3f | min ESS %.0f",
  ss$n_div[1], ss$pct_td[1], ss$n_lowebfmi[1], NCHAIN,
  max(ss$rhat, na.rm = TRUE), min(ss$ess_bulk, na.rm = TRUE)
))
say("  baseline corr by chain: ", paste(sprintf("%+.2f", bych), collapse = " "))
say(sprintf("  chains in the negative mode: %d of %d", sum(bych < 0), NCHAIN))
cc <- ss[grepl("cross_corr", ss$variable), ]
lab <- c("baseline", "peak", "peak-time", "decay-rate", "decay-shape")
say("")
say("  cross_corr recovery")
cat(sprintf(
  "    %-14s %7s %8s %8s %9s %8s\n",
  "parameter", "truth", "mean", "bias", "90% cov", "R-hat"
))
for (i in seq_len(nrow(cc))) {
  cat(sprintf(
    "    %-14s %7.3f %8.3f %+8.3f %9s %8.3f\n",
    lab[i], cc$truth[i], cc$mean[i], cc$bias[i],
    ifelse(cc$covered[i], "yes", "NO"), cc$rhat[i]
  ))
}
say("")
say("done -> ", file.path(OUTDIR, sprintf("res_seed%03d.rds", SEED)))
