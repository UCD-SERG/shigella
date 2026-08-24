#!/usr/bin/env Rscript
# =====================================================================
# chain_stacking.R -- weight chains by cross-validated predictive performance
# instead of by a log-density threshold.
#
# The authors' own script expects an rstan fit; ours are cmdstanr. Rather than
# adapt it, this uses the loo package directly, which implements both halves of
# the method: loo() gives the PSIS-LOO pointwise predictive densities per chain
# (their eq. 6), and loo_model_weights(method = "stacking") solves the simplex
# optimisation (their eq. 3). Treating each chain as a "model" is exactly what
# the paper proposes.
#
#   Yao, Vehtari & Gelman (2022) JMLR 23(79):1-45
#
# WHAT TO LOOK FOR. The lp filter and stacking should agree: chains sitting
# thousands of log-density units below the rest should receive weight near zero.
# Section 6.3 of the paper reports exactly that for a hierarchical model with a
# bad mode "orders of magnitude lower in posterior density".
#
#   Rscript chain_stacking.R fit_ch2_461_T.rds
#   Rscript chain_stacking.R sim_out_alt_prod/fit_seed003.rds
# =====================================================================

suppressPackageStartupMessages({
  library(posterior)
  library(loo)
})

path <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(path)) path <- "fit_ch2_461_T.rds"
LP_DROP <- as.numeric(Sys.getenv("LP_DROP", "100"))
LAMBDA <- as.numeric(Sys.getenv("LAMBDA", "1.0001"))
NM5 <- c("baseline", "peak", "peak-time", "decay-rate", "decay-shape")

cat(sprintf("=== %s ===\n\n", basename(path)))
fit <- readRDS(path)

lp <- apply(fit$draws("lp__", format = "draws_array"), 2, median)
K <- length(lp)
cc <- as.matrix(fit$draws("cross_corr", format = "draws_matrix"))
nit <- nrow(cc) / K
base_k <- apply(matrix(cc[, 1], nrow = nit, ncol = K), 2, median)

# ---- PSIS-LOO for each chain on its own -----------------------------------
# One chain at a time keeps peak memory at iter x nobs rather than
# iter x chain x nobs, which for 40 chains and 3,318 observations is the
# difference between 26 MB and 1 GB.
cat(sprintf("computing PSIS-LOO for %d chains", K))
loos <- vector("list", K)
kbad <- numeric(K)
for (k in seq_len(K)) {
  llk <- fit$draws("log_lik", format = "draws_array")[, k, , drop = FALSE]
  # relative_eff needs the chain structure; with one chain it is that chain's
  # own efficiency relative to independent draws.
  reff <- relative_eff(exp(llk), chain_id = rep(1, nit))
  loos[[k]] <- suppressWarnings(loo(llk, r_eff = reff))
  kbad[k] <- mean(loos[[k]]$diagnostics$pareto_k > 0.7)
  cat(".")
  if (k %% 10 == 0) cat(k)
}
cat(" done\n\n")

# ---- stacking weights -------------------------------------------------------
w <- loo_model_weights(loos, method = "stacking", optim_control = list(reltol = 1e-10))
w <- as.numeric(w)

# ---- compare with the lp rule ----------------------------------------------
keep <- lp > max(lp) - LP_DROP
tab <- data.frame(
  chain = seq_len(K),
  lp_gap = round(lp - max(lp)),
  baseline = round(base_k, 3),
  lp_keep = keep,
  stack_w = round(w, 5),
  pareto_bad = round(kbad, 2)
)
tab <- tab[order(-tab$stack_w), ]
print(tab, row.names = FALSE)

cat(sprintf("\n  lp rule keeps      %d of %d\n", sum(keep), K))
cat(sprintf("  stacking w > 0.001 %d of %d\n", sum(w > 0.001), K))
cat(sprintf("  effective chains   %.1f   (1 / sum w^2)\n", 1 / sum(w^2)))

agree <- sum((w > 0.001) == keep)
cat(sprintf("\n  the two rules agree on %d of %d chains\n", agree, K))
disc <- which((w > 0.001) != keep)
if (length(disc)) {
  for (i in disc) {
    cat(sprintf(
      "    chain %2d  lp %+6.0f  base %+.3f  lp_keep %s  stack_w %.5f\n",
      i, lp[i] - max(lp), base_k[i], keep[i], w[i]
    ))
  }
}

# ---- intervals under all three schemes -------------------------------------
# Importance resampling turns the weighted mixture into an unweighted sample, so
# ordinary quantiles apply (their Section 2.2).
set.seed(1)
S <- 20000
pick <- sample.int(K, S, replace = TRUE, prob = w)

cat(sprintf("\n  %-12s %26s %26s %26s\n", "", "all chains", "lp filter", "stacking"))
for (j in 1:5) {
  m <- matrix(cc[, j], nrow = nit, ncol = K)
  a <- as.vector(m)
  b <- as.vector(m[, keep, drop = FALSE])
  s <- m[cbind(sample.int(nit, S, replace = TRUE), pick)]
  f <- function(x) sprintf("%+.3f [%+.3f,%+.3f]", median(x), quantile(x, .05), quantile(x, .95))
  cat(sprintf("  %-12s %26s %26s %26s\n", NM5[j], f(a), f(b), f(s)))
}
cat(sprintf(
  "\n  %-12s %26.3f %26.3f %26.3f\n", "mean width",
  mean(sapply(1:5, function(j) {
    x <- cc[, j]
    diff(quantile(x, c(.05, .95)))
  })),
  mean(sapply(1:5, function(j) {
    m <- matrix(cc[, j], nrow = nit, ncol = K)
    x <- as.vector(m[, keep, drop = FALSE])
    diff(quantile(x, c(.05, .95)))
  })),
  mean(sapply(1:5, function(j) {
    m <- matrix(cc[, j], nrow = nit, ncol = K)
    x <- m[cbind(sample.int(nit, S, TRUE), pick)]
    diff(quantile(x, c(.05, .95)))
  }))
))

saveRDS(
  list(
    path = path, lp = lp, baseline = base_k, keep = keep,
    w = w, pareto_bad = kbad, lp_drop = LP_DROP
  ),
  sub("\\.rds$", "_stack.rds", basename(path))
)
cat(sprintf("\nsaved -> %s\n", sub("\\.rds$", "_stack.rds", basename(path))))

cat("\nPareto k: ", sprintf(
  "%.0f%% of observations exceed 0.7 on average across chains\n",
  100 * mean(kbad)
))
cat("The same diagnostic applies to the LOO model comparison; Yao et al. report\n")
cat("that stacking still outperforms other weightings when it is unreliable.\n")
