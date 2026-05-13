# ==========================================================================
# validate_fix_v2.R — Use Stan default init (no explicit Cholesky factors)
#
# Change from v1: removed explicit `L_Omega_B = diag(K)` init.
# Use init = 0.5 (Stan's default random init around 0 in unconstrained space).
#
# Rationale: passing identity matrix to a Cholesky factor parameter causes
# issues when Stan transforms back to unconstrained space — diagonal entries
# may map to 0, triggering lkj_corr_cholesky_lpdf exceptions.
#
# Run:
#   conda activate r_chapter2
#   cd ~/chapter2
#   rm -rf /tmp/$USER/cmdstan_bin*
#   Rscript scripts/validate_fix_v2.R 2>&1 | tee logs/validate_fix_v2.log
# ==========================================================================

setwd("~/chapter2")

cat("\n========================================================\n")
cat(" VALIDATE FIX v2: Stan default init (no explicit Cholesky)\n")
cat("========================================================\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(serodynamics)
  library(cmdstanr)
  library(posterior)
})

source("R/prep_data_stan.R")
source("R/prep_priors_stan.R")
source("R/postprocess_stan_output.R")
source("R/run_mod_stan.R")
source("R/sim_correlated_case_data.R")

make_omega_2x2 <- function(rho) {
  matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
}

scenarios <- list(
  list(name = "A", n = 48, rho_B = 0.6),
  list(name = "B", n = 11, rho_B = 0.6),
  list(name = "C", n = 48, rho_B = 0.0)
)

results <- list()

for (scn in scenarios) {
  cat(sprintf("\n=== Scenario %s: n=%d, true rho_B=%.1f ===\n",
              scn$name, scn$n, scn$rho_B))

  set.seed(2026 + which(sapply(scenarios, function(x) x$name) == scn$name))
  Omega_B_true <- make_omega_2x2(scn$rho_B)

  sim_dat <- sim_correlated_case_data(
    n                 = scn$n,
    Omega_B           = Omega_B_true,
    antigen_isos      = c("IgG", "IgA"),
    n_obs_per_subject = 5L
  )
  cat("  Simulated", nrow(sim_dat), "rows\n")

  t0 <- Sys.time()
  fit <- tryCatch({
    run_mod_stan(
      data            = sim_dat,
      model           = "model_2",
      chains          = 4,
      iter_warmup     = 1000,
      iter_sampling   = 1000,
      parallel_chains = 4,
      adapt_delta     = 0.95,
      max_treedepth   = 12,
      init            = 0.5,             # Stan default 
      with_post       = TRUE,
      stan_dir        = "inst/stan",
      refresh         = 250,
      show_messages   = FALSE
    )
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(fit)) {
    results[[scn$name]] <- list(scenario = scn$name, status = "FAILED")
    next
  }

  elapsed <- as.numeric(Sys.time() - t0, units = "mins")

  sf <- attr(fit, "stan_fit")[[1]]
  omega_B_draws <- posterior::as_draws_df(sf$draws(variables = "Omega_B[1,2]"))
  rho_B_post <- omega_B_draws[["Omega_B[1,2]"]]

  omega_B_21 <- posterior::as_draws_df(sf$draws(variables = "Omega_B[2,1]"))
  rho_B_21 <- omega_B_21[["Omega_B[2,1]"]]

  diag <- sf$diagnostic_summary()
  n_div <- sum(diag$num_divergent)
  n_td  <- sum(diag$num_max_treedepth)

  est_med <- median(rho_B_post)
  est_lo  <- quantile(rho_B_post, 0.025, names = FALSE)
  est_hi  <- quantile(rho_B_post, 0.975, names = FALSE)

  results[[scn$name]] <- list(
    scenario    = scn$name,
    status      = "OK",
    n           = scn$n,
    true        = scn$rho_B,
    est_12_med  = est_med,
    est_21_med  = median(rho_B_21),
    est_lo      = est_lo,
    est_hi      = est_hi,
    bias        = est_med - scn$rho_B,
    n_divergent = n_div,
    n_treedepth = n_td,
    elapsed_min = elapsed
  )

  cat(sprintf("  Estimate [1,2]: %+.3f [%.3f, %.3f]\n",
              est_med, est_lo, est_hi))
  cat(sprintf("  Bias: %+.3f\n", est_med - scn$rho_B))
  cat(sprintf("  Divergent: %d / 4000 (%.1f%%)\n",
              n_div, 100 * n_div / 4000))
  cat(sprintf("  Elapsed: %.1f min\n", elapsed))

  # Save intermediate after each scenario in case script is interrupted
  saveRDS(results, "outputs/validate_fix_v2_results.rds")
}

# Summary
cat("\n\n========================================================\n")
cat(" VALIDATION SUMMARY v2\n")
cat("========================================================\n\n")

cat(sprintf("%-10s %-6s %-12s %-25s %-15s\n",
            "Scenario", "n", "True rho_B", "Est rho_B [95% CI]",
            "Divergent"))
cat(strrep("-", 75), "\n")

for (s in names(results)) {
  r <- results[[s]]
  if (r$status == "FAILED") {
    cat(sprintf("%-10s FAILED\n", s))
    next
  }
  cat(sprintf("%-10s %-6d %-12.2f %+.3f [%+.3f, %+.3f]  %d (%.1f%%)\n",
              s, r$n, r$true,
              r$est_12_med, r$est_lo, r$est_hi,
              r$n_divergent, 100 * r$n_divergent / 4000))
}

cat("\nResults saved to outputs/validate_fix_v2_results.rds\n")
