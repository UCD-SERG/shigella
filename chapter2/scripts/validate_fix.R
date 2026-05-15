# ==========================================================================
# validate_fix.R — Verify model fixes work (v3 with init function)
#
# Goal: Run 3 small fits with FIXED model_2.stan + priors + explicit init
# function
#
# If rho_B is recovered near +0.6 for A and B, near 0 for C -> fix works.
#
# Run:
#   conda activate r_chapter2
#   cd ~/chapter2
#   rm -rf /tmp/$USER/cmdstan_bin*
#   Rscript scripts/validate_fix.R 2>&1 | tee logs/validate_fix.log
# ==========================================================================

setwd("~/chapter2")

cat("\n========================================================\n")
cat(" VALIDATE FIX v3: 3 quick fits with JAGS-aligned model\n")
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

# ==========================================================================
# Init function — start chains from a sensible point near JAGS truth
# ==========================================================================
make_init <- function(N, K = 2, P = 5) {
  function() {
    list(
      # M near the population mean (matches JAGS prior mean)
      M = matrix(c(1, 7, 1, -4, -1), K, P, byrow = TRUE) +
          matrix(rnorm(K * P, 0, 0.1), K, P),

      # Cholesky factors near identity (slight perturbation)
      L_Omega_B   = diag(K),
      L_Omega_P   = diag(P),
      L_Omega_eps = diag(K),

      # Moderate tau values
      tau_B   = rep(0.5, K),
      tau_P   = rep(0.3, P),
      tau_eps = rep(0.3, K),

      # Z close to zero — non-centered helper
      Z = matrix(rnorm(N * K * P, 0, 0.1), N, K * P)
    )
  }
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
    omega_B           = Omega_B_true,
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
      init            = make_init(scn$n),  # Use init function
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

  # Symmetry check
  omega_B_21 <- posterior::as_draws_df(sf$draws(variables = "Omega_B[2,1]"))
  rho_B_21 <- omega_B_21[["Omega_B[2,1]"]]

  diag <- sf$diagnostic_summary()
  n_div <- sum(diag$num_divergent)
  n_td  <- sum(diag$num_max_treedepth)

  est_med <- median(rho_B_post)
  est_lo  <- quantile(rho_B_post, 0.025, names = FALSE)
  est_hi  <- quantile(rho_B_post, 0.975, names = FALSE)
  est_med_21 <- median(rho_B_21)

  results[[scn$name]] <- list(
    scenario      = scn$name,
    status        = "OK",
    n             = scn$n,
    true          = scn$rho_B,
    est_12_med    = est_med,
    est_21_med    = est_med_21,
    est_lo        = est_lo,
    est_hi        = est_hi,
    bias          = est_med - scn$rho_B,
    n_divergent   = n_div,
    n_treedepth   = n_td,
    elapsed_min   = elapsed,
    symmetric_check = abs(est_med - est_med_21) < 0.01
  )

  cat(sprintf("  Estimate [1,2]: %+.3f [%.3f, %.3f]\n",
              est_med, est_lo, est_hi))
  cat(sprintf("  Estimate [2,1]: %+.3f (should ~equal [1,2])\n", est_med_21))
  cat(sprintf("  Bias: %+.3f\n", est_med - scn$rho_B))
  cat(sprintf("  Divergent: %d / 4000 (%.1f%%)\n",
              n_div, 100 * n_div / 4000))
  cat(sprintf("  Treedepth hits: %d / 4000 (%.1f%%)\n",
              n_td, 100 * n_td / 4000))
  cat(sprintf("  Elapsed: %.1f min\n", elapsed))
}

# ==========================================================================
# Summary
# ==========================================================================
cat("\n\n========================================================\n")
cat(" VALIDATION SUMMARY\n")
cat("========================================================\n\n")

cat(sprintf("%-10s %-6s %-12s %-25s %-15s\n",
            "Scenario", "n", "True rho_B", "Est rho_B [95% CI]",
            "Divergent (%)"))
cat(strrep("-", 75), "\n")

for (s in names(results)) {
  r <- results[[s]]
  if (r$status == "FAILED") {
    cat(sprintf("%-10s FAILED\n", s))
    next
  }
  cat(sprintf("%-10s %-6d %-12.2f %+.3f [%+.3f, %+.3f]   %d (%.1f%%)\n",
              s, r$n, r$true,
              r$est_12_med, r$est_lo, r$est_hi,
              r$n_divergent, 100 * r$n_divergent / 4000))
}

cat("\n=== DIAGNOSIS ===\n")
all_ok <- all(sapply(results, function(r) {
  if (r$status != "OK") return(FALSE)
  abs(r$est_12_med - r$true) < 0.3 && r$n_divergent < 1000
}))

if (all_ok) {
  cat("\n[OK] ALL SCENARIOS RECOVERED — Fix works!\n")
  cat("     Next: mini array test, then full R=200 run.\n\n")
} else {
  any_improved <- any(sapply(results, function(r) {
    if (r$status != "OK") return(FALSE)
    # Improvement from -0.83 bias
    r$est_12_med > -0.3
  }))

  if (any_improved) {
    cat("\n[PARTIAL] Some improvement vs original (-0.83 bias).\n")
    cat("          Check details. May need further tuning.\n\n")
  } else {
    cat("\n[FAIL] Still problematic.\n")
    cat("       Options:\n")
    cat("       1. Increase iter (warmup 2000 + sampling 2000)\n")
    cat("       2. Reparameterize (per-biomarker hierarchical)\n")
    cat("       3. Reduce complexity (fix Omega_P = identity)\n\n")
  }
}

saveRDS(results, "outputs/validate_fix_results.rds")
cat("Results saved to outputs/validate_fix_results.rds\n")
