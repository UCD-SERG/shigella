# ==========================================================================
# 02_run_scenarios.R — Shiva-compatible (cmdstanr backend)
#
# Phase 2: Simulation study to verify Stan Model 2 recovers known Sigma_B.
#
# Pipeline:
#   sim_correlated_case_data() -> case_data
#   run_mod_stan(model = "model_2") -> sr_model with Omega_B attribute
#   Compare Omega_B[1,2] posterior to true rho_B
#
# KEY DIFFERENCES from Mercury version:
#   1. Uses cmdstanr::cmdstan_model() instead of rstan::stan()
#   2. Compiled binaries written to /tmp (avoids /home noexec)
#   3. Argument names: iter_sampling/iter_warmup (not iter/warmup)
#   4. Posterior extracted via $draws() not rstan::extract()
#
# For the pilot (R=20), use this script directly 
# For the full study (R=500), use 02_run_array.R + slurm/run_phase2_array.sbatch
#   to parallelize across SLURM job array tasks.
# ==========================================================================

setwd("~/chapter2")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(cmdstanr)        # Not rstan
  library(posterior)       # for as_draws_df, as_draws_array
  library(cli)
  library(tibble)
  library(shigella)
})

source("R/make_omega_2x2.R")
source("R/run_one_replicate.R")

set.seed(2026)

# ==========================================================================
# 1. Scenarios — same as Mercury
# ==========================================================================
scenarios <- list(
  A = list(n = 48, rho_B = 0.6, label = "A (n=48, rho=0.6)"),
  B = list(n = 11, rho_B = 0.6, label = "B (n=11, rho=0.6)"),
  C = list(n = 48, rho_B = 0.0, label = "C (null, n=48, rho=0)")
)

# ----- Phase 2 settings -----
# Pilot:  R=20, 4 chains, 1500 warmup + 1500 sampling = 3000 iter
# Per Ezra's MC SE feedback (5/6 meeting), final paper will use R=500
# via 02_run_array.R (SLURM array). This script is for pilot only.
n_replicates  <- 20
n_chains      <- 4
n_iter_warmup <- 1500
n_iter_sample <- 1500

# ==========================================================================
# 2. Run scenarios
# ==========================================================================
all_results <- list()

for (s_name in names(scenarios)) {
  scn <- scenarios[[s_name]]
  cat(sprintf("\n=== Scenario %s: n=%d, rho_B=%.1f ===\n",
              s_name, scn$n, scn$rho_B))

  scenario_results <- list()

  for (rep in 1:n_replicates) {
    scenario_results[[rep]] <- run_one_replicate(
      s_name        = s_name,
      scn           = scn,
      rep           = rep,
      n_chains      = n_chains,
      n_iter_warmup = n_iter_warmup,
      n_iter_sample = n_iter_sample
    )
    saveRDS(scenario_results, sprintf("outputs/02_intermediate_%s.rds", s_name))
  }

  all_results[[s_name]] <- scenario_results
  saveRDS(all_results, "outputs/02_simulation_results.rds")
}

saveRDS(all_results, "outputs/02_simulation_results.rds")

cat("\n=== Phase 2 simulation complete ===\n")
cat("Output: outputs/02_simulation_results.rds\n")
cat("Run scripts/03_analyze_results.R next.\n")
