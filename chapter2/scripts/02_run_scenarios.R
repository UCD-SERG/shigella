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
  library(serodynamics)
  library(cmdstanr)        # Not rstan
  library(posterior)       # for as_draws_df, as_draws_array
  library(cli)
  library(tibble)
})

# Source local helpers (cmdstanr-compatible versions)
source("R/prep_data_stan.R")
source("R/prep_priors_stan.R")
source("R/postprocess_stan_output.R")
source("R/run_mod_stan.R")
source("R/sim_correlated_case_data.R")

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

make_omega_2x2 <- function(rho) {
  matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
}

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
    cat(sprintf("\n[Scenario %s rep %d/%d] %s\n",
                s_name, rep, n_replicates, format(Sys.time())))
    t0 <- Sys.time()

    Omega_B_true <- make_omega_2x2(scn$rho_B)

    set.seed(2026 * 100 + rep)
    sim_dat <- tryCatch({
      sim_correlated_case_data(
        n                 = scn$n,
        Omega_B           = Omega_B_true,
        antigen_isos      = c("IgG", "IgA"),
        n_obs_per_subject = 5L
      )
    }, error = function(e) {
      cat(sprintf("  SIM ERROR: %s\n", conditionMessage(e))); NULL
    })
    if (is.null(sim_dat)) {
      scenario_results[[rep]] <- list(scenario = s_name, rep = rep,
                                       status = "SIM_FAILED")
      next
    }

    fit <- tryCatch({
      run_mod_stan(
        data            = sim_dat,
        model           = "model_2",
        chains          = n_chains,
        iter_warmup     = n_iter_warmup,
        iter_sampling   = n_iter_sample,
        parallel_chains = n_chains,
        adapt_delta     = 0.99,
        max_treedepth   = 12,
        init            = 0.1,
        with_post       = TRUE,
        stan_dir        = "inst/stan",
        refresh         = 200,
        show_messages   = FALSE
      )
    }, error = function(e) {
      cat(sprintf("  FIT ERROR: %s\n", conditionMessage(e))); NULL
    })
    if (is.null(fit)) {
      scenario_results[[rep]] <- list(scenario = s_name, rep = rep,
                                       status = "FIT_FAILED")
      next
    }

    rho_B_post <- tryCatch({
      sf <- attr(fit, "stan_fit")[[1]]
      # cmdstanr API: $draws() returns posterior draws_array
      omega_B_draws <- posterior::as_draws_df(
        sf$draws(variables = "Omega_B[1,2]")
      )
      omega_B_draws[["Omega_B[1,2]"]]
    }, error = function(e) {
      cat(sprintf("  EXTRACT ERROR: %s\n", conditionMessage(e))); NULL
    })
    if (is.null(rho_B_post)) {
      scenario_results[[rep]] <- list(scenario = s_name, rep = rep,
                                       status = "EXTRACT_FAILED")
      next
    }

    # Diagnostics
    n_divergent <- tryCatch({
      sf <- attr(fit, "stan_fit")[[1]]
      sum(sf$diagnostic_summary()$num_divergent)
    }, error = function(e) NA_integer_)

    elapsed <- as.numeric(Sys.time() - t0, units = "mins")

    scenario_results[[rep]] <- list(
      scenario          = s_name,
      rep               = rep,
      status            = "OK",
      true_rho_B        = scn$rho_B,
      est_rho_B_median  = median(rho_B_post),
      est_rho_B_mean    = mean(rho_B_post),
      est_rho_B_lo      = quantile(rho_B_post, 0.025, names = FALSE),
      est_rho_B_hi      = quantile(rho_B_post, 0.975, names = FALSE),
      bias              = median(rho_B_post) - scn$rho_B,
      n_divergent       = n_divergent,
      elapsed_min       = elapsed
    )

    cat(sprintf("  rho_B = %.3f [%.3f, %.3f], bias = %+.3f, %d div, %.1f min\n",
                median(rho_B_post),
                quantile(rho_B_post, 0.025, names = FALSE),
                quantile(rho_B_post, 0.975, names = FALSE),
                median(rho_B_post) - scn$rho_B,
                n_divergent,
                elapsed))

    saveRDS(scenario_results, sprintf("outputs/02_intermediate_%s.rds", s_name))
  }

  all_results[[s_name]] <- scenario_results
  saveRDS(all_results, "outputs/02_simulation_results.rds")
}

saveRDS(all_results, "outputs/02_simulation_results.rds")

cat("\n=== Phase 2 simulation complete ===\n")
cat("Output: outputs/02_simulation_results.rds\n")
cat("Run scripts/03_analyze_results.R next.\n")
