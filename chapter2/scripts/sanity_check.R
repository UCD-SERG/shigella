# ==========================================================================
# sanity_check.R — Shiva-compatible
#
# Quick end-to-end test before running full Phase 2.
# Verifies cmdstanr works, simulation pipeline works, posterior extraction works.
#
# Run:
#   conda activate r_chapter2
#   cd ~/chapter2
#   Rscript scripts/sanity_check.R 2>&1 | tee logs/sanity_check.log
# ==========================================================================

setwd("~/chapter2")

cat("\n=========================================================\n")
cat(" SANITY CHECK: Chapter 2 Phase 2 Pipeline (Shiva)\n")
cat("=========================================================\n\n")

# ---- 1. Load packages ----
cat("[1/8] Loading packages...\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(serodynamics)
  library(cmdstanr)
  library(posterior)
  library(cli)
  library(tibble)
})
cat("      OK\n\n")

# ---- 2. Source local helpers ----
cat("[2/8] Sourcing R/ helpers...\n")
source("R/prep_data_stan.R")
source("R/prep_priors_stan.R")
source("R/postprocess_stan_output.R")
source("R/run_mod_stan.R")
source("R/sim_correlated_case_data.R")
cat("      OK\n\n")

# ---- 3. Verify Stan files + compile dir ----
cat("[3/8] Checking Stan files and compile directory...\n")
stan_files <- c("inst/stan/model_1.stan",
                "inst/stan/model_2.stan",
                "inst/stan/model_1_time_est.stan",
                "inst/stan/model_2_time_est.stan")
for (sf in stan_files) {
  if (file.exists(sf)) {
    cat("      [OK]   ", sf, "\n")
  } else {
    cat("      [MISS] ", sf, "\n")
  }
}

user <- Sys.getenv("USER", unset = "default")
compile_dir <- Sys.getenv("STAN_COMPILE_DIR",
                          unset = file.path("/tmp", user, "cmdstan_bin"))
if (!dir.exists(compile_dir)) {
  dir.create(compile_dir, recursive = TRUE, mode = "0755")
}
cat("      Compile dir:", compile_dir, "\n\n")

# ---- 4. Test simulation ----
cat("[4/8] Testing sim_correlated_case_data() ...\n")
test_dat <- sim_correlated_case_data(
  n              = 5,
  Omega_B        = matrix(c(1, 0.6, 0.6, 1), 2, 2),
  antigen_isos   = c("IgG", "IgA"),
  n_obs_per_subject = 3L
)
cat("      Class:", paste(class(test_dat), collapse = ", "), "\n")
cat("      Rows:", nrow(test_dat), "\n")
cat("      Subjects:", length(unique(test_dat$id)), "\n")
cat("      OK\n\n")

# ---- 5. Test prep_data + prep_data_stan ----
cat("[5/8] Testing prep_data + prep_data_stan ...\n")
prepped <- serodynamics::prep_data(test_dat)
stan_data <- prep_data_stan(prepped)
cat("      Stan data keys:", paste(names(stan_data), collapse = ", "), "\n")
cat("      N=", stan_data$N, ", K=", stan_data$K, ", P=", stan_data$P, "\n")
cat("      log_y dim:", paste(dim(stan_data$log_y), collapse = "x"), "\n")
cat("      OK\n\n")

# ---- 6. Test prep_priors_stan ----
cat("[6/8] Testing prep_priors_stan ...\n")
priors_2 <- prep_priors_stan(model = "model_2")
cat("      Model 2 prior keys:", paste(names(priors_2), collapse = ", "), "\n")
cat("      OK\n\n")

# ---- 7. Test full run_mod_stan with Model 2 (small chain) ----
cat("[7/8] Testing run_mod_stan(model='model_2') (this takes 3-5 min)...\n")
t0 <- Sys.time()
fit2 <- tryCatch({
  run_mod_stan(
    data            = test_dat,
    model           = "model_2",
    chains          = 2,
    iter_warmup     = 250,
    iter_sampling   = 250,
    parallel_chains = 2,
    with_post       = TRUE,
    stan_dir        = "inst/stan",
    compile_dir     = compile_dir,
    refresh         = 100,
    show_messages   = FALSE
  )
}, error = function(e) {
  cat("\n      [ERROR]:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(fit2)) {
  elapsed <- as.numeric(Sys.time() - t0, units = "mins")
  cat(sprintf("\n      Fit completed in %.1f min\n", elapsed))
  cat("      sr_model class:", paste(class(fit2), collapse = ", "), "\n")

  if (!is.null(attr(fit2, "Omega_B"))) {
    omega_B <- attr(fit2, "Omega_B")
    cat("      Omega_B[1,2] (should be near 0.6):",
        round(omega_B[1, 2], 3), "\n")
  } else {
    cat("      WARNING: Omega_B not in attributes!\n")
  }
  cat("      OK\n\n")
} else {
  cat("      FAILED — fix errors above before proceeding\n\n")
  stop("Sanity check failed at Step 7")
}

# ---- 8. Test direct posterior extraction (cmdstanr style) ----
cat("[8/8] Testing direct posterior extraction (Omega_B[1,2])...\n")
sf <- attr(fit2, "stan_fit")[[1]]
omega_B_draws <- posterior::as_draws_df(sf$draws(variables = "Omega_B[1,2]"))
rho_B_post <- omega_B_draws[["Omega_B[1,2]"]]
cat(sprintf("      rho_B posterior: median=%.3f [%.3f, %.3f]\n",
            median(rho_B_post),
            quantile(rho_B_post, 0.025),
            quantile(rho_B_post, 0.975)))
cat(sprintf("      n posterior draws: %d\n", length(rho_B_post)))
cat("      OK\n\n")

# ---- Summary ----
cat("=========================================================\n")
cat(" ALL SANITY CHECKS PASSED\n")
cat("=========================================================\n")
cat(" - cmdstanr backend works\n")
cat(" - Compile dir resolves (no /home noexec issue)\n")
cat(" - Pipeline ready for Phase 2 simulation\n")
cat("=========================================================\n")
cat("\nNext step: scripts/02_run_scenarios.R for pilot,\n")
cat("           or sbatch slurm/run_phase2_array.sbatch for full study.\n")
