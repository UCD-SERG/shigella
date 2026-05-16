# ==========================================================================
# debug_pipeline.R — Layer-by-layer pipeline debug
#
# Goal: Find where the bug is between sim_correlated_case_data() and the
# Stan posterior of Omega_B[1,2].
#
# Strategy: Check each layer sees CONSISTENT rho_B = 0.6.
#
# Run:
#   conda activate r_chapter2
#   cd ~/chapter2
#   Rscript scripts/debug_pipeline.R 2>&1 | tee logs/debug_pipeline.log
# ==========================================================================

setwd("~/chapter2")

cat("\n========================================================\n")
cat(" DEBUG: Pipeline truth-vs-recovery audit\n")
cat("========================================================\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(cmdstanr)
  library(posterior)
  library(tibble)
  library(shigella)
})

source("R/debug_utils.R")

# ==========================================================================
# LAYER 1: Simulation function — does it actually generate correlated data?
# ==========================================================================
cat("\n###############################################\n")
cat("### LAYER 1: sim_correlated_case_data() check\n")
cat("###############################################\n\n")

set.seed(42)
TRUE_RHO <- 0.6
N_BIG <- 200   # Big n for clean correlation estimate

Omega_B_true <- matrix(c(1, TRUE_RHO, TRUE_RHO, 1), 2, 2)
cat("Generating n =", N_BIG, "subjects with TRUE rho_B =", TRUE_RHO, "\n\n")

sim_dat <- sim_correlated_case_data(
  n                 = N_BIG,
  omega_B           = Omega_B_true,
  antigen_isos      = c("IgG", "IgA"),
  n_obs_per_subject = 5L
)

cat("Output structure:\n")
print(head(sim_dat, 10))
cat("\nColumns:", paste(colnames(sim_dat), collapse = ", "), "\n")
cat("Rows:", nrow(sim_dat), "\n")
cat("Subjects:", length(unique(sim_dat$id)), "\n")
cat("Iso types:", paste(unique(sim_dat$antigen_iso), collapse = ", "), "\n\n")

# Look at the TRUE PARAMETERS attached to sim_dat (if available)
if (!is.null(attr(sim_dat, "truth"))) {
  cat("=== TRUTH attribute keys ===\n")
  truth <- attr(sim_dat, "truth")
  cat(names(truth), "\n\n")

  if ("Theta" %in% names(truth) || "theta" %in% names(truth)) {
    theta <- truth[["Theta"]] %||% truth[["theta"]]
    cat("Theta shape:", paste(dim(theta), collapse = " x "), "\n")
    # Theta should be [N, K, P] or similar
  }
}

# ==========================================================================
# LAYER 2: Recover correlation from simulation truth
# ==========================================================================
cat("\n###############################################\n")
cat("### LAYER 2: Empirical recovery from sim truth\n")
cat("###############################################\n\n")

diagnose_layer2(sim_dat, TRUE_RHO)

# ==========================================================================
# LAYER 3: prep_data_stan — does it preserve the correlated structure?
# ==========================================================================
cat("\n###############################################\n")
cat("### LAYER 3: prep_data_stan check\n")
cat("###############################################\n\n")

# Use a smaller sim for Stan testing
set.seed(42)
sim_small <- sim_correlated_case_data(
  n                 = 30,
  omega_B           = Omega_B_true,
  antigen_isos      = c("IgG", "IgA"),
  n_obs_per_subject = 5L
)

prepped   <- serodynamics::prep_data(sim_small)
stan_data <- prep_data_stan(prepped)

cat("Stan data structure:\n")
cat("  N =", stan_data$N, "\n")
cat("  K =", stan_data$K, "\n")
cat("  P =", stan_data$P, "\n")
cat("  max_obs =", stan_data$max_obs, "\n")
cat("  log_y dim:", paste(dim(stan_data$log_y), collapse = " x "), "\n")
cat("  antigens attr:", attr(stan_data, "antigens"), "\n")

# CRITICAL: check biomarker ordering in log_y
# log_y is [N, max_obs, K]. Which slice is IgG vs IgA?
cat("\nFirst subject, first 3 obs, both biomarkers (log_y):\n")
print(stan_data$log_y[1, 1:3, ])
cat("\nNote: column 1 should be", attr(stan_data, "antigens")[1], "\n")
cat("      column 2 should be", attr(stan_data, "antigens")[2], "\n\n")

# ==========================================================================
# LAYER 4: One small Stan fit, check Omega_B extraction
# ==========================================================================
cat("\n###############################################\n")
cat("### LAYER 4: Single Stan fit, careful extraction\n")
cat("###############################################\n\n")

# Use generous settings for this debug fit
cat("Fitting Stan (n=30, generous settings, ~10 min)...\n")
t0 <- Sys.time()

fit <- run_mod_stan(
  data            = sim_small,
  model           = "model_2",
  chains          = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  parallel_chains = 4,
  adapt_delta     = 0.99,
  max_treedepth   = 14,
  init            = 0.1,
  with_post       = TRUE,
  stan_dir        = "inst/stan",
  refresh         = 200,
  show_messages   = FALSE
)
elapsed <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("Fit done in %.1f min\n\n", elapsed))

sf <- attr(fit, "stan_fit")[[1]]

# Diagnostics
diag <- sf$diagnostic_summary()
cat("Divergent transitions:", sum(diag$num_divergent), "/ 4000\n")
cat("Max treedepth hits:", sum(diag$num_max_treedepth), "/ 4000\n\n")

# Multiple ways to extract Omega_B[1,2]
omega_B_12_med <- diagnose_omega_B(sf, TRUE_RHO)

# ==========================================================================
# LAYER 5: Also extract OTHER posterior parts — what's going on overall?
# ==========================================================================
cat("\n###############################################\n")
cat("### LAYER 5: Other posterior diagnostics\n")
cat("###############################################\n\n")

print_posterior_diagnostics(sf)

# ==========================================================================
# SUMMARY
# ==========================================================================
cat("\n========================================================\n")
cat(" DEBUG SUMMARY\n")
cat("========================================================\n")
cat(sprintf(" TRUE rho_B:       %+.3f\n", TRUE_RHO))
cat(sprintf(" Recovered:        %+.3f\n", omega_B_12_med))
cat(sprintf(" Divergent rate:   %.1f%%\n",
            100 * sum(diag$num_divergent) / 4000))

cat("\n=== DIAGNOSIS GUIDE ===\n")
cat("If recovered is near +0.6 → no bug, prior R=200 run had bad settings\n")
cat("If recovered is NEGATIVE → bug in extraction OR Stan model\n")
cat("If extreme divergence → model identifiability problem\n")
cat("========================================================\n")
