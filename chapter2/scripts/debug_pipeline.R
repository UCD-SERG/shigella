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
  library(serodynamics)
  library(cmdstanr)
  library(posterior)
  library(tibble)
})

source("R/prep_data_stan.R")
source("R/prep_priors_stan.R")
source("R/postprocess_stan_output.R")
source("R/run_mod_stan.R")
source("R/sim_correlated_case_data.R")

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

# If sim_correlated_case_data() stores the TRUE per-subject params,
# we can compute the empirical correlation between IgG-params and IgA-params
# WITHOUT any Stan fitting.

truth <- attr(sim_dat, "truth")
if (!is.null(truth) && !is.null(truth$Theta_natural)) {
  Theta <- truth$Theta_natural   # likely [N, K, P]
  cat("Theta_natural dim:", paste(dim(Theta), collapse = " x "), "\n\n")

  if (length(dim(Theta)) == 3) {
    N_sim <- dim(Theta)[1]
    K_sim <- dim(Theta)[2]
    P_sim <- dim(Theta)[3]

    cat(sprintf("N=%d, K=%d, P=%d\n\n", N_sim, K_sim, P_sim))

    # For each of P parameters, what's the IgG-IgA correlation across subjects?
    param_names <- c("y0", "y1", "t1", "alpha", "shape")
    cat("Per-parameter empirical rho between biomarkers (should ALL be ~0.6):\n")
    for (p in 1:P_sim) {
      x <- Theta[, 1, p]
      y <- Theta[, 2, p]
      rho <- cor(x, y)
      cat(sprintf("  param %d (%s): rho = %+.3f\n",
                  p, param_names[p], rho))
    }

    # Average correlation
    rhos <- sapply(1:P_sim, function(p) cor(Theta[, 1, p], Theta[, 2, p]))
    cat(sprintf("\nMean rho across params: %+.3f (true should be %.1f)\n",
                mean(rhos), TRUE_RHO))
  }
}

# Even without truth attribute, we can check observed-data correlation
cat("\n--- Observed-data correlation (proxy) ---\n")
sim_wide <- sim_dat %>%
  select(id, antigen_iso, visit_num, value) %>%
  pivot_wider(names_from = antigen_iso, values_from = value)

if ("IgG" %in% colnames(sim_wide) && "IgA" %in% colnames(sim_wide)) {
  # Per-subject means (proxy for kinetic parameters)
  per_subj <- sim_wide %>%
    group_by(id) %>%
    summarise(
      mean_IgG = mean(log(IgG), na.rm = TRUE),
      mean_IgA = mean(log(IgA), na.rm = TRUE),
      max_IgG  = max(log(IgG), na.rm = TRUE),
      max_IgA  = max(log(IgA), na.rm = TRUE)
    )

  cat(sprintf("Cor(mean log IgG, mean log IgA): %+.3f\n",
              cor(per_subj$mean_IgG, per_subj$mean_IgA,
                  use = "complete.obs")))
  cat(sprintf("Cor(max log IgG, max log IgA):   %+.3f\n",
              cor(per_subj$max_IgG, per_subj$max_IgA,
                  use = "complete.obs")))
  cat("(These should be POSITIVE if rho_B = 0.6 is real)\n\n")
}

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
cat("=== Omega_B[1,2] extraction — multiple methods ===\n")

# Method 1: cmdstanr draws_df
m1 <- posterior::as_draws_df(sf$draws(variables = "Omega_B"))
omega_B_12 <- m1[["Omega_B[1,2]"]]
omega_B_21 <- m1[["Omega_B[2,1]"]]
cat(sprintf("Method 1 — Omega_B[1,2]: median = %+.3f, mean = %+.3f, n = %d\n",
            median(omega_B_12), mean(omega_B_12), length(omega_B_12)))
cat(sprintf("Method 1 — Omega_B[2,1]: median = %+.3f, mean = %+.3f\n",
            median(omega_B_21), mean(omega_B_21)))
# 1,2 and 2,1 should be IDENTICAL (it's a correlation matrix)

# Method 2: as_draws_array (more direct)
m2 <- posterior::as_draws_array(sf$draws(variables = "Omega_B"))
cat("\nMethod 2 — as_draws_array dims:", paste(dim(m2), collapse = " x "), "\n")
cat("Variables:", paste(dimnames(m2)$variable, collapse = ", "), "\n")

# Method 3: summary
omega_summary <- sf$summary(variables = "Omega_B")
cat("\nMethod 3 — summary:\n")
print(omega_summary)

# Compare to TRUE
cat(sprintf("\n** TRUE rho_B = %.3f **\n", TRUE_RHO))
cat(sprintf("** Recovered (method 1 median) = %+.3f **\n", median(omega_B_12)))
cat(sprintf("** Bias = %+.3f **\n", median(omega_B_12) - TRUE_RHO))

# ==========================================================================
# LAYER 5: Also extract OTHER posterior parts — what's going on overall?
# ==========================================================================
cat("\n###############################################\n")
cat("### LAYER 5: Other posterior diagnostics\n")
cat("###############################################\n\n")

# Population means M
M_summary <- sf$summary(variables = "M")
cat("=== Population means M (M[k, p]) ===\n")
print(M_summary)

# Sigma_B
SigB_summary <- sf$summary(variables = "Sigma_B")
cat("\n=== Sigma_B (covariance) ===\n")
print(SigB_summary)

# Omega_P
OmP_summary <- sf$summary(variables = "Omega_P")
cat("\n=== Omega_P (parameter correlation, top 5 rows) ===\n")
print(head(OmP_summary, 10))

# ==========================================================================
# SUMMARY
# ==========================================================================
cat("\n========================================================\n")
cat(" DEBUG SUMMARY\n")
cat("========================================================\n")
cat(sprintf(" TRUE rho_B:       %+.3f\n", TRUE_RHO))
cat(sprintf(" Recovered:        %+.3f\n", median(omega_B_12)))
cat(sprintf(" Divergent rate:   %.1f%%\n",
            100 * sum(diag$num_divergent) / 4000))

cat("\n=== DIAGNOSIS GUIDE ===\n")
cat("If recovered is near +0.6 → no bug, prior R=200 run had bad settings\n")
cat("If recovered is NEGATIVE → bug in extraction OR Stan model\n")
cat("If extreme divergence → model identifiability problem\n")
cat("========================================================\n")
