# ==========================================================================
# inspect_sim_function.R — Look at what sim_correlated_case_data() does
#
# Goal: Confirm the simulation function correctly implements the Kronecker
# structure with the requested rho_B.
#
# Run:
#   cd ~/chapter2
#   Rscript scripts/inspect_sim_function.R 2>&1 | tee logs/inspect_sim.log
# ==========================================================================

setwd("~/chapter2")

cat("\n========================================================\n")
cat(" INSPECT: sim_correlated_case_data() internals\n")
cat("========================================================\n\n")

library(shigella)

# Print the function body
cat("=== Function body ===\n")
print(sim_correlated_case_data)
cat("\n")

# Inspect arguments
cat("=== Arguments + defaults ===\n")
formals(sim_correlated_case_data)
cat("\n")

# ==========================================================================
# Run with explicit rho_B = 0.6 and large n
# ==========================================================================
cat("=== Generate n=500 with rho_B = 0.6 ===\n\n")
set.seed(123)

sim_dat <- sim_correlated_case_data(
  n                 = 500,
  omega_B           = matrix(c(1, 0.6, 0.6, 1), 2, 2),
  antigen_isos      = c("IgG", "IgA"),
  n_obs_per_subject = 5L
)

cat("Generated rows:", nrow(sim_dat), "\n")
cat("Unique subjects:", length(unique(sim_dat$id)), "\n\n")

# ==========================================================================
# Inspect the truth attribute if it exists
# ==========================================================================
cat("=== attr(sim_dat, 'truth') ===\n")
truth <- attr(sim_dat, "truth")
if (is.null(truth)) {
  cat("NO truth attribute. Will compute empirical correlation from data.\n\n")
} else {
  cat("Truth attribute keys:", paste(names(truth), collapse = ", "), "\n\n")
  for (k in names(truth)) {
    x <- truth[[k]]
    cat(sprintf("  %s: class=%s, ", k, class(x)[1]))
    if (is.null(dim(x))) {
      cat("length=", length(x), "\n")
    } else {
      cat("dim=", paste(dim(x), collapse = " x "), "\n")
    }
  }
  cat("\n")
}

# ==========================================================================
# Check theta (or Theta) if present
# ==========================================================================
if (!is.null(truth) && (!is.null(truth$Theta) || !is.null(truth$theta))) {
  Theta <- truth$Theta %||% truth$theta
  cat("=== True Theta inspection ===\n")
  cat("Dim:", paste(dim(Theta), collapse = " x "), "\n")

  if (length(dim(Theta)) == 3) {
    # Try [N, K, P] interpretation
    cat("\nInterpretation 1: [N, K, P] with K=biomarkers, P=params\n")
    cat("Per-parameter biomarker correlation (should be ~0.6):\n")
    for (p in 1:dim(Theta)[3]) {
      r <- cor(Theta[, 1, p], Theta[, 2, p])
      cat(sprintf("  param %d: cor(B1, B2) = %+.3f\n", p, r))
    }

    # Try [N, P, K] interpretation
    cat("\nInterpretation 2: [N, P, K]\n")
    cat("Per-biomarker cross-param correlation (should be ~0.6 if [N,P,K]):\n")
    for (k in 1:dim(Theta)[3]) {
      mat <- Theta[, , k]
      cat(sprintf("  bmk %d: cor(P1, P2) = %+.3f\n", k, cor(mat[, 1], mat[, 2])))
    }
  }
}

# ==========================================================================
# Just check raw data — make sure positive correlation is REALLY there
# ==========================================================================
cat("\n=== Raw observed data correlation check ===\n\n")

library(dplyr)
library(tidyr)

# Subject-level summaries — should reflect underlying correlation
subj_summary <- sim_dat %>%
  select(id, antigen_iso, visit_num, value) %>%
  pivot_wider(names_from = antigen_iso, values_from = value)

if (all(c("IgG", "IgA") %in% colnames(subj_summary))) {
  per_subj <- subj_summary %>%
    group_by(id) %>%
    summarise(
      log_IgG_max  = max(log(IgG + 1), na.rm = TRUE),
      log_IgA_max  = max(log(IgA + 1), na.rm = TRUE),
      log_IgG_mean = mean(log(IgG + 1), na.rm = TRUE),
      log_IgA_mean = mean(log(IgA + 1), na.rm = TRUE),
      .groups = "drop"
    )

  cat(sprintf("n_subjects analyzed: %d\n", nrow(per_subj)))
  cat(sprintf("Cor(log_IgG_max,  log_IgA_max):  %+.3f\n",
              cor(per_subj$log_IgG_max,  per_subj$log_IgA_max,
                  use = "complete.obs")))
  cat(sprintf("Cor(log_IgG_mean, log_IgA_mean): %+.3f\n",
              cor(per_subj$log_IgG_mean, per_subj$log_IgA_mean,
                  use = "complete.obs")))
  cat("\nINTERPRETATION:\n")
  cat("  - If rho_B = 0.6 is properly implemented in sim function,\n")
  cat("    these correlations should be POSITIVE (likely 0.3-0.7).\n")
  cat("  - If they are near 0 or NEGATIVE, the sim function is buggy.\n\n")
}

cat("========================================================\n")
cat(" DONE\n")
cat("========================================================\n")
