# Diagnostic helper functions for debug_pipeline.R

#' Check per-parameter empirical rho and observed-data correlation proxy
#'
#' Prints Layer 2 diagnostics: per-parameter recovery from truth attribute
#' and observed-data correlation as a proxy check.
#'
#' @param sim_dat case_data from sim_correlated_case_data()
#' @param TRUE_RHO the true biomarker correlation used in simulation
diagnose_layer2 <- function(sim_dat, TRUE_RHO) {
  truth <- attr(sim_dat, "truth")
  param_names <- c("y0", "y1", "t1", "alpha", "shape")

  if (!is.null(truth) && !is.null(attr(sim_dat, "theta_true"))) {
    Theta <- attr(sim_dat, "theta_true")
    cat("theta_true dim:", paste(dim(Theta), collapse = " x "), "\n\n")

    if (length(dim(Theta)) == 3) {
      N_sim <- dim(Theta)[1]
      P_sim <- dim(Theta)[3]

      cat(sprintf(
        "Per-parameter empirical rho between biomarkers (should ALL be ~%.1f):\n",
        TRUE_RHO
      ))
      rhos <- numeric(P_sim)
      for (p in seq_len(P_sim)) {
        x <- Theta[, 1, p]
        y <- Theta[, 2, p]
        rhos[p] <- cor(x, y)
        cat(sprintf("  param %d (%s): rho = %+.3f\n",
                    p, param_names[p], rhos[p]))
      }
      cat(sprintf("\nMean rho across params: %+.3f (true should be %.1f)\n",
                  mean(rhos), TRUE_RHO))
    }
  }

  cat("\n--- Observed-data correlation (proxy) ---\n")
  sim_wide <- sim_dat |>
    dplyr::select("id", "antigen_iso", "visit_num", "value") |>
    tidyr::pivot_wider(names_from = "antigen_iso", values_from = "value")

  if ("IgG" %in% colnames(sim_wide) && "IgA" %in% colnames(sim_wide)) {
    per_subj <- sim_wide |>
      dplyr::group_by(id) |>
      dplyr::summarise(
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
}


#' Diagnose Omega_B recovery using multiple extraction methods
#'
#' @param sf CmdStanMCMC object (raw fit from cmdstanr)
#' @param TRUE_RHO the true rho_B used in simulation
#' @return (invisible) the median of Omega_B[1,2] from Method 1
diagnose_omega_B <- function(sf, TRUE_RHO) {
  cat("=== Omega_B[1,2] extraction — multiple methods ===\n")

  m1 <- posterior::as_draws_df(sf$draws(variables = "Omega_B"))
  omega_B_12 <- m1[["Omega_B[1,2]"]]
  omega_B_21 <- m1[["Omega_B[2,1]"]]
  cat(sprintf(
    "Method 1 — Omega_B[1,2]: median = %+.3f, mean = %+.3f, n = %d\n",
    median(omega_B_12), mean(omega_B_12), length(omega_B_12)
  ))
  cat(sprintf("Method 1 — Omega_B[2,1]: median = %+.3f, mean = %+.3f\n",
              median(omega_B_21), mean(omega_B_21)))

  m2 <- posterior::as_draws_array(sf$draws(variables = "Omega_B"))
  cat("\nMethod 2 — as_draws_array dims:",
      paste(dim(m2), collapse = " x "), "\n")
  cat("Variables:", paste(dimnames(m2)$variable, collapse = ", "), "\n")

  cat("\nMethod 3 — summary:\n")
  print(sf$summary(variables = "Omega_B"))

  cat(sprintf("\n** TRUE rho_B = %.3f **\n", TRUE_RHO))
  cat(sprintf("** Recovered (method 1 median) = %+.3f **\n",
              median(omega_B_12)))
  cat(sprintf("** Bias = %+.3f **\n", median(omega_B_12) - TRUE_RHO))

  invisible(median(omega_B_12))
}


#' Print other posterior diagnostics (M, Sigma_B, Omega_P)
#'
#' @param sf CmdStanMCMC object
print_posterior_diagnostics <- function(sf) {
  cat("=== Population means M (M[k, p]) ===\n")
  print(sf$summary(variables = "M"))

  cat("\n=== Sigma_B (covariance) ===\n")
  print(sf$summary(variables = "Sigma_B"))

  cat("\n=== Omega_P (parameter correlation, top 10 rows) ===\n")
  print(head(sf$summary(variables = "Omega_P"), 10))
}
