# Helper: extract Omega_B, Sigma_B, Omega_P, Sigma_P from a cmdstanr fit.
# Returns a (possibly empty) named list. Model 2 only.
#' @keywords internal
#' @noRd
.extract_kronecker_matrices_stan <- function(stan_fit, K, param_names,
                                             antigens) {
  tryCatch({
    omega_B_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Omega_B")
    )
    sigma_B_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Sigma_B")
    )
    omega_P_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Omega_P")
    )
    sigma_P_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Sigma_P")
    )

    omega_B <- summarize_matrix_draws(omega_B_arr, "Omega_B", K, K)
    sigma_B <- summarize_matrix_draws(sigma_B_arr, "Sigma_B", K, K)
    omega_P <- summarize_matrix_draws(omega_P_arr, "Omega_P", 5L, 5L)
    sigma_P <- summarize_matrix_draws(sigma_P_arr, "Sigma_P", 5L, 5L)

    dimnames(omega_B) <- list(antigens, antigens)
    dimnames(sigma_B) <- list(antigens, antigens)
    dimnames(omega_P) <- list(param_names, param_names)
    dimnames(sigma_P) <- list(param_names, param_names)

    list(Omega_B = omega_B, Sigma_B = sigma_B,
         Omega_P = omega_P, Sigma_P = sigma_P)
  }, error = function(e) {
    cli::cli_warn("Kronecker matrices not extracted: {e$message}")
    list()
  })
}
