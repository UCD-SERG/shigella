# Helper: extract Omega_eps and Sigma_eps from a cmdstanr fit.
# Returns a (possibly empty) named list. Model 2 only.
#' @keywords internal
#' @noRd
.extract_residual_cov_stan <- function(stan_fit, K, antigens) {
  tryCatch({
    omega_eps_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Omega_eps")
    )
    sigma_eps_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Sigma_eps")
    )
    omega_eps <- .summarize_matrix_draws(omega_eps_arr, "Omega_eps", K, K)
    sigma_eps <- .summarize_matrix_draws(sigma_eps_arr, "Sigma_eps", K, K)
    dimnames(omega_eps) <- list(antigens, antigens)
    dimnames(sigma_eps) <- list(antigens, antigens)
    list(Omega_eps = omega_eps, Sigma_eps = sigma_eps)
  }, error = function(e) {
    cli::cli_warn("Omega_eps/Sigma_eps not extracted: {e$message}")
    list()
  })
}
