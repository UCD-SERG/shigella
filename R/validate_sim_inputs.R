# Helper: validate simulation inputs.
#' @keywords internal
#' @noRd
.validate_sim_inputs <- function(n_param, n_biomarker,
                                 tau_B, tau_eps, tau_P,
                                 omega_P, omega_B, omega_eps,
                                 time_grid, n_obs_per_subject) {
  if (n_biomarker != length(tau_B)) {
    cli::cli_abort("{.arg tau_B} must have length K.")
  }

  if (n_biomarker != length(tau_eps)) {
    cli::cli_abort("{.arg tau_eps} must have length K.")
  }

  if (n_param != length(tau_P)) {
    cli::cli_abort("{.arg tau_P} must have length P.")
  }

  if (!is.matrix(omega_P) || !identical(dim(omega_P), c(n_param, n_param))) {
    cli::cli_abort("{.arg omega_P} must be a P x P matrix.")
  }

  if (!is.matrix(omega_B) || !identical(dim(omega_B), c(n_biomarker, n_biomarker))) {
    cli::cli_abort("{.arg omega_B} must be a K x K matrix.")
  }

  if (!is.matrix(omega_eps) || !identical(dim(omega_eps), c(n_biomarker, n_biomarker))) {
    cli::cli_abort("{.arg omega_eps} must be a K x K matrix.")
  }

  if (length(time_grid) < n_obs_per_subject) {
    cli::cli_abort(
      "{.arg time_grid} must have at least {.arg n_obs_per_subject} entries."
    )
  }
}
