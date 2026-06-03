# Helper: sanity-check array sizes and zero-observation subjects.
#' @keywords internal
#' @noRd
.validate_stan_arrays <- function(nsmpl, max_obs) {
  if (any(nsmpl > max_obs)) {
    cli::cli_abort(
      "n_obs[{which(nsmpl > max_obs)}] > max_obs. Array sizes inconsistent."
    )
  }
  if (any(nsmpl == 0)) {
    cli::cli_warn(
      "Subject(s) with 0 observations detected; these contribute no likelihood."
    )
  }
}
