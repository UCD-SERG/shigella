# Helper: extract the log_lik draws matrix from a cmdstanr fit.
# Returns a (possibly empty) named list.
#' @keywords internal
#' @noRd
.extract_log_lik_stan <- function(stan_fit) {
  tryCatch({
    list(log_lik = posterior::as_draws_matrix(
      stan_fit$draws(variables = "log_lik")
    ))
  }, error = function(e) {
    cli::cli_warn("log_lik not extracted: {e$message}")
    list()
  })
}
