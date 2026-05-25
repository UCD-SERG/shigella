# Helper: compute log mu for each biomarker for subject i at time tt.
# Returns log_mu_k (length-K numeric vector).
#' @keywords internal
#' @noRd
.compute_log_mu_k <- function(theta_arr, i, n_biomarker, tt) {
  vapply(
    seq_len(n_biomarker),
    function(j) {
      .compute_kinetics_at_time( # nolint: object_usage_linter
        theta_arr[i, 1, j],
        theta_arr[i, 2, j],
        theta_arr[i, 3, j],
        theta_arr[i, 4, j],
        theta_arr[i, 5, j],
        tt
      )
    },
    numeric(1L)
  )
}
