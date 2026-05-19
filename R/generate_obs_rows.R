# Helper: generate all observation rows for all subjects.
# Returns a list of data frames (to be dplyr::bind_rows'd).
#' @keywords internal
#' @noRd
.generate_obs_rows <- function(n, n_obs_per_subject, time_grid,
                               n_biomarker, theta_arr, antigen_isos, l_eps) {
  rows <- list()
  for (i in seq_len(n)) {
    rows <- c(rows, .generate_obs_for_subject(
      i, n_obs_per_subject, time_grid,
      n_biomarker, theta_arr, antigen_isos, l_eps
    ))
  }
  rows
}
