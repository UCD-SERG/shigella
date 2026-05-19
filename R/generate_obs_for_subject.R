# Helper: generate all observation rows for one subject.
# Returns a flat list of data frames (one per biomarker x time point).
#' @keywords internal
#' @noRd
.generate_obs_for_subject <- function(i, n_obs_per_subject, time_grid,
                                      n_biomarker, theta_arr, antigen_isos,
                                      l_eps) {
  obs_times <- sort(
    sample(time_grid, size = n_obs_per_subject, replace = FALSE)
  )

  rows <- list()
  row_counter <- 1L

  for (tt_idx in seq_along(obs_times)) {
    tt        <- obs_times[tt_idx]
    log_mu_k  <- .compute_log_mu_k(theta_arr, i, n_biomarker, tt)
    z         <- rnorm(n_biomarker)
    log_y_obs <- log_mu_k + as.vector(t(l_eps) %*% z)

    for (j in seq_len(n_biomarker)) {
      rows[[row_counter]] <- data.frame(
        id           = as.character(i),
        visit_num    = tt_idx,
        timeindays   = tt,
        antigen_iso  = antigen_isos[j],
        value        = exp(log_y_obs[j]),
        stringsAsFactors = FALSE
      )
      row_counter <- row_counter + 1L
    }
  }

  rows
}
