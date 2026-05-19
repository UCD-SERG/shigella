# Helper: compute log mu for each biomarker for subject i at time tt.
# Returns log_mu_k (length-K numeric vector).
#' @keywords internal
#' @noRd
.compute_log_mu_k <- function(theta_arr, i, n_biomarker, tt) {
  log_mu_k <- numeric(n_biomarker)

  for (j in seq_len(n_biomarker)) {
    log_y0    <- theta_arr[i, 1, j]
    log_y1m0  <- theta_arr[i, 2, j]
    log_t1    <- theta_arr[i, 3, j]
    log_alpha <- theta_arr[i, 4, j]
    log_rm1   <- theta_arr[i, 5, j]

    y0    <- exp(log_y0)
    y1    <- y0 + exp(log_y1m0)
    t1_j  <- exp(log_t1)
    alpha <- exp(log_alpha)
    shape <- exp(log_rm1) + 1

    if (tt <= t1_j) {
      beta_growth <- (log(y1) - log(y0)) / t1_j
      log_mu_k[j] <- log(y0) + beta_growth * tt
    } else {
      term <- y1^(1 - shape) - (1 - shape) * alpha * (tt - t1_j)

      if (term <= 0) {
        log_mu_k[j] <- log(y0)
      } else {
        log_mu_k[j] <- log(term) / (1 - shape)
      }
    }
  }

  log_mu_k
}
