# Helper: draw subject parameters from the Kronecker covariance.
# Returns theta_arr (N x P x K).
#' @keywords internal
#' @noRd
.draw_subject_params <- function(n, mu_vec, sigma_full,
                                 n_param, n_biomarker, antigen_isos) {
  theta_vec <- MASS::mvrnorm(n = n, mu = mu_vec, Sigma = sigma_full)
  if (is.null(dim(theta_vec))) {
    theta_vec <- matrix(theta_vec, nrow = 1L)
  }

  # dim n x PK; reshape to N x P x K
  theta_arr <- array(NA_real_, dim = c(n, n_param, n_biomarker))

  for (i in seq_len(n)) {
    theta_arr[i, , ] <- matrix(
      theta_vec[i, ],
      nrow = n_param,
      ncol = n_biomarker
    )
  }

  dimnames(theta_arr) <- list(
    subject = as.character(seq_len(n)),
    param = c("log_y0", "log_y1m0", "log_t1", "log_alpha", "log_rm1"),
    biomarker = antigen_isos
  )

  theta_arr
}
