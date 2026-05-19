# Helper: build Sigma_P, Sigma_B, Sigma_eps, Sigma_full, and mu_vec.
# Returns a list with sigma_eps, sigma_full, mu_vec.
#' @keywords internal
#' @noRd
.build_sigma_matrices <- function(mu, tau_P, tau_B, tau_eps,
                                  omega_P, omega_B, omega_eps) {
  n_param <- length(mu)
  n_biomarker <- length(tau_B)

  sigma_P   <- diag(tau_P) %*% omega_P %*% diag(tau_P)
  sigma_B   <- diag(tau_B) %*% omega_B %*% diag(tau_B)
  sigma_eps <- diag(tau_eps) %*% omega_eps %*% diag(tau_eps)

  # Sigma_full = Sigma_B kron Sigma_P, dimension PK x PK
  sigma_full <- kronecker(sigma_B, sigma_P)

  # vec(M) where M is P x K (columns = biomarkers)
  # Assume same mu for all biomarkers (can be extended)
  mean_matrix <- matrix(
    mu,
    nrow = n_param,
    ncol = n_biomarker,
    byrow = FALSE
  )
  mu_vec <- as.vector(mean_matrix)

  list(sigma_eps = sigma_eps, sigma_full = sigma_full, mu_vec = mu_vec)
}
