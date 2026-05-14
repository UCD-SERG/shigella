#' @title Simulate correlated longitudinal case data (Chapter 2)
#' @description
#' Extends [serodynamics::sim_case_data()] to inject known correlation
#' structure at two levels:
#'
#' 1. **Parameter-level correlation** (Omega_B): "high IgG responder
#'    tends to be high IgA responder" — implemented via the Kronecker
#'    structure vec(Theta_i) ~ N(vec(M), Sigma_B kron Sigma_P).
#' 2. **Residual-level correlation** (Omega_eps): "IgG and IgA
#'    measurement errors co-vary within a time point" — implemented via
#'    multivariate log-normal observation noise with covariance
#'    Sigma_eps.
#'
#' This is the data-generating process for the Chapter 2 simulation
#' study.
#'
#' @param n [integer] number of individuals to simulate
#' @param mu [numeric] length-P vector of population means on log scale
#'   (defaults match JAGS prep_priors)
#' @param tau_P [numeric] length-P vector of SDs across kinetic
#'   parameters
#' @param tau_B [numeric] length-K vector of SDs across biomarkers
#' @param tau_eps [numeric] length-K vector of residual SDs
#' @param Omega_P [matrix] P x P parameter correlation matrix
#'   (default: identity — no within-biomarker parameter correlation)
#' @param Omega_B [matrix] K x K biomarker correlation matrix
#'   (default: identity — Scenario 2, residual correlation only)
#' @param Omega_eps [matrix] K x K residual correlation matrix
#'   (default: identity — no residual correlation)
#' @param antigen_isos [character] names for the K biomarkers
#' @param n_obs_per_subject [integer] number of observations per subject
#'   (default 5, matching the Shigella SOSAR cohort)
#' @param time_grid [numeric] follow-up times in days
#'   (default c(2, 7, 30, 90, 180) mimicking Chapter 1)
#' @param seed [integer] RNG seed
#'
#' @returns a `case_data` object plus attributes recording the truth:
#'   - `"truth"` — list with mu, tau_P, tau_B, tau_eps, Omega_P,
#'     Omega_B, Omega_eps
#'   - `"theta_true"` — N x K x P array of true subject parameters
#' @export
#' @example inst/examples/sim_correlated_case_data-examples.R
sim_correlated_case_data <- function(
    n               = 48,
    mu              = c(1.0, 7.0, 1.0, -4.0, -1.0),
    tau_P           = c(0.5, 0.7, 0.3, 1.0, 0.4),
    tau_B           = c(0.8, 0.8),
    tau_eps         = c(0.3, 0.3),
    Omega_P         = diag(5),
    Omega_B         = diag(2),
    Omega_eps       = diag(2),
    antigen_isos    = c("biomarker_1", "biomarker_2"),
    n_obs_per_subject = 5L,
    time_grid       = c(2, 7, 30, 90, 180),
    seed            = NULL) {

  if (!is.null(seed)) set.seed(seed)

  P <- length(mu)
  K <- length(antigen_isos)

  if (K != length(tau_B))    cli::cli_abort("length(tau_B) must equal K")
  if (K != length(tau_eps))  cli::cli_abort("length(tau_eps) must equal K")
  if (P != length(tau_P))    cli::cli_abort("length(tau_P) must equal P (5)")
  if (any(dim(Omega_P)   != c(P, P))) cli::cli_abort("Omega_P must be PxP")
  if (any(dim(Omega_B)   != c(K, K))) cli::cli_abort("Omega_B must be KxK")
  if (any(dim(Omega_eps) != c(K, K))) cli::cli_abort("Omega_eps must be KxK")
  if (length(time_grid) < n_obs_per_subject) {
    cli::cli_abort("time_grid must have at least n_obs_per_subject entries")
  }

  # --- Build Kronecker covariance on parameters ---
  Sigma_P   <- diag(tau_P)   %*% Omega_P   %*% diag(tau_P)
  Sigma_B   <- diag(tau_B)   %*% Omega_B   %*% diag(tau_B)
  Sigma_eps <- diag(tau_eps) %*% Omega_eps %*% diag(tau_eps)

  # Σ_full = Σ_B ⊗ Σ_P, dimension PK × PK
  Sigma_full <- kronecker(Sigma_B, Sigma_P)

  # vec(M) where M is P × K (columns = biomarkers)
  # Assume same mu for all biomarkers (can be extended)
  M <- matrix(mu, nrow = P, ncol = K, byrow = FALSE)
  mu_vec <- as.vector(M)   # column-major stack

  # Draw θ_i for each subject
  theta_vec <- MASS::mvrnorm(n = n, mu = mu_vec, Sigma = Sigma_full)
  # dim n × PK; reshape to N × P × K
  theta_arr <- array(NA, dim = c(n, P, K))
  for (i in seq_len(n)) {
    theta_arr[i, , ] <- matrix(theta_vec[i, ], nrow = P, ncol = K)
  }
  dimnames(theta_arr) <- list(
    subject = as.character(seq_len(n)),
    param   = c("log_y0", "log_y1m0", "log_t1", "log_alpha", "log_rm1"),
    biomarker = antigen_isos
  )

  # --- Generate observations ---
  L_eps <- chol(Sigma_eps)  # upper triangular; use t() for lower

  rows <- list()
  row_counter <- 1L
  for (i in seq_len(n)) {
    obs_times <- sort(sample(time_grid, size = n_obs_per_subject, replace = FALSE))
    for (tt_idx in seq_along(obs_times)) {
      tt <- obs_times[tt_idx]

      # Compute log mu for each biomarker
      log_mu_k <- numeric(K)
      for (j in seq_len(K)) {
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
            log_mu_k[j] <- log(y0)  # floor
          } else {
            log_mu_k[j] <- log(term) / (1 - shape)
          }
        }
      }

      # Add correlated residual noise
      z <- rnorm(K)
      log_y_obs <- log_mu_k + as.vector(t(L_eps) %*% z)

      for (j in seq_len(K)) {
        rows[[row_counter]] <- data.frame(
          id          = as.character(i),
          visit_num   = tt_idx,
          timeindays  = tt,
          antigen_iso = antigen_isos[j],
          value       = exp(log_y_obs[j]),
          stringsAsFactors = FALSE
        )
        row_counter <- row_counter + 1L
      }
    }
  }

  sim_df <- dplyr::bind_rows(rows)

  # Convert to case_data
  case <- sim_df |>
    serodynamics::as_case_data(
      id_var        = "id",
      biomarker_var = "antigen_iso",
      time_in_days  = "timeindays",
      value_var     = "value"
    )

  # Attach ground truth
  attr(case, "truth") <- list(
    mu        = mu,
    tau_P     = tau_P,
    tau_B     = tau_B,
    tau_eps   = tau_eps,
    Omega_P   = Omega_P,
    Omega_B   = Omega_B,
    Omega_eps = Omega_eps,
    Sigma_P   = Sigma_P,
    Sigma_B   = Sigma_B,
    Sigma_eps = Sigma_eps
  )
  attr(case, "theta_true") <- theta_arr

  return(case)
}
