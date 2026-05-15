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
#' @param omega_P [matrix] P x P parameter correlation matrix
#'   (default: identity — no within-biomarker parameter correlation)
#' @param omega_B [matrix] K x K biomarker correlation matrix
#'   (default: identity — Scenario 2, residual correlation only)
#' @param omega_eps [matrix] K x K residual correlation matrix
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
#'   - `"theta_true"` — N x P x K array of true subject parameters
#' @export
#' @example inst/examples/sim_correlated_case_data-examples.R
sim_correlated_case_data <- function(
    n = 48,
    mu = c(1.0, 7.0, 1.0, -4.0, -1.0),
    tau_P = c(0.5, 0.7, 0.3, 1.0, 0.4),
    tau_B = c(0.8, 0.8),
    tau_eps = c(0.3, 0.3),
    omega_P = diag(5),
    omega_B = diag(2),
    omega_eps = diag(2),
    antigen_isos = c("biomarker_1", "biomarker_2"),
    n_obs_per_subject = 5L,
    time_grid = c(2, 7, 30, 90, 180),
    seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_param <- length(mu)
  n_biomarker <- length(antigen_isos)
  
  if (n_biomarker != length(tau_B)) {
    cli::cli_abort("{.arg tau_B} must have length K.")
  }
  
  if (n_biomarker != length(tau_eps)) {
    cli::cli_abort("{.arg tau_eps} must have length K.")
  }
  
  if (n_param != length(tau_P)) {
    cli::cli_abort("{.arg tau_P} must have length P.")
  }
  
  if (any(dim(omega_P) != c(n_param, n_param))) {
    cli::cli_abort("{.arg omega_P} must be a P x P matrix.")
  }
  
  if (any(dim(omega_B) != c(n_biomarker, n_biomarker))) {
    cli::cli_abort("{.arg omega_B} must be a K x K matrix.")
  }
  
  if (any(dim(omega_eps) != c(n_biomarker, n_biomarker))) {
    cli::cli_abort("{.arg omega_eps} must be a K x K matrix.")
  }
  
  if (length(time_grid) < n_obs_per_subject) {
    cli::cli_abort(
      "{.arg time_grid} must have at least {.arg n_obs_per_subject} entries."
    )
  }
  
  # --- Build Kronecker covariance on parameters ---
  sigma_P <- diag(tau_P) %*% omega_P %*% diag(tau_P)
  sigma_B <- diag(tau_B) %*% omega_B %*% diag(tau_B)
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
  
  # Draw theta_i for each subject
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
  
  # --- Generate observations ---
  l_eps <- chol(sigma_eps)
  
  rows <- list()
  row_counter <- 1L
  
  for (i in seq_len(n)) {
    obs_times <- sort(
      sample(time_grid, size = n_obs_per_subject, replace = FALSE)
    )
    
    for (tt_idx in seq_along(obs_times)) {
      tt <- obs_times[tt_idx]
      
      # Compute log mu for each biomarker
      log_mu_k <- numeric(n_biomarker)
      
      for (j in seq_len(n_biomarker)) {
        log_y0 <- theta_arr[i, 1, j]
        log_y1m0 <- theta_arr[i, 2, j]
        log_t1 <- theta_arr[i, 3, j]
        log_alpha <- theta_arr[i, 4, j]
        log_rm1 <- theta_arr[i, 5, j]
        
        y0 <- exp(log_y0)
        y1 <- y0 + exp(log_y1m0)
        t1_j <- exp(log_t1)
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
      
      # Add correlated residual noise
      z <- rnorm(n_biomarker)
      log_y_obs <- log_mu_k + as.vector(t(l_eps) %*% z)
      
      for (j in seq_len(n_biomarker)) {
        rows[[row_counter]] <- data.frame(
          id = as.character(i),
          visit_num = tt_idx,
          timeindays = tt,
          antigen_iso = antigen_isos[j],
          value = exp(log_y_obs[j]),
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
      id_var = "id",
      biomarker_var = "antigen_iso",
      time_in_days = "timeindays",
      value_var = "value"
    )
  
  # Attach ground truth
  attr(case, "truth") <- list(
    mu = mu,
    tau_P = tau_P,
    tau_B = tau_B,
    tau_eps = tau_eps,
    omega_P = omega_P,
    omega_B = omega_B,
    omega_eps = omega_eps,
    sigma_P = sigma_P,
    sigma_B = sigma_B,
    sigma_eps = sigma_eps
  )
  attr(case, "theta_true") <- theta_arr
  
  return(case)
}
