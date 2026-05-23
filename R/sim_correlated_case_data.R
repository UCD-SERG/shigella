#' @title Simulate correlated longitudinal case data
#' @description
#' Extends [serodynamics::sim_case_data()] to inject known correlation
#' structure at two levels:
#'
#' 1. **Parameter-level (between-biomarker) correlation** via a Kronecker
#'    covariance on the vectorised per-subject parameter matrix.
#' 2. **Residual-level correlation** via multivariate log-normal
#'    observation noise.
#'
#' **Data-generating process**
#'
#' *Subject parameters.*  Let \eqn{\Theta_i} be the \eqn{P \times K}
#' matrix of log-scale kinetic parameters for subject \eqn{i}
#' (rows = parameters, columns = biomarkers).  Draw
#' \deqn{
#'   \mathrm{vec}(\Theta_i) \sim
#'   \mathcal{N}\!\bigl(\mathrm{vec}(M),\;
#'                        \Sigma_B \otimes \Sigma_P\bigr),
#' }
#' where \eqn{M} is a \eqn{P \times K} population-mean matrix,
#' \eqn{\Sigma_P = \mathrm{diag}(\tau_P)\,\Omega_P\,\mathrm{diag}(\tau_P)}
#' is the \eqn{P \times P} within-biomarker parameter covariance, and
#' \eqn{\Sigma_B = \mathrm{diag}(\tau_B)\,\Omega_B\,\mathrm{diag}(\tau_B)}
#' is the \eqn{K \times K} between-biomarker covariance.
#'
#' *Observation model.*  For each subject \eqn{i}, time \eqn{t}, and
#' biomarker \eqn{k},
#' \deqn{
#'   \log y_{i,t,k} = \log \mu_{i,t,k} + \varepsilon_{i,t,k},
#'   \quad
#'   \boldsymbol{\varepsilon}_{i,t} \sim
#'   \mathcal{N}(\mathbf{0},\, \Sigma_\varepsilon),
#' }
#' where
#' \eqn{\Sigma_\varepsilon =
#'   \mathrm{diag}(\tau_\varepsilon)\,\Omega_\varepsilon\,
#'   \mathrm{diag}(\tau_\varepsilon)}.
#'
#' *Two-phase kinetics.*  The deterministic log-mean follows
#' \deqn{
#'   \log \mu_{i,t,k} =
#'   \begin{cases}
#'     \log(y0_{ik}) + \beta_{ik}\,t, & t \le t1_{ik},\\[4pt]
#'     \dfrac{1}{1-s_{ik}}
#'       \log\!\bigl(y1_{ik}^{1-s_{ik}}
#'              - (1-s_{ik})\,\alpha_{ik}(t - t1_{ik})\bigr),
#'     & t > t1_{ik},
#'   \end{cases}
#' }
#' with growth rate
#' \eqn{\beta_{ik} = (\log y1_{ik} - \log y0_{ik})\,/\,t1_{ik}}
#' and shape \eqn{s_{ik} = \exp(\mathtt{log\_rm1}_{ik}) + 1 > 1}.
#'
#' This is the data-generating process for a Kronecker-correlated simulation 
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
#'   (default: identity - no within-biomarker parameter correlation)
#' @param omega_B [matrix] K x K biomarker correlation matrix
#'   (default: identity - Scenario 2, residual correlation only)
#' @param omega_eps [matrix] K x K residual correlation matrix
#'   (default: identity - no residual correlation)
#' @param antigen_isos [character] names for the K biomarkers
#' @param n_obs_per_subject [integer] number of observations per subject
#'   (default 5, matching the Shigella SOSAR cohort)
#' @param time_grid [numeric] follow-up times in days
#'   (default c(2, 7, 30, 90, 180))
#' @param seed [integer] RNG seed
#'
#' @returns a `case_data` object plus attributes recording the truth:
#'   - `"truth"` - list with mu, tau_P, tau_B, tau_eps, omega_P,
#'     omega_B, omega_eps, sigma_P, sigma_B, sigma_eps
#'   - `"theta_true"` - N x P x K array of true subject parameters in
#'     **Stan's internal log-scale parameterisation**:
#'     `log_y0`, `log_y1m0`, `log_t1`, `log_alpha`, `log_rm1`.
#'     Fitted posterior summaries use natural-scale names
#'     (`y0`, `y1`, `t1`, `alpha`, `shape`); apply `exp()` before comparing.
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

  .validate_sim_inputs(
    n_param, n_biomarker,
    tau_B, tau_eps, tau_P,
    omega_P, omega_B, omega_eps,
    time_grid, n_obs_per_subject
  )

  mats <- .build_sigma_matrices(mu, tau_P, tau_B, tau_eps,
                                omega_P, omega_B, omega_eps)
  sigma_eps <- mats$sigma_eps
  sigma_full <- mats$sigma_full
  mu_vec <- mats$mu_vec

  theta_arr <- .draw_subject_params(n, mu_vec, sigma_full,
                                    n_param, n_biomarker, antigen_isos)

  # --- Generate observations ---
  l_eps <- chol(sigma_eps)

  rows <- .generate_obs_rows(n, n_obs_per_subject, time_grid,
                             n_biomarker, theta_arr, antigen_isos, l_eps)

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
    sigma_P = diag(tau_P) %*% omega_P %*% diag(tau_P),
    sigma_B = diag(tau_B) %*% omega_B %*% diag(tau_B),
    sigma_eps = sigma_eps
  )
  attr(case, "theta_true") <- theta_arr

  return(case)
}
