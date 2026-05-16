#' @title Prepare priors for Stan backend
#' @description
#' Translates the JAGS prior specification into Stan's LKJ + half-Cauchy
#' decomposition.
#'
#' Defaults match the JAGS Chapter 1 model (which works), with two
#' adjustments for Stan compatibility:
#'   - mu_hyp_sd capped at 10 (was 316 in JAGS) — Stan's HMC sampler
#'     handles weakly-informative priors better with more reasonable scales.
#'     JAGS Gibbs sampling tolerates wider priors, but Stan's gradient-based
#'     HMC explores the tails too aggressively when sd is huge.
#'   - tau scales = 1.0 (was 2.5) — keeps initial steps reasonable
#'
#'
#' @param mu_hyp_mean [numeric] length-5 prior mean for population params
#' @param mu_hyp_sd [numeric] length-5 prior SD for population params.
#'   Weakly informative, Stan-friendly. JAGS original was c(1, 316, 1, 32, 1).
#'   5.0 on log-scale params covers ~5 orders of magnitude — plenty wide.
#' @param tau_P_scale half-Cauchy scale for parameter SDs
#' @param tau_B_scale half-Cauchy scale for biomarker SDs (Model 2 only)
#' @param tau_eps_scale half-Cauchy scale for residual SDs
#' @param lkj_P_eta LKJ shape for parameter correlation
#' @param lkj_B_eta LKJ shape for biomarker correlation (Model 2 only)
#' @param lkj_eps_eta LKJ shape for residual correlation
#' @param model character: "model_1", "model_2"
#'
#' @returns named list with priors for the Stan data block
#' @example inst/examples/prep_priors_stan-examples.R
#' @export
prep_priors_stan <- function(
    mu_hyp_mean   = c(1.0, 7.0, 1.0, -4.0, -1.0),
    mu_hyp_sd     = c(5.0, 5.0, 5.0, 5.0, 5.0),
    tau_P_scale   = 1.0,
    tau_B_scale   = 1.0,
    tau_eps_scale = 1.0,
    lkj_P_eta     = 2.0,
    lkj_B_eta     = 1.0,
    lkj_eps_eta   = 2.0,
    model         = c("model_2", "model_1")) {

  model <- match.arg(model)

  has_kron <- identical(model, "model_2")

  if (length(mu_hyp_mean) != 5) {
    cli::cli_abort("{.arg mu_hyp_mean} must be length 5.")
  }
  if (length(mu_hyp_sd) != 5) {
    cli::cli_abort("{.arg mu_hyp_sd} must be length 5.")
  }

  priors <- list(
    mu_hyp_mean   = mu_hyp_mean,
    mu_hyp_sd     = mu_hyp_sd,
    tau_P_scale   = tau_P_scale,
    tau_eps_scale = tau_eps_scale,
    lkj_P_eta     = lkj_P_eta,
    lkj_eps_eta   = lkj_eps_eta
  )

  if (has_kron) {
    priors$tau_B_scale <- tau_B_scale
    priors$lkj_B_eta   <- lkj_B_eta
  }

  attr(priors, "model")            <- model
  attr(priors, "used_stan_priors") <- priors
  class(priors) <- c("curve_params_priors_stan", "list")
  return(priors)
}
