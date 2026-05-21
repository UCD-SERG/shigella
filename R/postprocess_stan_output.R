#' @title Post-process Stan output to sr_model format (cmdstanr version)
#' @description
#' Converts a CmdStanMCMC object (from cmdstanr's mod$sample()) into the
#' long-format tibble produced by run_mod(), so downstream plotting/summary
#' functions work without modification.
#'
#' @param stan_fit CmdStanMCMC object from cmdstanr (Not rstan stanfit)
#' @param ids subject IDs from attr(stan_data, "ids")
#' @param antigens biomarker names from attr(stan_data, "antigens")
#' @param model "model_1", "model_2"
#' @param stratification label for this stratum
#' @return list with sr_tibble and cov_summaries
#' @example inst/examples/postprocess_stan_output-examples.R
#' @export
postprocess_stan_output <- function(stan_fit,
                                    ids,
                                    antigens,
                                    model = c("model_2", "model_1"),
                                    stratification = "None") {

  model <- match.arg(model)
  has_kron <- identical(model, "model_2")

  if (!requireNamespace("posterior", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg posterior} required for cmdstanr
                   postprocessing.")
  }

  # Names must match `generated quantities` block in both model_1.stan and 
  # model_2.stan.
  param_names <- c("y0", "y1", "t1", "alpha", "shape")
  N <- length(ids)
  K <- length(antigens)

  draws_df <- tibble::as_tibble(posterior::as_draws_df(
    stan_fit$draws(variables = param_names)
  ))
  n_chain <- max(draws_df$.chain)

  sr_tibble <- .extract_param_draws( # nolint: object_usage_linter
    param_names    = param_names,
    draws_df       = draws_df,
    N              = N,
    K              = K,
    ids            = ids,
    antigens       = antigens,
    stratification = stratification,
    n_chain        = n_chain
  )

  if (has_kron) {
    cov_summaries <- c(
      .extract_residual_cov_stan(stan_fit, K, antigens), # nolint: object_usage_linter
      .extract_kron_matrices(stan_fit, K, param_names, antigens) # nolint: object_usage_linter
    )
  } else {
    cov_summaries <- .extract_model1_omega_p_stan( # nolint: object_usage_linter
      stan_fit, K, param_names, antigens
    )
  }

  cov_summaries <- c(cov_summaries, .extract_log_lik_stan(stan_fit)) # nolint: object_usage_linter

  return(list(
    sr_tibble     = sr_tibble,
    cov_summaries = cov_summaries
  ))
}
