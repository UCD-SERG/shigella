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

  param_names <- c("y0", "y1", "t1", "alpha", "shape")
  N <- length(ids)
  K <- length(antigens)

  # cmdstanr returns draws via $draws() which is a draws_array
  # Convert to data frame format for processing
  draws_df <- tibble::as_tibble(posterior::as_draws_df(
    stan_fit$draws(variables = param_names)
  ))

  n_chain <- max(draws_df$.chain)

  sr_tibble <- .extract_param_draws(
    param_names    = param_names,
    draws_df       = draws_df,
    N              = N,
    K              = K,
    ids            = ids,
    antigens       = antigens,
    stratification = stratification,
    n_chain        = n_chain
  )

  cov_summaries <- list()

  # ---- Residual covariance (model_2 only — model_1 uses independent residuals)
  if (has_kron) {
    tryCatch({
      omega_eps_arr <- posterior::as_draws_array(
        stan_fit$draws(variables = "Omega_eps")
      )
      sigma_eps_arr <- posterior::as_draws_array(
        stan_fit$draws(variables = "Sigma_eps")
      )
      # Compute median across iterations and chains for each cell
      omega_eps_mat <- summarize_matrix_draws(omega_eps_arr, "Omega_eps", K, K)
      sigma_eps_mat <- summarize_matrix_draws(sigma_eps_arr, "Sigma_eps", K, K)
      cov_summaries$Omega_eps <- omega_eps_mat
      cov_summaries$Sigma_eps <- sigma_eps_mat
      dimnames(cov_summaries$Omega_eps) <- list(antigens, antigens)
      dimnames(cov_summaries$Sigma_eps) <- list(antigens, antigens)
    }, error = function(e) {
      cli::cli_warn("Omega_eps/Sigma_eps not extracted: {e$message}")
    })
  }

  # ---- Kronecker matrices (Model 2 only) ----
  if (has_kron) {
    tryCatch({
      omega_B_arr <- posterior::as_draws_array(
        stan_fit$draws(variables = "Omega_B")
      )
      sigma_B_arr <- posterior::as_draws_array(
        stan_fit$draws(variables = "Sigma_B")
      )
      omega_P_arr <- posterior::as_draws_array(
        stan_fit$draws(variables = "Omega_P")
      )
      sigma_P_arr <- posterior::as_draws_array(
        stan_fit$draws(variables = "Sigma_P")
      )
      cov_summaries$Omega_B <- summarize_matrix_draws(omega_B_arr,
                                                      "Omega_B", K, K)
      cov_summaries$Sigma_B <- summarize_matrix_draws(sigma_B_arr,
                                                      "Sigma_B", K, K)
      cov_summaries$Omega_P <- summarize_matrix_draws(omega_P_arr,
                                                      "Omega_P", 5L, 5L)
      cov_summaries$Sigma_P <- summarize_matrix_draws(sigma_P_arr,
                                                      "Sigma_P", 5L, 5L)
      dimnames(cov_summaries$Omega_B) <- list(antigens, antigens)
      dimnames(cov_summaries$Sigma_B) <- list(antigens, antigens)
      dimnames(cov_summaries$Omega_P) <- list(param_names, param_names)
      dimnames(cov_summaries$Sigma_P) <- list(param_names, param_names)
    }, error = function(e) {
      cli::cli_warn("Kronecker matrices not extracted: {e$message}")
    })
  }

  # ---- Parameter correlation (model_1 only — per-biomarker Omega_P[k]) ----
  # model_1 generates array[K] corr_matrix[P] Omega_P; cmdstanr names
  # cells Omega_P[k,p,q]. model_2's single Omega_P is handled above.
  if (!has_kron) {
    P <- length(param_names)
    tryCatch({
      omega_P_arr <- posterior::as_draws_array(
        stan_fit$draws(variables = "Omega_P")
      )
      omega_P_list <- summarize_matrix_array(
        omega_P_arr, "Omega_P", K, P, P
      )
      for (k in seq_len(K)) {
        dimnames(omega_P_list[[k]]) <- list(param_names, param_names)
      }
      names(omega_P_list) <- antigens
      cov_summaries$Omega_P <- omega_P_list
    }, error = function(e) {
      cli::cli_warn("model_1 Omega_P not extracted: {e$message}")
    })
  }

  # ---- log_lik for LOO ----
  tryCatch({
    cov_summaries$log_lik <- posterior::as_draws_matrix(
      stan_fit$draws(variables = "log_lik")
    )
  }, error = function(e) {
    cli::cli_warn("log_lik not extracted: {e$message}")
  })

  return(list(
    sr_tibble     = sr_tibble,
    cov_summaries = cov_summaries
  ))
}
