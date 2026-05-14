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
#' @examples
#' \dontrun{
#' source(system.file(
#'   "examples",
#'   "postprocess_stan_output-examples.R",
#'   package = "shigella"
#' ))
#' }
#' @export
postprocess_stan_output <- function(stan_fit,
                                     ids,
                                     antigens,
                                     model = c("model_2", "model_1"),
                                     stratification = "None") {

  model <- match.arg(model)
  has_kron <- model %in% c("model_2", "model_1")

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' required for cmdstanr postprocessing.")
  }

  param_names <- c("y0", "y1", "t1", "alpha", "shape")
  N <- length(ids)
  K <- length(antigens)

  # cmdstanr returns draws via $draws() which is a draws_array
  # Convert to data frame format for processing
  draws_df <- posterior::as_draws_df(
    stan_fit$draws(variables = param_names)
  )
  n_iter  <- max(draws_df$.iteration)
  n_chain <- max(draws_df$.chain)

  out_list <- list()
  row_counter <- 1L

  for (p in seq_along(param_names)) {
    pname <- param_names[p]
    # Find columns matching pname[i,k]
    matching_cols <- grep(paste0("^", pname, "\\["), colnames(draws_df), value = TRUE)
    if (length(matching_cols) != N * K) {
      stop(sprintf("Expected %d %s draws; got %d", N * K, pname,
                   length(matching_cols)))
    }

    for (col_name in matching_cols) {
      m <- regmatches(col_name, regexec("\\[(\\d+),(\\d+)\\]", col_name))[[1]]
      subj_idx <- as.integer(m[2])
      iso_idx  <- as.integer(m[3])

      # Extract draws for this parameter index
      sub_df <- draws_df[, c(".chain", ".iteration", col_name)]

      for (ch in seq_len(n_chain)) {
        chain_data <- sub_df[sub_df$.chain == ch, ]
        out_list[[row_counter]] <- tibble::tibble(
          Iteration      = chain_data$.iteration,
          Chain          = ch,
          Parameter      = pname,
          Iso_type       = antigens[iso_idx],
          Stratification = stratification,
          Subject        = ids[subj_idx],
          value          = chain_data[[col_name]]
        )
        row_counter <- row_counter + 1L
      }
    }
  }

  sr_tibble <- dplyr::bind_rows(out_list)

  cov_summaries <- list()

  # ---- Residual covariance (all models) ----
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
      cov_summaries$Omega_B <- summarize_matrix_draws(omega_B_arr, "Omega_B", K, K)
      cov_summaries$Sigma_B <- summarize_matrix_draws(sigma_B_arr, "Sigma_B", K, K)
      cov_summaries$Omega_P <- summarize_matrix_draws(omega_P_arr, "Omega_P", 5L, 5L)
      cov_summaries$Sigma_P <- summarize_matrix_draws(sigma_P_arr, "Sigma_P", 5L, 5L)
      dimnames(cov_summaries$Omega_B) <- list(antigens, antigens)
      dimnames(cov_summaries$Sigma_B) <- list(antigens, antigens)
      dimnames(cov_summaries$Omega_P) <- list(param_names, param_names)
      dimnames(cov_summaries$Sigma_P) <- list(param_names, param_names)
    }, error = function(e) {
      cli::cli_warn("Kronecker matrices not extracted: {e$message}")
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

# Helper: summarize a draws_array of a matrix variable to a single matrix
# by taking the median across all draws.
summarize_matrix_draws <- function(draws_arr, var_name, nrow, ncol) {
  result <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  var_dim <- dimnames(draws_arr)$variable
  for (i in seq_len(nrow)) {
    for (j in seq_len(ncol)) {
      cell_name <- sprintf("%s[%d,%d]", var_name, i, j)
      if (cell_name %in% var_dim) {
        cell_draws <- as.numeric(draws_arr[, , cell_name])
        result[i, j] <- median(cell_draws, na.rm = TRUE)
      }
    }
  }
  result
}
