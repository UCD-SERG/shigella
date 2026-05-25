# ─── Internal extraction helpers ─────────────────────────────────────────────
# Consolidated here (from separate extract_*.R files) to eliminate cross-file
# object_usage_linter false positives on dotted-prefix internal helper calls.

# Helper: extract all chain draws for one parameter cell `pname[subj, iso]`.
#' @keywords internal
#' @noRd
.extract_draws_for_cell <- function(col_name, draws_df, pname,
                                    ids, antigens, stratification, n_chain) {
  m        <- regmatches(col_name, regexec("\\[(\\d+),(\\d+)\\]",
                                           col_name))[[1]]
  subj_idx <- as.integer(m[2])
  iso_idx  <- as.integer(m[3])
  sub_df   <- draws_df[, c(".chain", ".iteration", col_name)]

  lapply(seq_len(n_chain), function(ch) {
    chain_data <- sub_df[sub_df$.chain == ch, ]
    tibble::tibble(
      Iteration      = chain_data$.iteration,
      Chain          = ch,
      Parameter      = pname,
      Iso_type       = antigens[iso_idx],
      Stratification = stratification,
      Subject        = ids[subj_idx],
      value          = chain_data[[col_name]]
    )
  })
}

# Helper: extract draws for all param_names into a single tibble.
#' @keywords internal
#' @noRd
.extract_param_draws <- function(param_names,
                                 draws_df,
                                 N,
                                 K,
                                 ids,
                                 antigens,
                                 stratification,
                                 n_chain) {
  draws_per_param <- lapply(param_names, function(pname) {
    matching_cols <- grep(paste0("^", pname, "\\["), colnames(draws_df),
                          value = TRUE)
    if (length(matching_cols) != N * K) {
      cli::cli_abort(
        "Expected {N * K} {.var {pname}} draws; got {length(matching_cols)}."
      )
    }
    draws_per_col <- lapply(
      matching_cols, .extract_draws_for_cell,
      draws_df       = draws_df,
      pname          = pname,
      ids            = ids,
      antigens       = antigens,
      stratification = stratification,
      n_chain        = n_chain
    )
    dplyr::bind_rows(unlist(draws_per_col, recursive = FALSE))
  })

  dplyr::bind_rows(draws_per_param)
}

# Helper: extract Omega_eps and Sigma_eps from a cmdstanr fit. Model 2 only.
#' @keywords internal
#' @noRd
.extract_residual_cov_stan <- function(stan_fit, K, antigens) {
  tryCatch({
    omega_eps_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Omega_eps")
    )
    sigma_eps_arr <- posterior::as_draws_array(
      stan_fit$draws(variables = "Sigma_eps")
    )
    omega_eps <- .summarize_matrix_draws(omega_eps_arr, "Omega_eps", K, K)
    sigma_eps <- .summarize_matrix_draws(sigma_eps_arr, "Sigma_eps", K, K)
    dimnames(omega_eps) <- list(antigens, antigens)
    dimnames(sigma_eps) <- list(antigens, antigens)
    list(Omega_eps = omega_eps, Sigma_eps = sigma_eps)
  }, error = function(e) {
    cli::cli_warn("Omega_eps/Sigma_eps not extracted: {e$message}")
    list()
  })
}

# Helper: extract Omega_B, Sigma_B, Omega_P, Sigma_P from a cmdstanr fit.
# Model 2 only.
#' @keywords internal
#' @noRd
.extract_kron_matrices <- function(stan_fit, K, param_names, antigens) {
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

    omega_B <- .summarize_matrix_draws(omega_B_arr, "Omega_B", K, K)
    sigma_B <- .summarize_matrix_draws(sigma_B_arr, "Sigma_B", K, K)
    P <- length(param_names)
    omega_P <- .summarize_matrix_draws(omega_P_arr, "Omega_P", P, P)
    sigma_P <- .summarize_matrix_draws(sigma_P_arr, "Sigma_P", P, P)

    dimnames(omega_B) <- list(antigens, antigens)
    dimnames(sigma_B) <- list(antigens, antigens)
    dimnames(omega_P) <- list(param_names, param_names)
    dimnames(sigma_P) <- list(param_names, param_names)

    list(Omega_B = omega_B, Sigma_B = sigma_B,
         Omega_P = omega_P, Sigma_P = sigma_P)
  }, error = function(e) {
    cli::cli_warn("Kronecker matrices not extracted: {e$message}")
    list()
  })
}

# Helper: extract per-biomarker Omega_P list from a model_1 cmdstanr fit.
# model_1 generates array[K] corr_matrix[P] Omega_P; cmdstanr names
# cells Omega_P[k,p,q]. Returns a named list of K matrices.
#' @keywords internal
#' @noRd
.extract_model1_omega_p_stan <- function(stan_fit, K, param_names, antigens) {
  P <- length(param_names)
  tryCatch({
    omega_P_arr  <- posterior::as_draws_array(
      stan_fit$draws(variables = "Omega_P")
    )
    omega_P_list <- .summarize_matrix_array(omega_P_arr, "Omega_P", K, P, P)
    for (k in seq_len(K)) {
      dimnames(omega_P_list[[k]]) <- list(param_names, param_names)
    }
    names(omega_P_list) <- antigens
    list(Omega_P = omega_P_list)
  }, error = function(e) {
    cli::cli_warn("model_1 Omega_P not extracted: {e$message}")
    list()
  })
}

# Helper: extract the log_lik draws matrix from a cmdstanr fit.
#' @keywords internal
#' @noRd
.extract_log_lik_stan <- function(stan_fit) {
  tryCatch({
    list(log_lik = posterior::as_draws_matrix(
      stan_fit$draws(variables = "log_lik")
    ))
  }, error = function(e) {
    cli::cli_warn("log_lik not extracted: {e$message}")
    list()
  })
}

# ─── Main function ────────────────────────────────────────────────────────────

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
#' @param param_names Character vector of parameter names to extract from the
#'   Stan model's \code{generated quantities} block. Must match the variable
#'   names exactly as declared in both \file{inst/stan/model_1.stan} and
#'   \file{inst/stan/model_2.stan}. Defaults to
#'   \code{c("y0", "y1", "t1", "alpha", "shape")}.
#' @return list with sr_tibble and cov_summaries
#' @example inst/examples/postprocess_stan_output-examples.R
#' @export
postprocess_stan_output <- function(stan_fit,
                                    ids,
                                    antigens,
                                    model = c("model_2", "model_1"),
                                    stratification = "None",
                                    param_names = c("y0", "y1", "t1",
                                                    "alpha", "shape")) {

  model <- match.arg(model)
  has_kron <- identical(model, "model_2")

  if (!requireNamespace("posterior", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg posterior} required for cmdstanr
                   postprocessing.")
  }

  avail <- stan_fit$metadata()$stan_variables
  missing_vars <- setdiff(param_names, avail)
  if (length(missing_vars) > 0) {
    cli::cli_abort(c(
      "{.arg param_names} contains variables absent from the Stan model.",
      "x" = "Missing: {.val {missing_vars}}",
      "i" = "Available Stan variables: {.val {avail}}"
    ))
  }

  N <- length(ids)
  K <- length(antigens)

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

  if (has_kron) {
    cov_summaries <- c(
      .extract_residual_cov_stan(stan_fit, K, antigens),
      .extract_kron_matrices(stan_fit, K, param_names, antigens)
    )
  } else {
    cov_summaries <- .extract_model1_omega_p_stan(
      stan_fit, K, param_names, antigens
    )
  }

  cov_summaries <- c(cov_summaries, .extract_log_lik_stan(stan_fit))

  return(list(
    sr_tibble     = sr_tibble,
    cov_summaries = cov_summaries
  ))
}
