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
