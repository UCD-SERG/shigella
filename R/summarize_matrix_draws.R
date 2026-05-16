#' @keywords internal
#' @noRd
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

# For `array[n_arr] matrix[nrow, ncol]` Stan variables (e.g. model_1's
# `array[K] corr_matrix[P] Omega_P`), cmdstanr names cells `var[k,i,j]`.
# Returns a list of n_arr matrices, each nrow x ncol.
#' @keywords internal
#' @noRd
summarize_array_of_matrix_draws <- function(draws_arr, var_name,
                                             n_arr, nrow, ncol) {
  var_dim <- dimnames(draws_arr)$variable
  lapply(seq_len(n_arr), function(k) {
    mat <- matrix(NA_real_, nrow = nrow, ncol = ncol)
    for (i in seq_len(nrow)) {
      for (j in seq_len(ncol)) {
        cell_name <- sprintf("%s[%d,%d,%d]", var_name, k, i, j)
        if (cell_name %in% var_dim) {
          mat[i, j] <- median(as.numeric(draws_arr[, , cell_name]),
                              na.rm = TRUE)
        }
      }
    }
    mat
  })
}
