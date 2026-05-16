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
