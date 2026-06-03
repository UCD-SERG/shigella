#' Validate that a matrix is a correlation matrix
#'
#' @param M Numeric matrix to validate.
#' @param name Character; name to use in error messages.
#' @param tol Numeric tolerance for symmetry and unit-diagonal checks.
#' @keywords internal
#' @noRd
.validate_corr_matrix <- function(M, name, tol = 1e-8) {
  if (!is.matrix(M) || !is.numeric(M)) {
    cli::cli_abort("{.arg {name}} must be a numeric matrix.")
  }
  if (nrow(M) != ncol(M)) {
    cli::cli_abort("{.arg {name}} must be square; got {nrow(M)} x {ncol(M)}.")
  }
  if (max(abs(M - t(M))) > tol) {
    cli::cli_abort("{.arg {name}} must be symmetric.")
  }
  if (max(abs(diag(M) - 1)) > tol) {
    cli::cli_abort("{.arg {name}} must have unit diagonal.")
  }
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < -tol) {
    cli::cli_abort(c(
      "{.arg {name}} must be positive semi-definite.",
      "i" = "Smallest eigenvalue: {min(eig)}"
    ))
  }
  invisible(M)
}
