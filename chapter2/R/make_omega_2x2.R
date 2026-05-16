#' Build a 2x2 correlation matrix with off-diagonal rho
#' @param rho numeric off-diagonal correlation
#' @return 2x2 matrix
make_omega_2x2 <- function(rho) {
  matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
}
