#' Format median and credible interval as a single string
#'
#' @param med Median value.
#' @param lo Lower bound.
#' @param hi Upper bound.
#' @param digits Number of digits.
#' @param sci Logical; if TRUE uses scientific notation.
#'
#' @return A string like "1.23 (0.50–2.00)".
#'
#' @export
fmt_mci <- function(med, lo, hi, digits = 2, sci = FALSE) {
  f <- function(x) {
    if (sci) formatC(x, format = "e", digits = digits) else formatC(x, format = "f", digits = digits)
  }
  sprintf("%s (%s–%s)", f(med), f(lo), f(hi))
}
