#' Format a median (lo--hi) triple for display
#'
#' @param med,lo,hi Numeric vectors.
#' @param digits Decimal places (fixed notation) or significant digits (`sci`).
#' @param sci If `TRUE`, use scientific notation.
#' @return Character vector `"med (lo-hi)"`.
#' @export
fmt_mci <- function(med, lo, hi, digits = 2, sci = FALSE) {
  f <- function(x) {
    if (sci) formatC(x, format = "e", digits = digits)
    else     formatC(x, format = "f", digits = digits)
  }
  sprintf("%s (%s\u2013%s)", f(med), f(lo), f(hi))
}
