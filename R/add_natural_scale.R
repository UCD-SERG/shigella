#' Back-transform log-scale `mu.par` draws to the natural parameter scale
#'
#' Applies the manuscript back-transformation to each draw:
#' \itemize{
#'   \item `y0_natural    = exp(y0)`
#'   \item `y1_natural    = exp(y1) + exp(y0)` (the `y1` column is
#'     `log(y1 - y0)`)
#'   \item `t1_natural    = exp(t1)`
#'   \item `alpha_natural = exp(alpha) * k`
#'     (k = 365 for per-year, 1 for per-day)
#'   \item `rho_natural   = exp(shape) + 1`
#' }
#'
#' @param df Output of [extract_mu_draws()].
#' @param alpha_unit `"per_year"` (default, for display tables) multiplies the
#'   decay rate by 365; `"per_day"` leaves it in day^-1.
#' @return `df` with the five `*_natural` columns added.
#' @export
add_natural_scale <- function(df, alpha_unit = c("per_year", "per_day")) {
  alpha_unit <- match.arg(alpha_unit)
  k <- if (alpha_unit == "per_year") 365 else 1

  dplyr::mutate(
    df,
    y0_natural    = exp(.data$y0),
    y1_natural    = exp(.data$y1) + exp(.data$y0),
    t1_natural    = exp(.data$t1),
    alpha_natural = exp(.data$alpha) * k,
    rho_natural   = exp(.data$shape) + 1
  )
}
