#' Summarise natural-scale draws to posterior median + 95% credible interval
#'
#' @param df Output of [add_natural_scale()].
#' @param group_vars Character vector of grouping columns
#'   (e.g. `c("antigen", "Iso_type")`, or
#'   `c("age_group", "antigen", "Iso_type")`).
#' @return One row per group with `{y0,y1,t1,alpha,rho}_{med,lo,hi}` columns.
#'   (`alpha_*` follow the unit chosen in [add_natural_scale()].)
#' @export
summarise_pop_params <- function(df, group_vars = c("antigen", "Iso_type")) {
  q_lo <- function(x) stats::quantile(x, 0.025, na.rm = TRUE)
  q_hi <- function(x) stats::quantile(x, 0.975, na.rm = TRUE)

  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      y0_med = stats::median(.data$y0_natural, na.rm = TRUE),
      y0_lo = q_lo(.data$y0_natural), y0_hi = q_hi(.data$y0_natural),
      y1_med = stats::median(.data$y1_natural, na.rm = TRUE),
      y1_lo = q_lo(.data$y1_natural), y1_hi = q_hi(.data$y1_natural),
      t1_med = stats::median(.data$t1_natural, na.rm = TRUE),
      t1_lo = q_lo(.data$t1_natural), t1_hi = q_hi(.data$t1_natural),
      alpha_med = stats::median(.data$alpha_natural, na.rm = TRUE),
      alpha_lo = q_lo(.data$alpha_natural), alpha_hi = q_hi(.data$alpha_natural), # nolint: line_length_linter.
      rho_med = stats::median(.data$rho_natural, na.rm = TRUE),
      rho_lo = q_lo(.data$rho_natural), rho_hi = q_hi(.data$rho_natural),
      .groups = "drop"
    )
}
