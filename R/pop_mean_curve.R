#' Population-mean trajectory (median + 95% CrI) over a time grid
#'
#' Evaluates the two-phase curve at the population-mean draws for each isotype.
#' The decay rate is kept in **day^-1** here because the time grid and
#' serodynamics:::ab() operate in days (contrast with the per-year display in
#' the parameter tables).
#'
#' @param model_obj Fitted `sr_model`.
#' @param t_grid Time grid in days.
#' @return Tibble `Iso_type, t, med, lo, hi, iso` (`iso` an ordered factor).
#' @export
pop_mean_curve <- function(model_obj, t_grid = seq(0, 210, by = 5)) {
  pop <- attr(model_obj, "population_params")
  purrr::map_dfr(c("IgG", "IgA"), function(iso) {
    pop |>
      dplyr::filter(.data$Population_Parameter == "mu.par",
                    .data$Iso_type == iso) |>
      dplyr::select("Chain", "Iteration", "Iso_type", "Parameter", "value") |>
      tidyr::pivot_wider(names_from = "Parameter", values_from = "value") |>
      dplyr::mutate(
        y0_nat    = exp(.data$y0),
        y1_nat    = exp(.data$y0) + exp(.data$y1),
        t1_nat    = exp(.data$t1),
        alpha_nat = exp(.data$alpha),          # day^-1
        shape_nat = exp(.data$shape) + 1
      ) |>
      # each retained draw x each grid point; ~2000 draws x length(t_grid) rows per isotype # nolint: line_length_linter.
      tidyr::crossing(t = t_grid) |>
      # TODO: replace with exported serodynamics API once ab() is exported (tracked separately).
      dplyr::mutate(res = serodynamics:::ab(.data$t, .data$y0_nat, .data$y1_nat,
                                            .data$t1_nat, .data$alpha_nat,
                                            .data$shape_nat)) |>
      dplyr::filter(is.finite(.data$res), .data$res > 0) |>
      dplyr::summarise(
        med = stats::quantile(.data$res, 0.50, na.rm = TRUE),
        lo  = stats::quantile(.data$res, 0.025, na.rm = TRUE),
        hi  = stats::quantile(.data$res, 0.975, na.rm = TRUE),
        .by = c("Iso_type", "t"))
  }) |>
    dplyr::mutate(iso = factor(.data$Iso_type, levels = c("IgG", "IgA")))
}
