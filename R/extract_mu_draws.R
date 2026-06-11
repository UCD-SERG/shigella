#' Extract population-mean (`mu.par`) posterior draws in wide form
#'
#' Pulls the `mu.par` rows from a fitted model's `population_params` attribute
#' and pivots them so each row is one posterior draw (Iteration x Chain x
#' isotype) with one column per kinetic parameter on the **log/JAGS scale**
#' (`y0`, `y1`, `t1`, `alpha`, `shape`).
#'
#' @param fit_object An `sr_model` object from [serodynamics::run_mod_pop()].
#' @param antigen_label Optional antigen label added as an `antigen` column,
#'   used when stacking several antigens before summarising.
#' @return A tibble with columns `Iteration`, `Chain`, `Iso_type`,
#'   `y0`, `y1`, `t1`, `alpha`, `shape` (and `antigen` if `antigen_label` given).
#' @export
extract_mu_draws <- function(fit_object, antigen_label = NULL) {
  w <- attr(fit_object, "population_params") |>
    dplyr::filter(.data$Population_Parameter == "mu.par") |>
    dplyr::select("Iteration", "Chain", "Iso_type", "Parameter", "value") |>
    tidyr::pivot_wider(names_from = "Parameter", values_from = "value")

  if (!is.null(antigen_label)) {
    w <- dplyr::mutate(w, antigen = antigen_label)
  }
  w
}
