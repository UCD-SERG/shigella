#' Convenience: fitted models -> tidy population-parameter summary table
#'
#' Chains [extract_mu_draws()], [add_natural_scale()] and
#' [summarise_pop_params()] over a named list of fitted models. This is the
#' single entry point used by the parameter tables (Table 2, S4 Table, forest
#' plot).
#'
#' @param model_list Named list `antigen_label = sr_model`.
#' @param antigen_levels Factor order for the `antigen` column.
#' @param alpha_unit Passed to [add_natural_scale()].
#' @return A `summarise_pop_params()` table over all supplied antigens.
#' @export
pop_param_summary <- function(model_list,
                              antigen_levels = c(
                                "IpaB", "Sf2a", "Sf3a",
                                "Sf6", "Sonnei"
                              ),
                              alpha_unit = "per_year") {
  purrr::imap(model_list, ~ extract_mu_draws(.x, antigen_label = .y)) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      antigen  = factor(.data$antigen, levels = antigen_levels),
      Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA"))
    ) |>
    add_natural_scale(alpha_unit = alpha_unit) |>
    summarise_pop_params(group_vars = c("antigen", "Iso_type"))
}
