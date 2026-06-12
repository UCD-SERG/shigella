#' Assemble the full population-trajectory figure (Fig 5, panels A-E)
#'
#' @param models Named list: `IpaB`, `IpaB_under5`, `IpaB_plus5`,
#'   `serotype_sonnei`, `serotype_sf2a`, `combined_sf3a`, `overall_Sf6`.
#' @param datasets Named list: `dL_clean_Ipab`, `dL_serotype_sonnei`,
#'   `dL_serotype_sf2a`, `dL_combined_sf3a`, `dL_clean_sf6`.
#' @return A patchwork figure (5 stacked panels).
#' @export
figure_population_trajectories <- function(models, datasets) {
  p_ipab <- figure_ipab_with_age(
    models$IpaB, models$IpaB_under5, models$IpaB_plus5, datasets$dL_clean_Ipab
  )
  p_sonnei <- figure_with_individuals(models$serotype_sonnei, datasets$dL_serotype_sonnei, # nolint: line_length_linter.
    "B) S. sonnei (Serotype-specific, n=11)",
    individual_alpha = 0.22
  )
  p_sf2a <- figure_with_individuals(models$serotype_sf2a, datasets$dL_serotype_sf2a, # nolint: line_length_linter.
    "C) S. flexneri 2a (Serotype-specific, n=17)",
    individual_alpha = 0.22
  )
  p_sf3a <- figure_with_individuals(models$combined_sf3a, datasets$dL_combined_sf3a, # nolint: line_length_linter.
    "D) S. flexneri 3a (Combined, n=25)",
    individual_alpha = 0.22
  )
  p_sf6 <- figure_with_individuals(models$overall_Sf6, datasets$dL_clean_sf6,
    "E) S. flexneri 6 (Overall, n=48)",
    individual_alpha = 0.22
  )

  (p_ipab / p_sonnei / p_sf2a / p_sf3a / p_sf6) +
    patchwork::plot_layout(ncol = 1, heights = c(1.2, 1, 1, 1, 1)) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5))
    )
}
