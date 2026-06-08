#' Main-text individual model-comparison figure (Sf2a / Sonnei / Sf3a)
#'
#' @param models Named list with `overall_Sf2a`, `serotype_sf2a`,
#'   `overall_Sonnei`, `serotype_sonnei`, `overall_Sf3a`, `serotype_sf3a`,
#'   `combined_sf3a` fitted models.
#' @param datasets Named list with `dL_clean_sf2a`, `dL_serotype_sf2a`,
#'   `dL_clean_sonnei`, `dL_serotype_sonnei`, `dL_clean_sf3a`, `dL_serotype_sf3a`,
#'   `dL_combined_sf3a` case data.
#' @param sids Length-3 character vector of representative subject ids
#'   (Sf2a, Sonnei, Sf3a). Defaults match the manuscript.
#' @return A patchwork figure.
#' @export
figure_individual_comparison <- function(models, datasets,
                                         sids = c("SOSAR-22008", "SOSAR-12005",
                                                  "SOSAR-22007")) {
  row1 <- build_figure4_row(sids[1], models$overall_Sf2a, models$serotype_sf2a,
    datasets$dL_clean_sf2a, datasets$dL_serotype_sf2a, NULL,
    "A) S. flexneri 2a (n = 17): 2-way comparison", "SF2a")
  row2 <- build_figure4_row(sids[2], models$overall_Sonnei, models$serotype_sonnei,
    datasets$dL_clean_sonnei, datasets$dL_serotype_sonnei, NULL,
    "B) S. sonnei (n = 11): 2-way comparison", "Sonnei")
  row3 <- build_figure4_row(sids[3], models$overall_Sf3a, models$serotype_sf3a,
    datasets$dL_clean_sf3a, datasets$dL_serotype_sf3a, models$combined_sf3a,
    "C) S. flexneri 3a (n = 8; combined n = 25): 3-way comparison", "SF3a")
  assemble_comparison_rows(list(row1, row2, row3))
}
