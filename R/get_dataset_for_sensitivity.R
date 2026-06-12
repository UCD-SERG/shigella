#' Map an (antigen, model) pair to its case-data object
#'
#' @param antigen `"Sonnei"`, `"Sf3a"`, or `"Sf2a"`.
#' @param model `"overall"` or `"serotype"`.
#' @param datasets Named list holding the case-data objects:
#'   `dL_clean_sonnei_new`, `dL_clean_sf3a_new`, `dL_clean_sf2a_new`
#'   (overall) and `dL_serotype_sonnei`, `dL_serotype_sf3a`,
#'   `dL_serotype_sf2a` (serotype).
#' @return The matching case-data object, or `NULL`.
#' @export
get_dataset_for_sensitivity <- function(antigen, model, datasets) {
  key <- if (model == "overall") {
    c(
      Sonnei = "dL_clean_sonnei_new", Sf3a = "dL_clean_sf3a_new",
      Sf2a = "dL_clean_sf2a_new"
    )[antigen]
  } else if (model == "serotype") {
    c(
      Sonnei = "dL_serotype_sonnei", Sf3a = "dL_serotype_sf3a",
      Sf2a = "dL_serotype_sf2a"
    )[antigen]
  } else {
    NA_character_
  }
  if (is.na(key)) {
    return(NULL)
  }
  datasets[[key]]
}
