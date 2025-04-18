#' @title Process Shigella Antibody Data for Longitudinal Modeling
#' @author Kwan Ho Lee
#' @description
#' `process_shigella_data()` filters, restructures, and formats longitudinal
#' Shigella antibody data for downstream modeling and visualization.
#' The function subsets by study, reshapes the data with harmonized column names,
#' and converts it to a standardized case data format.
#'
#' This function is especially useful for preparing input for models such as
#' those in the `shigella` and `serodynamics` packages.
#'
#' @param data A data.frame or tibble with Shigella antibody data.
#' @param study_filter A character string indicating the study name to filter.
#' @param antigen The unquoted name of the column representing antigen-specific
#' fluorescence intensity (e.g., `n_ipab_MFI`).
#'
#' @return A case data object compatible with `serodynamics::as_case_data()`,
#' containing harmonized fields: `index_id`, `antigen_iso`, `timeindays`, `result`.
#'
#' @examples
#' \dontrun{
#' df <- readxl::read_excel("example.xlsx", sheet = "Compiled")
#' processed <- process_shigella_data(df, study_filter = "SOSAR", antigen = n_ipab_MFI)
#' }
#'
#' @export
process_shigella_data <- function(data, study_filter, antigen) {
  # 1. Filter the data for the specific study
  filtered_data <- data |>
    dplyr::filter(study_name == study_filter)
  
  # 2. Capture the column name of the antigen
  antigen_col <- rlang::ensym(antigen)
  
  # 3. Manipulate and restructure the data
  processed_data <- filtered_data |>
    dplyr::select(isotype_name, sid, timepoint, `Actual day`, !!antigen_col) |>
    dplyr::mutate(
      index_id    = sid,
      antigen_iso = isotype_name,
      visit       = timepoint,
      timeindays  = `Actual day`,
      result      = !!antigen_col
    ) |>
    dplyr::group_by(index_id, antigen_iso) |>
    dplyr::arrange(visit) |>
    dplyr::mutate(visit_num = rank(visit, ties.method = "first")) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(timeindays))
  
  # 4. Convert directly into "case_data" object and return
  processed_data |>
    serodynamics::as_case_data(
      id_var         = "index_id",
      biomarker_var  = "antigen_iso",
      time_in_days   = "timeindays",
      value_var      = "result"
    )
}
