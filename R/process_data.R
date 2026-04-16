#' Process Raw Shigella Luminex Data into Case-Data Format
#'
#' Filters a compiled Shigella dataset to a single study, selects a given
#' antigen MFI column, restructures it into long per-visit format, and
#' converts it to a case-data object.
#'
#' @param data A data frame containing the compiled Shigella Luminex data
#'   (e.g., loaded from the compiled spreadsheet, sheet `"Compiled"`).
#'   Must contain columns `study_name`, `isotype_name`, `sid`, `timepoint`,
#'   `Actual day`, `cohort_name`, and the antigen column specified by
#'   `antigen_col`.
#' @param study_filter A character string identifying the study to retain
#'   (e.g., `"SOSAR"`).
#' @param antigen_col An unquoted column name (tidy-select style) for the
#'   MFI variable to extract (e.g., `n_ipab_MFI`).
#'
#' @return A `case_data` object with columns `index_id`, `antigen_iso`,
#'   `visit`, `timeindays`, `result`, and `cohort_name`.
#'
#' @examples
#' \dontrun{
#' library(readxl)
#' df <- read_excel("path/to/compiled_data.xlsx", sheet = "Compiled")
#' dL_ipab <- process_shigella_data(df, "SOSAR", n_ipab_MFI)
#' }
#'
#' @importFrom dplyr filter select mutate group_by arrange ungroup
#' @importFrom rlang ensym
#' @export
process_shigella_data <- function(data, study_filter, antigen_col) {
  antigen_col_sym <- rlang::ensym(antigen_col)
  out <- data |>
    dplyr::filter(.data$study_name == study_filter) |>
    dplyr::select(
      "isotype_name", "sid", "timepoint",
      "Actual day", "cohort_name",
      !!antigen_col_sym
    ) |>
    dplyr::mutate(
      index_id    = .data$sid,
      antigen_iso = .data$isotype_name,
      visit       = .data$timepoint,
      timeindays  = .data[["Actual day"]],
      result      = !!antigen_col_sym
    ) |>
    dplyr::group_by(.data$index_id, .data$antigen_iso) |>
    dplyr::arrange(.data$visit) |>
    dplyr::mutate(
      visit_num = rank(.data$visit, ties.method = "first")
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$timeindays))

  .as_case_data(
    out,
    id_var        = "index_id",
    biomarker_var = "antigen_iso",
    time_in_days  = "timeindays",
    value_var     = "result"
  )
}


#' Load All Shigella Antigen Case-Data Objects from a Compiled Spreadsheet
#'
#' Convenience wrapper that calls [process_shigella_data()] for each of the
#' five standard antigens and returns a named list of case-data objects.
#'
#' @param data A compiled Shigella data frame (same format as required by
#'   [process_shigella_data()]).
#' @param study_filter A character string for study filtering. Default
#'   `"SOSAR"`.
#'
#' @return A named list with elements `IpaB`, `Sf2a`, `Sf3a`, `Sf6`,
#'   `Sonnei`, each a case-data object.
#'
#' @examples
#' \dontrun{
#' df  <- readxl::read_excel("path/to/compiled_data.xlsx", sheet = "Compiled")
#' dL  <- load_all_antigens(df)
#' dL$IpaB
#' }
#'
#' @export
load_all_antigens <- function(data, study_filter = "SOSAR") {
  list(
    IpaB   = process_shigella_data(data, study_filter, n_ipab_MFI),
    Sf2a   = process_shigella_data(data, study_filter, n_sf2aospbsa_MFI),
    Sf3a   = process_shigella_data(data, study_filter, n_sf3aospbsa_MFI),
    Sf6    = process_shigella_data(data, study_filter, n_sf6ospbsa_MFI),
    Sonnei = process_shigella_data(data, study_filter, n_sonneiospbsa_MFI)
  )
}
