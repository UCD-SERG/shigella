#' Reshape Shigella longitudinal data for serodynamics workflows
#'
#' Filters a dataset to a given study and antigen column, standardizes column
#' names, and returns a visit-ordered long dataset suitable for conversion to
#' `serodynamics::as_case_data()`.
#'
#' @param data A data frame containing longitudinal measurements.
#' @param study_filter Character scalar. Value of `study_name` to keep
#'   (e.g. "SOSAR").
#' @param antigen Unquoted column name for the antigen measurement
#'   (e.g. n_ipab_MFI).
#'
#' @return A tibble with standardized columns:
#' \describe{
#'   \item{index_id}{Participant ID (copied from `sid`).}
#'   \item{antigen_iso}{Isotype label (copied from `isotype_name`).}
#'   \item{visit}{Visit label (copied from `timepoint`).}
#'   \item{timeindays}{Time since infection (copied from `Actual day`).}
#'   \item{result}{Antibody measurement (from `antigen`).}
#' }
#'
#' @examples
#' \dontrun{
#' dat_long <- process_shigella_data(df, "SOSAR", n_ipab_MFI)
#' dL <- serodynamics::as_case_data(dat_long,
#'   id_var = "index_id", biomarker_var = "antigen_iso",
#'   time_in_days = "timeindays", value_var = "result"
#' )
#' }
#'
#' @importFrom rlang ensym
#' @export
process_shigella_data <- function(data, study_filter, antigen) {
  antigen_col <- rlang::ensym(antigen)

  data |>
    dplyr::filter(.data$study_name == study_filter) |>
    dplyr::select(
      .data$isotype_name,
      .data$sid,
      .data$timepoint,
      .data$`Actual day`,
      !!antigen_col
    ) |>
    dplyr::mutate(
      index_id   = .data$sid,
      antigen_iso = .data$isotype_name,
      visit      = .data$timepoint,
      timeindays = .data$`Actual day`,
      result     = !!antigen_col
    ) |>
    dplyr::group_by(.data$index_id, .data$antigen_iso) |>
    dplyr::arrange(.data$visit, .by_group = TRUE) |>
    dplyr::mutate(visit_num = rank(.data$visit, ties.method = "first")) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$timeindays))
}
