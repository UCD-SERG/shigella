#' Build single-antigen case data from the compiled Shigella table
#'
#' Filters the compiled MFI table to one study and one antigen and reshapes it
#' into the long "case data" layout expected by `serodynamics`/`serocalculator`.
#' Mirrors the original manuscript `process_shigella_data()`, and additionally
#' **keeps `cohort_name` (infecting serotype) and `age`**, which live in the
#' compiled sheet (columns J and M; see the raw-data screenshot). Retaining them
#' lets the `subset_*()` helpers below build the serotype-specific and
#' age-stratified datasets directly -- no separate metadata join is needed.
#'
#' @details
#' **ID column.** The original helper created `index_id = sid` and attached
#' `ids_varname = "index_id"`. The shipped case-data objects (`dL_clean_sf2a`
#' etc.) and the current working data instead carry an `id` column, identical in
#' value to `sid`. To reproduce those objects exactly, this function
#' standardises the id column to `id` via `serodynamics::as_case_data()`
#' (`id_var = "id"`), rather than relying on any internal renaming. If a
#' downstream step needs `index_id`, change `id_var` here -- but the shipped
#' `.rda` inputs and the manuscript figures use `id`.
#'
#' **cohort_name / age.** The original bare object `dL_clean_sf2a` did not carry
#' these two columns, so the base (non-`_new`) datasets built here are a harmless
#' superset of the originals (two extra metadata columns). This is cosmetic:
#' the model fit (`prep_data_new()` uses only `id`/`antigen_iso`/`timeindays`/
#' `result`) and the figures are unaffected, and keeping the columns is what
#' makes the serotype / age subsets reproducible from one function.
#'
#' @param data Compiled MFI table (the `Compiled` sheet): one row per
#'   sample x isotype, containing `study_name`, `sid`, `timepoint`,
#'   `Actual day`, `isotype_name`, `cohort_name`, `age`, and the antigen MFI
#'   columns.
#' @param antigen <[`data-masking`][dplyr::dplyr_data_masking]> The antigen MFI
#'   column to extract, e.g. `n_ipab_MFI`.
#' @param study_filter Study to keep (default `"SOSAR"`, the Dhaka cohort).
#' @param as_case If `TRUE` (default) return a `serocalculator`
#'   case-data object;
#'   if `FALSE` return a plain tibble (useful for testing the reshape alone).
#'
#' @return Long case data with columns `id`, `antigen_iso`, `visit`,
#'   `timeindays`, `result`, `visit_num`, `cohort_name`, `age` (plus `sid`,
#'   `Actual day`, and the antigen MFI column). `timeindays` equals `Actual day`
#'   at this stage; it is shifted later by [adjust_timeindays_durdia()] for the
#'   `_new` datasets.
#' @export
process_shigella_data <- function(data,
                                  antigen,
                                  study_filter = "SOSAR",
                                  as_case = TRUE) {
  antigen_col <- rlang::ensym(antigen)

  out <- data |>
    dplyr::filter(.data$study_name == study_filter) |>
    dplyr::select(
      "isotype_name", "sid", "timepoint", "Actual day",
      "cohort_name", "age", !!antigen_col
    ) |>
    dplyr::mutate(
      id          = .data$sid,
      antigen_iso = .data$isotype_name,
      visit       = .data$timepoint,
      timeindays  = .data[["Actual day"]],
      result      = !!antigen_col
    ) |>
    dplyr::group_by(.data$id, .data$antigen_iso) |>
    dplyr::arrange(.data$visit, .by_group = TRUE) |>
    dplyr::mutate(visit_num = rank(.data$visit, ties.method = "first")) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$timeindays))

  if (!as_case) {
    return(out)
  }

  # as_case_data() is exported by serodynamics
  serodynamics::as_case_data(
    out,
    id_var        = "id",
    biomarker_var = "antigen_iso",
    time_in_days  = "timeindays",
    value_var     = "result"
  )
}
