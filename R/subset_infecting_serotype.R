#' Subset case data to a set of infecting serotypes
#'
#' Keeps only participants whose infecting serotype (`cohort_name`) is in
#' `serotypes`. Used to build the serotype-specific datasets, e.g.
#' `serotypes = "Sf2a"` -> the n = 17 matched-infection set.
#'
#' @param case_data Case data carrying a `cohort_name` column.
#' @param serotypes Character vector of `cohort_name` values to keep.
#' @return Subset case data with attributes preserved.
#' @export
subset_infecting_serotype <- function(case_data, serotypes) {
  if (!"cohort_name" %in% names(case_data)) {
    cli::cli_abort(
      "{.arg case_data} must contain a {.field cohort_name} column; \\
       ensure {.fn process_shigella_data} kept it (it selects it from the \\
       compiled sheet)."
    )
  }
  out <- dplyr::filter(case_data, .data$cohort_name %in% serotypes)
  .reattach_case_attrs(out, case_data) # nolint: object_usage_linter
}
