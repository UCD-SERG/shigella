#' Subset case data to an age group
#'
#' @param case_data Case data carrying an `age` column.
#' @param group `"under5"` keeps `age < 5`; `"plus5"` keeps `age > 5`.
#'   **Both groups exclude `age == 5`** (children aged exactly 5 are in neither
#'   stratum), matching how the shipped `_under5` / `_plus5`
#'   datasets were built.
#' @return Subset case data with attributes preserved.
#' @export
subset_age_group <- function(case_data, group = c("under5", "plus5")) {
  group <- match.arg(group)
  if (!"age" %in% names(case_data)) {
    cli::cli_abort(
      "{.arg case_data} must contain an {.field age} column; \\
       ensure {.fn process_shigella_data} kept it (it selects it from the \\
       compiled sheet)."
    )
  }
  n_boundary <- dplyr::n_distinct(
    dplyr::filter(case_data, .data$age == 5)$id
  )
  if (n_boundary > 0) {
    cli::cli_inform(
      "{n_boundary} participant{?s} with age == 5 excluded from both age
      strata."
    )
  }
  out <- if (group == "under5") {
    dplyr::filter(case_data, .data$age < 5)
  } else {
    dplyr::filter(case_data, .data$age > 5)
  }
  .reattach_case_attrs(out, case_data) # nolint: object_usage_linter
}
