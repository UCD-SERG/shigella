#' Shift observation times by pre-presentation diarrhea duration ("_new")
#'
#' Implements the time correction used to build the `dL_*_new` objects. Each
#' participant's whole trajectory is shifted forward by the duration of diarrhea
#' that preceded hospital presentation, moving `t = 0` closer to the true onset
#' of infection.
#'
#' The transformation (validated against `dL_clean_sf2a` vs
#'   `dL_clean_sf2a_new`):
#' \enumerate{
#'   \item merge `DurDia_hours` onto the case data by participant id (constant
#'         per participant);
#'   \item `DurDia_days = DurDia_hours / 24`;
#'   \item `timeindays_new = timeindays + DurDia_days` (added to every visit);
#'   \item overwrite `timeindays` with `timeindays_new`.
#' }
#'
#' @param case_data Case data from [process_shigella_data()].
#' @param durdia A data frame with one row per participant containing the id and
#'   diarrhea-duration columns.
#' @param id_col Name of the participant id column in `case_data` (default
#'   `"id"`, the standardised id from [process_shigella_data()]; `"sid"` is
#'   equivalent since the two columns are identical).
#' @param durdia_id_col Name of the id column in `durdia`. Default `"CaseID"`,
#'   confirmed against the "Duration of symptoms" sheet.
#' @param hours_col Name of the diarrhea-duration-in-hours column in `durdia`
#'   (default `"DurDia_hours"`, confirmed).
#'
#' @return `case_data` with added columns `DurDia_hours`, `DurDia_days`,
#'   `timeindays_new`, and with `timeindays` overwritten by `timeindays_new`.
#'   Case-data attributes are preserved.
#' @export
adjust_timeindays_durdia <- function(case_data,
                                     durdia,
                                     id_col = "id",
                                     durdia_id_col = "CaseID",
                                     hours_col = "DurDia_hours") {

  durdia_small <- durdia |>
    dplyr::transmute(
      !!id_col := .data[[durdia_id_col]],
      DurDia_hours = suppressWarnings(as.numeric(.data[[hours_col]]))
    ) |>
    dplyr::distinct()

  # Guard: DurDia must be a single value per participant.
  dup <- durdia_small |>
    dplyr::count(.data[[id_col]]) |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup) > 0) {
    cli::cli_abort(
      "{.arg durdia} has multiple {.field {hours_col}} values for {nrow(dup)} \\
       participant{?s}; it must be constant per participant."
    )
  }

  out <- case_data |>
    dplyr::left_join(durdia_small, by = id_col) |>
    dplyr::mutate(
      DurDia_days    = .data$DurDia_hours / 24,
      timeindays_new = .data$timeindays + .data$DurDia_days,
      timeindays     = .data$timeindays_new
    )

  .reattach_case_attrs(out, case_data) # nolint: object_usage_linter
}
