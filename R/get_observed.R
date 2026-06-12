#' Observed antibody series for one subject x isotype
#'
#' Column-name tolerant accessor (handles `id`/`sid`, the several isotype and
#' time column names, and `result` vs an `*_MFI` column) used by the
#' individual-trajectory figures.
#'
#' @param dataset Case data or raw data frame.
#' @param sid Subject id.
#' @param iso Isotype (`"IgG"`/`"IgA"`).
#' @return Tibble `t, value` sorted by time.
#' @export
get_observed <- function(dataset, sid, iso) {
  id_col <- if ("id" %in% names(dataset)) "id" else if ("sid" %in% names(dataset)) "sid" else NA_character_ # nolint: line_length_linter.
  if (is.na(id_col)) cli::cli_abort("Dataset has no {.code id} or {.code sid} column.") # nolint: line_length_linter.

  iso_col <- if ("antigen_iso"  %in% names(dataset)) "antigen_iso" else
             if ("isotype_name" %in% names(dataset)) "isotype_name" else
             if ("Iso_type"     %in% names(dataset)) "Iso_type" else NA_character_
  if (is.na(iso_col)) cli::cli_abort("Dataset has no isotype column.")

  time_col <- if ("timeindays" %in% names(dataset)) "timeindays" else
              if ("Actual day" %in% names(dataset)) "Actual day" else
              if ("timepoint"  %in% names(dataset)) "timepoint" else NA_character_
  if (is.na(time_col)) cli::cli_abort("Dataset has no time column.")

  val_col <- if ("result" %in% names(dataset)) {
    "result"
  } else {
    mfi <- names(dataset)[grepl("_MFI$", names(dataset))]
    if (length(mfi) == 1) mfi else cli::cli_abort("Dataset has no {.code result} or {.code _MFI} column.") # nolint: line_length_linter.
  }

  out <- dataset |>
    dplyr::filter(.data[[id_col]] == .env$sid, .data[[iso_col]] == .env$iso) |>
    dplyr::transmute(t = .data[[time_col]], value = .data[[val_col]]) |>
    dplyr::arrange(.data$t) |>
    dplyr::filter(!is.na(.data$value))

  if (nrow(out) == 0) {
    cli::cli_warn("No observed data for sid={.val {sid}}, iso={.val {iso}}.")
  }
  out
}
