#' Observed individual series (long) for grey background trajectories
#' @param raw_dataset Case/raw data.
#' @return Tibble `id, t, value, iso`.
#' @export
extract_individual_obs <- function(raw_dataset) {
  time_col <- if ("timeindays" %in% names(raw_dataset)) "timeindays" else "Actual day" # nolint: line_length_linter.
  value_col <- if ("result" %in% names(raw_dataset)) {
    "result"
  } else {
    mfi <- names(raw_dataset)[grepl("_MFI$", names(raw_dataset))]
    if (length(mfi) == 1) mfi else "result"
  }
  id_col <- if ("sid" %in% names(raw_dataset)) "sid" else "id"
  iso_col <- if ("isotype_name" %in% names(raw_dataset)) {
    "isotype_name"
  } else if
  ("Iso_type" %in% names(raw_dataset)) {
    "Iso_type"
  } else {
    "antigen_iso"
  }

  raw_dataset |>
    dplyr::transmute(
      id = .data[[id_col]], t = .data[[time_col]],
      value = .data[[value_col]], iso = .data[[iso_col]]
    ) |>
    dplyr::filter(!is.na(.data$value), !is.na(.data$t)) |>
    dplyr::mutate(iso = factor(.data$iso, levels = c("IgG", "IgA")))
}
