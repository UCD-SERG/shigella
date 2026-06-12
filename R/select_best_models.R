# Internal helpers

# Internal: among models sharing individuals, which wins most per-individual?
#' @keywords internal
#' @noRd
.common_id_winners <- function(mae_long) {
  models <- unique(mae_long$model)
  if (length(models) == 1) return(models[1])
  ids_list <- mae_long |>
    dplyr::group_by(.data$model) |>
    dplyr::summarise(ids = list(unique(.data$sid)), .groups = "drop")
  common <- Reduce(intersect, ids_list$ids)
  if (length(common) == 0) return(models[1])
  mae_long |>
    dplyr::filter(.data$sid %in% common) |>
    dplyr::group_by(.data$sid) |>
    dplyr::slice_min(.data$mae, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::count(.data$model) |>
    dplyr::arrange(dplyr::desc(.data$n)) |>
    dplyr::slice(1) |>
    dplyr::pull(.data$model)
}

# Main function

#' Pick the best-performing model per antigen-isotype (lowest per-individual MAE)
#'
#' @param mae_overall,mae_serospec,mae_combined `get_mae()` tibbles
#'   (combine the per-antigen calls upstream). `NULL` for an unavailable
#'   model class.
#' @return Tibble `antigen, Iso_type, best_model`.
#' @export
select_best_models <- function(mae_overall, mae_serospec = NULL, mae_combined = NULL) { # nolint: line_length_linter.
  dplyr::bind_rows(
    dplyr::mutate(mae_overall,  model = "Overall"),
    if (!is.null(mae_serospec)) dplyr::mutate(mae_serospec, model = "Serotype-specific"), # nolint: line_length_linter.
    if (!is.null(mae_combined)) dplyr::mutate(mae_combined, model = "Combined flexneri") # nolint: line_length_linter.
  ) |>
    dplyr::group_by(.data$antigen, .data$Iso_type) |>
    dplyr::group_modify(~ tibble::tibble(best_model = .common_id_winners(.x))) |> # nolint: line_length_linter.
    dplyr::ungroup()
}
