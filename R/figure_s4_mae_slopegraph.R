# ── Internal helpers ────────────────────────────────

# Step 1: stack the three MAE tables into one long frame (biomarker + model).
#' @keywords internal
#' @noRd
.mae_slope_long <- function(mae_overall, mae_serospec, mae_combined) {
  dplyr::bind_rows(
    dplyr::mutate(mae_overall,  model = "Overall"),
    dplyr::mutate(mae_serospec, model = "Serotype-specific"),
    dplyr::mutate(mae_combined, model = "Combined")
  ) |>
    dplyr::mutate(
      biomarker = paste0(.data$antigen, " ", .data$Iso_type),
      model = factor(.data$model, levels = c("Overall", "Serotype-specific", "Combined")) # nolint: line_length_linter.
    )
}

# Step 2: keep only biomarkers with >1 model, and within each only the
# individuals shared across all of that biomarker's models.
#' @keywords internal
#' @noRd
.mae_slope_shared <- function(mae_slope) {
  multi <- mae_slope |>
    dplyr::group_by(.data$biomarker) |>
    dplyr::summarise(n_models = dplyr::n_distinct(.data$model), .groups = "drop") |> # nolint: line_length_linter.
    dplyr::filter(.data$n_models > 1) |>
    dplyr::pull(.data$biomarker)

  mae_slope |>
    dplyr::filter(.data$biomarker %in% multi) |>
    dplyr::group_by(.data$biomarker) |>
    dplyr::group_modify(~{
      ids_by_model <- .x |>
        dplyr::group_by(.data$model) |>
        dplyr::summarise(ids = list(unique(.data$sid)), .groups = "drop")
      common <- Reduce(intersect, ids_by_model$ids)
      dplyr::filter(.x, .data$sid %in% common)
    }) |>
    dplyr::ungroup()
}

# Step 3: per individual, flag whether the best alternative beat overall; tidy labels. # nolint: line_length_linter.
#' @keywords internal
#' @noRd
.mae_slope_flagged <- function(shared) {
  shared |>
    dplyr::group_by(.data$biomarker, .data$sid) |>
    dplyr::mutate(
      overall_mae  = .data$mae[.data$model == "Overall"][1],
      best_alt_mae = min(.data$mae[.data$model != "Overall"], na.rm = TRUE),
      improved = dplyr::if_else(.data$best_alt_mae < .data$overall_mae,
                                "Alternative better", "Overall better")
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      biomarker = dplyr::recode(.data$biomarker,
        "sf2a_osp IgA"   = "Sf2a IgA",   "sf2a_osp IgG"   = "Sf2a IgG",
        "sf3a_osp IgA"   = "Sf3a IgA",   "sf3a_osp IgG"   = "Sf3a IgG",
        "sonnei_osp IgA" = "Sonnei IgA", "sonnei_osp IgG" = "Sonnei IgG")
    )
}

# Step 4: the slopegraph itself.
#' @keywords internal
#' @noRd
.mae_slope_plot <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$model, y = .data$mae, group = .data$sid)) + # nolint: line_length_linter.
    ggplot2::geom_line(ggplot2::aes(color = .data$improved), alpha = 0.6, linewidth = 0.45) + # nolint: line_length_linter.
    ggplot2::geom_point(
      ggplot2::aes(shape = .data$model, fill = .data$model),
      size = 2.4, stroke = 0.2, color = "black", alpha = 0.9) +
    ggplot2::scale_color_manual(
      values = c("Overall better" = "grey70", "Alternative better" = "#5BBE7A"),
      name = "Lower MAE") +
    ggplot2::scale_fill_manual(
      values = c("Overall" = "#3B82F6", "Serotype-specific" = "#EF4444",
                 "Combined" = "#10B981"), name = "Model") +
    ggplot2::scale_shape_manual(
      values = c("Overall" = 22, "Serotype-specific" = 21, "Combined" = 24), name = "Model") + # nolint: line_length_linter.
    ggplot2::facet_wrap(~ biomarker, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = "Per-individual MAE across modeling approaches",
      subtitle = "Each line connects the same individual across available models within each antigen-isotype panel", # nolint: line_length_linter.
      x = "Modeling approach", y = "Per-individual MAE (log-MFI scale)") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey20", color = "grey20"), # nolint: line_length_linter.
      strip.text = ggplot2::element_text(colour = "white", face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "bottom", legend.box = "vertical",
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

# ── Main function ─────────────────────────────────

#' S4 Fig: per-individual MAE across modeling approaches (slopegraph)
#'
#' Composes `.mae_slope_long()` -> `.mae_slope_shared()` -> `.mae_slope_flagged()`
#' -> `.mae_slope_plot()`.
#'
#' @param mae_overall,mae_serospec,mae_combined `get_mae()` tibbles
#'   (`sid, antigen, Iso_type, mae`) per model class.
#' @return A ggplot.
#' @export
figure_s4_mae_slopegraph <- function(mae_overall, mae_serospec, mae_combined) {
  .mae_slope_long(mae_overall, mae_serospec, mae_combined) |>
    .mae_slope_shared() |>
    .mae_slope_flagged() |>
    .mae_slope_plot()
}
