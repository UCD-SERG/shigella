# ── Internal helpers ────────────────────────────────

# One isotype panel for a Figure-4 comparison row: the subject's observed points
# plus overall / sero-specific / (optional) combined posterior-median bands.
#' @keywords internal
#' @noRd
.fig4_isotype_panel <- function(iso, sid, model_overall, model_sero,
                                data_overall, data_sero, model_combined,
                                antigen_label, xlim, times) {
  obs <- get_observed(data_sero, sid, iso)
  if (nrow(obs) == 0) obs <- get_observed(data_overall, sid, iso)

  pred_ovr  <- get_prediction_summary(model_overall, sid, iso, times)
  pred_sero <- get_prediction_summary(model_sero, sid, iso, times)

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = pred_ovr,
      ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
      fill = .model_colours[["Overall"]], alpha = 0.12) + # nolint: object_usage_linter
    ggplot2::geom_line(data = pred_ovr,
      ggplot2::aes(x = .data$t, y = .data$med,
                   color = "Overall", linetype = "Overall"), linewidth = 0.9) +
    ggplot2::geom_ribbon(data = pred_sero,
      ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
      fill = .model_colours[["Serotype-specific"]], alpha = 0.12) + # nolint: object_usage_linter
    ggplot2::geom_line(data = pred_sero,
      ggplot2::aes(x = .data$t, y = .data$med,
                   color = "Serotype-specific", linetype = "Serotype-specific"),
      linewidth = 0.9)

  if (!is.null(model_combined)) {
    pred_cmb <- get_prediction_summary(model_combined, sid, iso, times)
    if (nrow(pred_cmb) > 0) {
      p <- p +
        ggplot2::geom_ribbon(data = pred_cmb,
          ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
          fill = .model_colours[["Combined flexneri"]], alpha = 0.10) + # nolint: object_usage_linter
        ggplot2::geom_line(data = pred_cmb,
          ggplot2::aes(x = .data$t, y = .data$med,
                       color = "Combined flexneri", linetype = "Combined flexneri"),
          linewidth = 0.9)
    }
  }

  if (nrow(obs) > 0) {
    p <- p + ggplot2::geom_point(data = obs,
      ggplot2::aes(x = .data$t, y = .data$value),
      size = 2.8, shape = 21, fill = "white", color = "black", stroke = 0.8)
  }

  p +
    ggplot2::scale_color_manual(name = "Model", values = .model_colours, # nolint: object_usage_linter
                                breaks = names(.model_colours)) + # nolint: object_usage_linter
    ggplot2::scale_linetype_manual(name = "Model", values = .model_linetypes, # nolint: object_usage_linter
                                   breaks = names(.model_linetypes)) + # nolint: object_usage_linter
    ggplot2::scale_y_log10(labels = scales::label_comma()) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = "Days since symptom onset", y = "MFI (log scale)",
                  title = paste0(antigen_label, " - ", iso, " - ", sid)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 10, hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey20"),
      strip.text       = ggplot2::element_text(colour = "white", face = "bold"),
      legend.position  = "none")
}

# ── Main function ─────────────────────────────────

#' One IgG|IgA model-comparison row for a single subject
#'
#' Overlays overall (blue solid), serotype-specific (red dashed) and optional
#' combined-flexneri (green dotted) posterior-median trajectories with 95% CrI
#' ribbons and the subject's observed points. Each panel is built by
#' `.fig4_isotype_panel()`.
#'
#' @param sid Subject id.
#' @param model_overall,model_sero Fitted `sr_model` objects.
#' @param data_overall,data_sero Case data for observed-point lookup.
#' @param model_combined Optional combined-flexneri `sr_model` (`NULL` to omit).
#' @param row_label Row title (e.g. "A) S. flexneri 2a (n = 17)").
#' @param antigen_label Antigen label shown in each panel title.
#' @param xlim,times x-axis limits and prediction time grid.
#' @return A patchwork row (IgG | IgA).
#' @export
build_figure4_row <- function(sid, model_overall, model_sero,
                              data_overall, data_sero, model_combined = NULL,
                              row_label = "", antigen_label = "Sf2a",
                              xlim = c(0, 200), times = seq(0, 200, by = 1)) {
  panel <- function(iso) {
    .fig4_isotype_panel(iso, sid, model_overall, model_sero,
                        data_overall, data_sero, model_combined,
                        antigen_label, xlim, times)
  }
  (panel("IgG") | panel("IgA")) +
    patchwork::plot_annotation(
      title = row_label,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0)))
}
