#' Shared model-comparison legend (drawn once at the bottom of the figure)
#' @return A grob.
#' @export
model_comparison_legend <- function() {
  legend_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = 1, y = 1,
      color = "Overall", linetype = "Overall")) +
    ggplot2::geom_line(ggplot2::aes(x = 1, y = 1,
      color = "Sero-specific", linetype = "Sero-specific")) +
    ggplot2::geom_line(ggplot2::aes(x = 1, y = 1,
      color = "Combined flexneri", linetype = "Combined flexneri")) +
    ggplot2::geom_point(ggplot2::aes(x = 1, y = 1, shape = "Observed"),
      size = 3, fill = "white", color = "black") +
    ggplot2::scale_color_manual(name = "Model", values = .model_colours, # nolint: object_usage_linter
                                breaks = names(.model_colours)) + # nolint: object_usage_linter
    ggplot2::scale_linetype_manual(name = "Model", values = .model_linetypes, # nolint: object_usage_linter
                                   breaks = names(.model_linetypes)) + # nolint: object_usage_linter
    ggplot2::scale_shape_manual(name = NULL, values = c("Observed" = 21)) +
    ggplot2::guides(
      color    = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
      linetype = ggplot2::guide_legend(order = 1),
      shape    = ggplot2::guide_legend(order = 2, override.aes = list(size = 3))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position  = "bottom", legend.box = "horizontal",
      legend.text      = ggplot2::element_text(size = 10),
      legend.title     = ggplot2::element_text(face = "bold", size = 10),
      legend.key.width = ggplot2::unit(2, "cm"))

  cowplot::get_legend(legend_plot)
}
