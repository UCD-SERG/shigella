#' One isotype panel of the cross-reactivity heatmap
#'
#' @param data Output of [prep_heatmap_data()] filtered to one isotype.
#' @param iso_name Isotype label for the title/legend.
#' @param show_y_axis Whether to draw participant labels on the y-axis.
#' @param y_label_lookup Named vector mapping `row_id` -> display label.
#' @return A ggplot.
#' @export
make_heatmap <- function(data, iso_name, show_y_axis = TRUE,
                         y_label_lookup = attr(data, "y_label_lookup")) {
  ggplot2::ggplot(data,
    ggplot2::aes(x = .data$timepoint_label, y = .data$row_id,
                 fill = .data$result + 0.01)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.05) +
    ggplot2::facet_grid(serotype ~ antigen_clean, scales = "free_y",
                        space = "free_y", switch = "y") +
    ggplot2::scale_fill_viridis_c(option = "plasma",
                                  name = paste0(iso_name, "\nMFI"),
                                  labels = scales::label_number()) +
    ggplot2::scale_y_discrete(labels = y_label_lookup) +
    ggplot2::labs(title = iso_name, x = "Timepoint",
                  y = if (show_y_axis) "Participants (age)" else NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5), # nolint: line_length_linter.
      axis.text.x        = ggplot2::element_text(size = 8),
      axis.text.y        = if (show_y_axis) ggplot2::element_text(size = 7) else ggplot2::element_blank(), # nolint: line_length_linter.
      axis.ticks.y       = if (show_y_axis) ggplot2::element_line() else ggplot2::element_blank(), # nolint: line_length_linter.
      axis.title         = ggplot2::element_text(face = "bold"),
      strip.background.x = ggplot2::element_rect(fill = "grey15", colour = "grey15"), # nolint: line_length_linter.
      strip.background.y = ggplot2::element_rect(fill = "grey85", colour = "grey85"), # nolint: line_length_linter.
      strip.text.x       = ggplot2::element_text(colour = "white", face = "bold", size = 9), # nolint: line_length_linter.
      strip.text.y.left  = ggplot2::element_text(angle = 0, face = "bold", size = 9, hjust = 1), # nolint: line_length_linter.
      strip.placement    = "outside",
      panel.grid         = ggplot2::element_blank(),
      legend.position    = "right",
      panel.spacing.x    = ggplot2::unit(0.15, "lines"),
      panel.spacing.y    = ggplot2::unit(0.2, "lines")
    )
}
