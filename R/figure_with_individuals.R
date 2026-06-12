#' Population-trajectory panel for a single antigen (B-E of Fig 5)
#'
#' @param model_obj Fitted `sr_model`.
#' @param raw_dataset Case data for the grey individual lines.
#' @param antigen_label Panel title.
#' @param t_grid,line_color,ribbon_alpha,individual_alpha,individual_color,log_y,xlim
#'   Styling controls.
#' @return A ggplot.
#' @export
figure_with_individuals <- function(model_obj, raw_dataset, antigen_label,
                                    t_grid = seq(0, 210, by = 5),
                                    line_color = "#1f77b4", ribbon_alpha = 0.15,
                                    individual_alpha = 0.12, individual_color = "grey70", # nolint: line_length_linter.
                                    log_y = TRUE, xlim = c(0, 210)) {
  pop_curves <- pop_mean_curve(model_obj, t_grid)
  indiv_data <- extract_individual_obs(raw_dataset)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = indiv_data,
      ggplot2::aes(x = .data$t, y = .data$value, group = .data$id),
      color = individual_color, alpha = individual_alpha, linewidth = 0.3) +
    ggplot2::geom_ribbon(data = pop_curves,
      ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
      fill = "grey50", alpha = ribbon_alpha) +
    ggplot2::geom_line(data = pop_curves,
      ggplot2::aes(x = .data$t, y = .data$med),
      color = line_color, linewidth = 1.0, lineend = "round") +
    ggplot2::facet_wrap(~ iso, ncol = 2) +
    ggplot2::labs(x = "Days since symptom onset", y = "Normalized MFI",
                  title = antigen_label) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      axis.title = ggplot2::element_text(size = 11),
      axis.text  = ggplot2::element_text(size = 10),
      legend.position = "none")

  if (log_y) p <- p + ggplot2::scale_y_log10(labels = scales::label_comma())
  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  p
}
