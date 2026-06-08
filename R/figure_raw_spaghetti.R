#' Raw longitudinal IgG/IgA trajectories, coloured by infecting serotype
#'
#' @param case_list Named list of the five **raw** single-antigen case datasets
#'   (i.e. unadjusted times), as in the main-text Figure 2.
#' @param numeric_y If `TRUE` (default, Figure 2) use friendly numeric y labels;
#'   if `FALSE` (S1/S2) use the plain log axis.
#' @return A ggplot.
#' @export
figure_raw_spaghetti <- function(case_list, numeric_y = TRUE) {
  dL_all <- stack_antigen_series(case_list)

  p <- ggplot2::ggplot(
    dL_all,
    ggplot2::aes(x = .data$timeindays, y = .data$result,
                 group = .data$sid, color = .data$infecting_serotype)) +
    ggplot2::geom_line(
      ggplot2::aes(linewidth = .data$infecting_serotype,
                   alpha = .data$infecting_serotype), lineend = "round") +
    ggplot2::geom_point(
      ggplot2::aes(alpha = .data$infecting_serotype),
      size = 1.2, shape = 16, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = serotype_palette()) +
    ggplot2::scale_linewidth_manual(values = serotype_linewidth(), guide = "none") +
    ggplot2::scale_alpha_manual(values = serotype_line_alpha(), guide = "none") +
    ggplot2::facet_grid(antigen ~ isotype_name, scales = "free_y") +
    ggplot2::labs(x = "Days since symptom onset",
                  y = "Antibody level (normalized MFI, log scale)",
                  color = "Infecting serotype") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey15", colour = "grey15"),
      strip.text       = ggplot2::element_text(colour = "white", face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(linewidth = 1.2, alpha = 0.9)))

  if (numeric_y) {
    p + ggplot2::scale_y_log10(
      breaks = c(10, 100, 1e3, 1e4, 1e5),
      labels = scales::comma_format(accuracy = 1))
  } else {
    p + ggplot2::scale_y_log10()
  }
}
