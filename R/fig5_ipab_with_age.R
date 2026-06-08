# ── Internal helpers ────────────────────────────────

# Population mean curves for overall + the two age strata, factor-ordered.
#' @keywords internal
#' @noRd
.age_stratified_curves <- function(model_overall, model_under5, model_plus5, t_grid) {
  dplyr::bind_rows(
    dplyr::mutate(pop_mean_curve(model_overall, t_grid), group = "Overall (n=48)"),
    dplyr::mutate(pop_mean_curve(model_under5,  t_grid), group = "<5 years"),
    dplyr::mutate(pop_mean_curve(model_plus5,   t_grid), group = "\u22655 years")
  ) |>
    dplyr::mutate(group = factor(.data$group,
      levels = c("<5 years", "\u22655 years", "Overall (n=48)")))
}

# The Fig 5A panel: grey individual lines + age-group ribbons and median curves.
#' @keywords internal
#' @noRd
.fig5_age_plot <- function(all_curves, indiv_data, antigen_label, log_y, xlim) {
  color_vals  <- c("<5 years" = "#E64A19", "\u22655 years" = "#2E7D32",
                   "Overall (n=48)" = "#1f77b4")
  lt_vals     <- c("<5 years" = "dashed", "\u22655 years" = "dashed",
                   "Overall (n=48)" = "solid")
  alpha_vals  <- c("<5 years" = 0.10, "\u22655 years" = 0.10,
                   "Overall (n=48)" = 0.08)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = indiv_data,
      ggplot2::aes(x = .data$t, y = .data$value, group = .data$id),
      color = "grey75", alpha = 0.18, linewidth = 0.3)

  for (grp in names(color_vals)) {
    p <- p + ggplot2::geom_ribbon(
      data = dplyr::filter(all_curves, .data$group == grp),
      ggplot2::aes(x = .data$t, ymin = .data$lo, ymax = .data$hi),
      fill = color_vals[[grp]], alpha = alpha_vals[[grp]])
  }

  p <- p +
    ggplot2::geom_line(data = all_curves,
      ggplot2::aes(x = .data$t, y = .data$med,
                   color = .data$group, linetype = .data$group),
      linewidth = 1.0, lineend = "round") +
    ggplot2::scale_color_manual(values = color_vals, name = "Age group") +
    ggplot2::scale_linetype_manual(values = lt_vals, name = "Age group") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2)),
      linetype = ggplot2::guide_legend()) +
    ggplot2::facet_wrap(~ iso, ncol = 2) +
    ggplot2::labs(x = "Days since symptom onset", y = "Normalized MFI",
                  title = antigen_label) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      strip.text      = ggplot2::element_text(face = "bold", size = 12),
      axis.title      = ggplot2::element_text(size = 11),
      axis.text       = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(face = "bold", size = 10),
      legend.text     = ggplot2::element_text(size = 9),
      legend.key.width = ggplot2::unit(1.5, "cm"))

  if (log_y) p <- p + ggplot2::scale_y_log10(labels = scales::label_comma())
  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  p
}

# ── Main function ─────────────────────────────────

#' IpaB population-trajectory panel with age-stratified curves (Fig 5A)
#'
#' Builds the curves with [.age_stratified_curves()] and renders via
#' [.fig5_age_plot()] over grey individual lines.
#'
#' @param model_overall,model_under5,model_plus5 Fitted IpaB `sr_model`s.
#' @param raw_overall Case data for the grey individual lines.
#' @param antigen_label Panel title.
#' @param t_grid,log_y,xlim Styling controls.
#' @return A ggplot.
#' @export
fig5_ipab_with_age <- function(model_overall, model_under5, model_plus5,
                               raw_overall,
                               antigen_label = "A) IpaB (Overall + age-stratified)",
                               t_grid = seq(0, 210, by = 5),
                               log_y = TRUE, xlim = c(0, 210)) {
  all_curves <- .age_stratified_curves(model_overall, model_under5, model_plus5, t_grid)
  indiv_data <- extract_individual_obs(raw_overall)
  .fig5_age_plot(all_curves, indiv_data, antigen_label, log_y, xlim)
}
