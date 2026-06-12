# Internal helpers

# Reshape one model's wide parameter summary into long (one row per parameter).
#' @keywords internal
#' @noRd
.forest_long <- function(pop_table, model_label) {
  spec <- tibble::tribble(
    ~param,                 ~order, ~med,        ~lo,        ~hi,
    "y[0]~(baseline)",      1L,     "y0_med",    "y0_lo",    "y0_hi",
    "y[1]~(peak)",          2L,     "y1_med",    "y1_lo",    "y1_hi",
    "t[1]~(days~to~peak)",  3L,     "t1_med",    "t1_lo",    "t1_hi",
    "alpha~(decay~rate)",   4L,     "alpha_med", "alpha_lo", "alpha_hi",
    "rho~(shape)",          5L,     "rho_med",   "rho_lo",   "rho_hi"
  )
  purrr::pmap_dfr(spec, function(param, order, med, lo, hi) {
    pop_table |>
      dplyr::transmute(
        .data$antigen, .data$Iso_type,
        med = .data[[med]], lo = .data[[lo]], hi = .data[[hi]],
        param = param, param_order = order
      )
  }) |>
    dplyr::mutate(model = model_label)
}

# Stack the (2 or 3) model summaries into one long, factor-ordered frame.
#' @keywords internal
#' @noRd
.forest_params_long <- function(pop_overall, pop_serospec, pop_combined = NULL) { # nolint: line_length_linter.
  out <- dplyr::bind_rows(
    .forest_long(pop_overall, "Overall"),
    .forest_long(pop_serospec, "Serotype-specific"),
    if (!is.null(pop_combined)) .forest_long(pop_combined, "Combined flexneri")
  )
  out |>
    dplyr::mutate(
      antigen = factor(.data$antigen, levels = c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei")), # nolint: line_length_linter.
      Iso_type = factor(.data$Iso_type, levels = c("IgG", "IgA")),
      model = factor(.data$model,
        levels = c("Overall", "Serotype-specific", "Combined flexneri")
      ), # nolint: line_length_linter.
      param = stats::reorder(.data$param, .data$param_order)
    )
}

# One parameter panel (antigen on y, faceted by isotype); optional log x-axis.
#' @keywords internal
#' @noRd
.forest_panel <- function(df_sub, colors, point_size, line_width, dodge_width,
                          use_log = TRUE) {
  p <- ggplot2::ggplot(
    df_sub,
    ggplot2::aes(x = .data$med, y = .data$antigen, colour = .data$model)
  ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
      linewidth = line_width,
      position = ggplot2::position_dodge(width = dodge_width), 
      show.legend = FALSE
    ) + # nolint: line_length_linter.
    ggplot2::geom_point(
      size = point_size,
      position = ggplot2::position_dodge(width = dodge_width)
    ) +
    ggplot2::scale_colour_manual(values = colors, name = NULL) +
    ggplot2::facet_wrap(~Iso_type, ncol = 2) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey20"),
      strip.text = ggplot2::element_text(colour = "white", face = "bold", 
                                         size = 12), # nolint: line_length_linter.
      axis.text.y = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 8, 4, 4)
    )
  if (use_log) {
    p <- p + ggplot2::scale_x_log10(labels = function(x) {
      ifelse(x >= 1,
        formatC(x, format = "g", digits = 3), formatC(x, format = "e", 
                                                      digits = 0)
      )
    }) # nolint: line_length_linter.
  }
  p
}

# Main function

#' S6 Fig: population kinetic-parameter forest plot (model comparison)
#'
#' Composes `.forest_params_long()` and per-parameter `.forest_panel()`s, then
#' stacks them with patchwork. y0, y1, t1, alpha are drawn on the log scale.
#'
#' @param pop_overall,pop_serospec `pop_param_summary()` tables for the overall
#'   and serotype-specific models.
#' @param pop_combined Optional `pop_param_summary()` table for the combined
#'   flexneri model (`NULL` to omit).
#' @param colors Named model -> colour vector.
#' @param point_size,line_width,dodge_width Styling controls.
#' @return A patchwork figure.
#' @export
plot_model_comparison_forest <- function(pop_overall, pop_serospec,
                                         pop_combined = NULL,
                                         colors = c(
                                           "Overall"           = "#2563EB",
                                           "Serotype-specific" = "#DC2626",
                                           "Combined flexneri" = "#16A34A"
                                         ),
                                         point_size = 3, line_width = 0.8,
                                         dodge_width = 0.5) {
  params_long <- .forest_params_long(pop_overall, pop_serospec, pop_combined)
  log_params <- c(1, 2, 3, 4) # y0, y1, t1, alpha on log scale

  panels <- params_long |>
    dplyr::group_by(.data$param, .data$param_order) |>
    dplyr::group_split() |>
    purrr::map(function(df_sub) {
      .forest_panel(df_sub, colors, point_size, line_width, dodge_width,
        use_log = unique(df_sub$param_order) %in% log_params
      ) +
        ggplot2::ggtitle(parse(text = unique(as.character(df_sub$param))))
    })
  panels[[length(panels)]] <- panels[[length(panels)]] +
    ggplot2::theme(legend.position = "bottom")

  patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(
      title = paste(
        "Population-level antibody kinetic parameters:",
        "Overall vs. Serotype-specific models"
      ),
      subtitle = paste(
        "Posterior median and 95% credible interval.",
        "Blue = overall (n=48); Red = serotype-specific;",
        "Green = combined S. flexneri (n=25, where applicable)"
      ),
      caption = paste(
        "Panels faceted by isotype (IgG | IgA).\n",
        "Parameters y\u2080, y\u2081, t\u2081, and \u03B1 plotted on log scale.\n", # nolint: line_length_linter.
        "IpaB and Sf6: overall model only (no serotype-specific model).\n",
        "MCMC: 10,000,000 iterations, 50,000 burn-in, 4 chains per model."
      ),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5), # nolint: line_length_linter.
        plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, colour = "grey40"), # nolint: line_length_linter.
        plot.caption  = ggplot2::element_text(size = 10, colour = "grey50", 
                                              hjust = 0.5)
      )
    ) # nolint: line_length_linter.
}
