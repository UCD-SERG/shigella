# Plotting helpers for Phase 1 empirical correlation analysis

#' Scatter grid: IgG vs IgA posterior medians, faceted by antigen × parameter
#'
#' @param all_scatter_df data.frame with columns IgG, IgA, Antigen, Parameter,
#'   panel_label (rho with CI string)
#' @return ggplot object
plot_param_scatter_grid <- function(all_scatter_df) {
  ggplot2::ggplot(all_scatter_df, ggplot2::aes(x = IgG, y = IgA)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.8) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "steelblue",
                         linewidth = 0.7, alpha = 0.15) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Antigen),
      cols = ggplot2::vars(Parameter),
      scales = "free",
      labeller = ggplot2::labeller(
        Parameter = c(y0 = "log(y0)", y1 = "log(y1)", t1 = "log(t1)",
                      alpha = "log(alpha)", shape = "log(shape-1)")
      )
    ) +
    ggplot2::geom_text(
      data = all_scatter_df |>
        dplyr::group_by(Antigen, Parameter) |>
        dplyr::slice(1),
      ggplot2::aes(label = panel_label),
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
      size = 2.8, fontface = "bold", color = "#B2182B"
    ) +
    ggplot2::labs(
      title    = "Panel B Supplement - IgG vs IgA posterior medians per subject",
      subtitle = "Each point = one individual; rho [95% CI] shown in each panel",
      x        = "IgG parameter (posterior median)",
      y        = "IgA parameter (posterior median)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey20"),
      strip.text       = ggplot2::element_text(color = "white", face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text        = ggplot2::element_text(size = 7),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Forest plot of parameter correlations with Fisher z CIs
#'
#' @param forest_df data.frame with columns rho, ci_fisher_lo, ci_fisher_hi,
#'   param_label, antigen_color
#' @return ggplot object
plot_param_forest <- function(forest_df) {
  ggplot2::ggplot(forest_df,
                  ggplot2::aes(x = rho, y = param_label,
                               color = antigen_color)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey50") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted",
                        color = "grey40") +
    ggplot2::geom_point(size = 3,
                        position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_fisher_lo, xmax = ci_fisher_hi),
      height   = 0.2,
      position = ggplot2::position_dodge(width = 0.5),
      linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = c("IpaB" = "#2166AC", "Sonnei" = "#4393C3",
                 "Sf2a" = "#92C5DE"),
      name   = "Antigen"
    ) +
    ggplot2::scale_x_continuous(limits = c(-0.5, 1),
                                breaks = seq(-0.5, 1, 0.25)) +
    ggplot2::labs(
      title    = "Parameter correlation 95% CI (Fisher z; subject-level, n_subj>=11)",
      subtitle = "Dashed = 0 (independence); dotted = 0.5 (Cohen large)",
      x        = "rho_parameter",
      y        = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title           = ggplot2::element_text(face = "bold"),
      legend.position      = "bottom",
      panel.grid.major.y   = ggplot2::element_blank()
    )
}
