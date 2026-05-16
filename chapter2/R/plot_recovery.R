#' Recovery plot for Phase 2 simulation study
#'
#' @param results_df data.frame of successful simulation results
#'   (already filtered to status == "OK")
#' @return ggplot object
plot_recovery <- function(results_df) {
  ggplot2::ggplot(results_df,
                  ggplot2::aes(x = true_rho_B, y = est_rho_B_median,
                               color = scenario)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = est_rho_B_lo,
                                        ymax = est_rho_B_hi),
                           width = 0.02, alpha = 0.4) +
    ggplot2::geom_jitter(width = 0.01, height = 0, size = 2.5, alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("A" = "#2166AC",
                                           "B" = "#92C5DE",
                                           "C" = "#B2182B")) +
    ggplot2::labs(
      title    = "Phase 2: Sigma_B Parameter Recovery",
      subtitle = "Each point = 1 simulation replicate, 95% CrI as error bars",
      x        = "True rho_B",
      y        = "Estimated rho_B (posterior median)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::facet_wrap(~ scenario, scales = "fixed")
}
