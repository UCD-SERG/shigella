utils::globalVariables(c("antigen", "iso"))

#' Summarize population-level ("newperson") antibody trajectories
#' from overall models
#'
#' Computes median and credible interval bands of the antibody trajectory for a
#' hypothetical new individual drawn from the population distribution
#' ("newperson").
#'
#' @param overall_models Named list of posterior draws in long format
#'   (one per antigen),
#'   with columns at least: Subject, Iso_type, Chain, Iteration,
#'   Parameter, value.
#' @param osps Character vector of antigen names
#'   (must match names in `overall_models`).
#' @param ids Character vector of Subject IDs to include (default: "newperson").
#' @param isotypes Character vector of isotypes
#'   (default: c("IgG", "IgA")).
#' @param t_grid Numeric vector of time points (days) to evaluate.
#' @param cred Credible level (default 0.95).
#' @param log_y Logical; if TRUE, applies log10 scale to y when returning plot.
#' @param xlim Optional numeric length-2 vector for x-axis limits when
#'   returning plot.
#' @param ylab Y-axis label for plot.
#' @param line_color Line color for plot.
#' @param ribbon_alpha Alpha for credible ribbon in plot.
#' @param facet_scales Passed to ggplot facet scales.
#' @param return_data If TRUE, returns a list with `plot` and `data`.
#'
#' @return By default, a ggplot. If `return_data = TRUE`, returns
#'   list(plot = p, data = df).
#'
#' @details
#' Requires that each model draw can be pivoted to wide parameters including:
#' y0, y1, t1, alpha, shape.
#'
#' @export
fig2_overall_newperson <- function(
  overall_models,
  osps = c("IpaB", "Sf2a", "Sf3a", "Sf6", "Sonnei"),
  ids = "newperson",
  isotypes = c("IgG", "IgA"),
  t_grid = seq(0, 210, by = 5),
  cred = 0.95,
  log_y = TRUE,
  xlim = c(0, 210),
  ylab = "Normalized MFI",
  line_color = "#1f77b4",
  ribbon_alpha = 0.20,
  facet_scales = "fixed",
  return_data = FALSE
) {
  q_lo <- (1 - cred) / 2
  q_hi <- 1 - q_lo

  get_sum_overall <- function(model_df, osp, iso) {
    model_df |>
      dplyr::filter(.data$Subject %in% ids, .data$Iso_type == iso) |>
      dplyr::select(.data$Chain, .data$Iteration, .data$Iso_type,
                    .data$Parameter, .data$value, .data$Subject) |>
      tidyr::pivot_wider(
        names_from = .data$Parameter,
        values_from = .data$value
      ) |>
      dplyr::mutate(
        antigen = osp,
        iso = factor(.data$Iso_type, levels = c("IgG", "IgA"))
      ) |>
      tidyr::crossing(t = t_grid) |>
      dplyr::mutate(
        res = ab(
          t = .data$t,
          y0 = .data$y0,
          y1 = .data$y1,
          t1 = .data$t1,
          alpha = .data$alpha,
          shape = .data$shape
        )
      ) |>
      dplyr::filter(is.finite(.data$res), .data$res > 0) |>
      dplyr::summarise(
        res.med  = stats::quantile(.data$res, 0.50, na.rm = TRUE),
        res.low  = stats::quantile(.data$res, q_lo, na.rm = TRUE),
        res.high = stats::quantile(.data$res, q_hi, na.rm = TRUE),
        .by = c("antigen", "iso", "t")
      )
  }

  sum_all <- lapply(osps, function(osp) {
    dplyr::bind_rows(lapply(isotypes, function(iso) {
      get_sum_overall(overall_models[[osp]], osp, iso)
    }))
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(antigen = factor(.data$antigen, levels = osps))

  p <- ggplot2::ggplot(sum_all, ggplot2::aes(x = .data$t)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$res.low, ymax = .data$res.high),
      alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$res.med),
      linewidth = 1.1,
      color = line_color,
      lineend = "round"
    ) +
    ggplot2::labs(x = "Days since fever onset", y = ylab) +
    ggplot2::facet_grid(rows = ggplot2::vars(antigen),
                        cols = ggplot2::vars(iso),
                        scales = facet_scales) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  if (log_y) p <- p + ggplot2::scale_y_log10()
  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)

  if (return_data) return(list(plot = p, data = sum_all))
  p
}
