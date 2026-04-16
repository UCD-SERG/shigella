#' Posterior Predictive Check for a Fitted Serodynamics Model
#'
#' Simulates replicated datasets from the joint posterior distribution and
#' overlays their density on the observed data distribution.
#'
#' @param data A fitted individual-level model object with a
#'   `"population_params"` attribute containing `prec.logy` draws.
#' @param raw_dat A case_data object containing the observed measurements.
#' @param by_antigen Logical. If `TRUE` (default `FALSE`), facet by isotype.
#' @param n_sim Integer. Number of simulated datasets. Default `4`.
#' @param antigen_name Optional character scalar for the plot title.
#' @param ... Currently unused.
#'
#' @return A `ggplot` object. Blue = simulated; black = observed.
#'
#' @examples
#' \dontrun{
#' p <- posterior_pred(
#'   data         = overall_IpaB_pop_6,
#'   raw_dat      = dL_clean_Ipab_new,
#'   by_antigen   = TRUE,
#'   n_sim        = 25,
#'   antigen_name = "IpaB"
#' )
#' print(p)
#' }
#'
#' @importFrom dplyr filter select rename mutate left_join slice_sample
#'   bind_rows rowwise ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_density scale_color_manual
#'   scale_linewidth_manual scale_alpha_manual theme_bw theme
#'   element_blank labs scale_x_log10 facet_wrap
#' @importFrom stats rnorm
#' @export
posterior_pred <- function(data,
                           raw_dat,
                           by_antigen   = FALSE,
                           n_sim        = 4,
                           antigen_name = NULL,
                           ...) {

  # --- Extract measurement precision draws --------------------------------
  mod_prec_logy <- attr(data, "population_params") |>
    dplyr::filter(.data$Population_Parameter == "prec.logy") |>
    dplyr::select(
      "Iteration", "Chain", "value", "Iso_type", "Stratification"
    ) |>
    dplyr::rename(prec_logy = "value")

  # --- Pivot individual-level draws to wide -------------------------------
  mod_dat <- data |>
    tidyr::pivot_wider(names_from = "Parameter", values_from = "value")

  # --- Observed data (using internal helper) ------------------------------
  obs_dat <- .use_att_names(raw_dat) |>
    dplyr::select("Subject", "Iso_type", "t", "result")

  # --- Match draws to observations ----------------------------------------
  matched_dat <- dplyr::right_join(
    mod_dat, obs_dat,
    by = c("Subject", "Iso_type")
  )

  fitted_dat <- matched_dat |>
    dplyr::mutate(
      mu_hat = .ab(
        .data$t, .data$y0, .data$y1,
        .data$t1, .data$alpha, .data$shape
      )
    )

  # --- Simulate n_sim replicated datasets ---------------------------------
  prepare_plot_tab <- tibble::tibble(
    Iso_type   = character(),
    value      = numeric(),
    estimate   = character(),
    simulation = character()
  )

  for (j in seq_len(n_sim)) {
    sampled_posterior <- dplyr::slice_sample(
      fitted_dat,
      by = c("Subject", "Iso_type", "t")
    )

    sampled_posterior <- dplyr::left_join(
      sampled_posterior, mod_prec_logy,
      by = c("Iteration", "Chain", "Iso_type", "Stratification")
    ) |>
      dplyr::mutate(sd = 1 / sqrt(.data$prec_logy)) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        value = exp(stats::rnorm(1, mean = log(.data$mu_hat), sd = .data$sd))
      ) |>
      dplyr::ungroup() |>
      dplyr::select("Iso_type", "value") |>
      dplyr::mutate(estimate = "simulated", simulation = as.character(j))

    prepare_plot_tab <- dplyr::bind_rows(prepare_plot_tab, sampled_posterior)
  }

  # --- Observed density layer ---------------------------------------------
  obs_dat_prep <- obs_dat |>
    dplyr::select("Iso_type", "result") |>
    dplyr::rename(value = "result") |>
    dplyr::mutate(estimate = "observed", simulation = "observed")

  plot_dat <- dplyr::bind_rows(prepare_plot_tab, obs_dat_prep)

  # --- Build plot ---------------------------------------------------------
  plot_title <- if (!is.null(antigen_name)) {
    paste("Posterior Predictive Check \u2014", antigen_name)
  } else {
    "Posterior Predictive Check"
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = plot_dat,
      ggplot2::aes(
        x         = .data$value,
        color     = .data$estimate,
        group     = .data$simulation,
        linewidth = .data$estimate,
        alpha     = .data$estimate
      ),
      fill = NA
    ) +
    ggplot2::scale_color_manual(
      values = c(observed = "black", simulated = "dodgerblue")
    ) +
    ggplot2::scale_linewidth_manual(
      values = c(observed = 0.7, simulated = 0.3)
    ) +
    ggplot2::scale_alpha_manual(
      values = c(observed = 0.7, simulated = 0.4)
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::labs(title = plot_title, x = "Assay value (MFI)")

  if (by_antigen) {
    p <- p + ggplot2::facet_wrap(~ Iso_type)
  }

  p
}
