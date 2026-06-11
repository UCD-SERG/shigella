# ── Internal helpers ────────────────────────────────

# Simulate `n_sim` replicate datasets from the joint posterior and return a long
# frame of simulated + observed assay values. For each retained draw the fitted
# mean is serodynamics:::ab(); a simulated value is exp(rnorm(1, log(mu_hat), sd)).
#' @keywords internal
#' @noRd
.ppc_simulate <- function(data, raw_dat, n_sim) {
  mod_prec_logy <- attr(data, "population_params") |>
    dplyr::filter(.data$Population_Parameter == "prec.logy") |>
    dplyr::select("Iteration", "Chain", "value", "Iso_type", "Stratification") |>
    dplyr::rename(prec_logy = "value")

  obs_dat <- raw_dat |>
    serodynamics:::use_att_names() |>
    dplyr::select("Subject", "Iso_type", "t", "result")

  fitted_dat <- data |>
    tidyr::pivot_wider(names_from = "Parameter", values_from = "value") |>
    dplyr::right_join(obs_dat, by = c("Subject", "Iso_type")) |>
    dplyr::mutate(mu_hat = serodynamics:::ab(.data$t, .data$y0, .data$y1,
                                             .data$t1, .data$alpha, .data$shape))

  sim_tab <- purrr::map_dfr(seq_len(n_sim), function(j) {
    dplyr::slice_sample(fitted_dat, by = c("Subject", "Iso_type", "t")) |>
      dplyr::left_join(mod_prec_logy,
                       by = c("Iteration", "Chain", "Iso_type", "Stratification")) |>
      dplyr::mutate(sd = 1 / sqrt(.data$prec_logy)) |>
      dplyr::rowwise() |>
      dplyr::mutate(value = exp(stats::rnorm(1, mean = log(.data$mu_hat), sd = .data$sd))) |>
      dplyr::ungroup() |>
      dplyr::select("Iso_type", "value") |>
      dplyr::mutate(estimate = "simulated", simulation = as.character(j))
  })

  obs_prep <- obs_dat |>
    dplyr::select("Iso_type", "result") |>
    dplyr::rename(value = "result") |>
    dplyr::mutate(estimate = "observed", simulation = "observed")

  dplyr::bind_rows(sim_tab, obs_prep)
}

# Overlay simulated vs observed densities (optionally faceted by isotype).
#' @keywords internal
#' @noRd
.ppc_density_plot <- function(tab, by_antigen, antigen_name) {
  title <- if (!is.null(antigen_name)) {
    paste("Posterior Predictive Check \u2014", antigen_name)
  } else {
    "Posterior Predictive Check"
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = tab,
      ggplot2::aes(x = .data$value, color = .data$estimate, group = .data$simulation,
                   linewidth = .data$estimate, alpha = .data$estimate), fill = NA) +
    ggplot2::scale_color_manual(values = c(observed = "black", simulated = "dodgerblue")) +
    ggplot2::scale_linewidth_manual(values = c(observed = 0.7, simulated = 0.3)) +
    ggplot2::scale_alpha_manual(values = c(observed = 0.7, simulated = 0.4)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = title, x = "Assay value")
  if (by_antigen) p <- p + ggplot2::facet_wrap(~ Iso_type)
  p
}

# ── Main function ─────────────────────────────────

#' S5 Fig: posterior predictive check for one fitted model
#'
#' Composes `.ppc_simulate()` (the Monte-Carlo replicate generation) and
#' `.ppc_density_plot()` (the overlaid densities).
#'
#' @param data A fitted `sr_model` (long draws with a `population_params`
#'   attribute), e.g. `overall_IpaB_pop_6`.
#' @param raw_dat The case data the model was fit to, e.g. `dL_clean_Ipab_new`.
#' @param by_antigen If `TRUE`, facet by isotype.
#' @param n_sim Number of simulated datasets (supplement uses 25).
#' @param antigen_name Optional antigen label for the title.
#' @param ... Unused (kept for call compatibility).
#' @return A ggplot.
#' @export
posterior_pred <- function(data = NULL, raw_dat = NULL, by_antigen = FALSE,
                           n_sim = 4, antigen_name = NULL, ...) {
  if (is.null(data)) cli::cli_abort("{.arg data} must be supplied.")
  if (is.null(raw_dat)) cli::cli_abort("{.arg raw_dat} must be supplied.")
  .ppc_simulate(data, raw_dat, n_sim) |>
    .ppc_density_plot(by_antigen = by_antigen, antigen_name = antigen_name)
}
