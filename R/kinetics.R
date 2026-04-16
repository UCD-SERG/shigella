#' Extract and Reshape Population-Level Kinetic Parameters from a Fitted Model
#'
#' Pulls the `"mu.par"` population-level draws from the `population_params`
#' attribute of a fitted serodynamics model object, pivots them to wide
#' format, and attaches an `antigen` label column.
#'
#' @param fit_object A fitted model object produced by the serodynamics
#'   hierarchical MCMC workflow. Must have an attribute `"population_params"`
#'   with columns `Population_Parameter`, `Iso_type`, `Parameter`, `value`,
#'   `Iteration`, and `Chain`.
#' @param antigen_label A character scalar used to label this antigen in the
#'   resulting data frame (e.g., `"IpaB"`, `"Sf2a"`).
#'
#' @return A tibble with one row per MCMC draw (chain × iteration), and
#'   columns `Iteration`, `Chain`, `antigen`, `Iso_type`, `y0`, `y1`, `t1`,
#'   `alpha`, `rho`. All parameter values are on the **log scale** as stored
#'   internally by the model; use [transform_pop_params()] to convert to the
#'   natural scale.
#'
#' @examples
#' \dontrun{
#' draws <- prep_pop_params(overall_IpaB_pop_6, "IpaB")
#' head(draws)
#' }
#'
#' @importFrom dplyr filter mutate select rename
#' @importFrom tidyr pivot_wider
#' @export
prep_pop_params <- function(fit_object, antigen_label) {
  pop_draws <- attr(fit_object, "population_params")

  pop_draws |>
    dplyr::filter(.data$Population_Parameter == "mu.par") |>
    dplyr::mutate(antigen = antigen_label) |>
    dplyr::select(
      "Iteration", "Chain", "antigen", "Iso_type", "Parameter", "value"
    ) |>
    tidyr::pivot_wider(names_from = "Parameter", values_from = "value") |>
    dplyr::rename(rho = "shape") |>
    dplyr::select("Iteration", "Chain", "antigen", "Iso_type",
                  "y0", "y1", "t1", "alpha", "rho")
}


#' Transform Population Parameters from Log Scale to Natural Scale
#'
#' Converts the log-scale MCMC draws produced by [prep_pop_params()] to
#' natural-scale kinetic parameters suitable for biological interpretation
#' and plotting.
#'
#' @param pop_params A tibble as returned by [prep_pop_params()], with
#'   columns `y0`, `y1`, `t1`, `alpha`, `rho` on the log scale.
#'
#' @return The same tibble with five additional columns:
#' \describe{
#'   \item{`y0_natural`}{Baseline antibody level (MFI): `exp(y0)`}
#'   \item{`y1_natural`}{Peak antibody level (MFI): `exp(y1) + exp(y0)`}
#'   \item{`t1_natural`}{Time to peak (days): `exp(t1)`}
#'   \item{`alpha_year`}{Decay rate (year\eqn{^{-1}}): `exp(alpha) * 365`}
#'   \item{`rho_natural`}{Decay shape: `exp(rho) + 1`}
#' }
#'
#' @details
#' The model stores parameters as:
#' \itemize{
#'   \item `y0`    = log(baseline)
#'   \item `y1`    = log(peak - baseline)  → peak = exp(y1) + exp(y0)
#'   \item `t1`    = log(time to peak in days)
#'   \item `alpha` = log(decay rate in days\eqn{^{-1}}) → convert to year\eqn{^{-1}}
#'   \item `rho`   = log(shape - 1)  → shape = exp(rho) + 1
#' }
#'
#' @examples
#' \dontrun{
#' draws <- prep_pop_params(overall_IpaB_pop_6, "IpaB")
#' draws_nat <- transform_pop_params(draws)
#' summary(draws_nat$y1_natural)
#' }
#'
#' @importFrom dplyr mutate
#' @export
transform_pop_params <- function(pop_params) {
  pop_params |>
    dplyr::mutate(
      y0_natural  = exp(.data$y0),
      y1_natural  = exp(.data$y1) + exp(.data$y0),
      t1_natural  = exp(.data$t1),
      alpha_year  = exp(.data$alpha) * 365,
      rho_natural = exp(.data$rho) + 1
    )
}


#' Summarise Population Kinetic Parameters to Posterior Medians and CrIs
#'
#' Aggregates MCMC draws (on the natural scale) to posterior median and 95%
#' credible intervals for each antigen–isotype combination.
#'
#' @param pop_params_natural A tibble as returned by [transform_pop_params()],
#'   containing `y0_natural`, `y1_natural`, `t1_natural`, `alpha_year`,
#'   `rho_natural`, `antigen`, and `Iso_type`.
#'
#' @return A summarised tibble with one row per antigen × isotype, containing
#'   posterior medians (`*_med`) and 2.5%/97.5% quantiles (`*_lo`, `*_hi`)
#'   for all five kinetic parameters.
#'
#' @examples
#' \dontrun{
#' draws_nat <- prep_pop_params(overall_IpaB_pop_6, "IpaB") |>
#'   transform_pop_params()
#' summarise_pop_params(draws_nat)
#' }
#'
#' @importFrom dplyr group_by summarise
#' @importFrom stats median quantile
#' @export
summarise_pop_params <- function(pop_params_natural) {
  pop_params_natural |>
    dplyr::group_by(.data$antigen, .data$Iso_type) |>
    dplyr::summarise(
      y0_med    = stats::median(.data$y0_natural,  na.rm = TRUE),
      y0_lo     = stats::quantile(.data$y0_natural,  0.025, na.rm = TRUE),
      y0_hi     = stats::quantile(.data$y0_natural,  0.975, na.rm = TRUE),
      y1_med    = stats::median(.data$y1_natural,  na.rm = TRUE),
      y1_lo     = stats::quantile(.data$y1_natural,  0.025, na.rm = TRUE),
      y1_hi     = stats::quantile(.data$y1_natural,  0.975, na.rm = TRUE),
      t1_med    = stats::median(.data$t1_natural,  na.rm = TRUE),
      t1_lo     = stats::quantile(.data$t1_natural,  0.025, na.rm = TRUE),
      t1_hi     = stats::quantile(.data$t1_natural,  0.975, na.rm = TRUE),
      alpha_med = stats::median(.data$alpha_year,  na.rm = TRUE),
      alpha_lo  = stats::quantile(.data$alpha_year,  0.025, na.rm = TRUE),
      alpha_hi  = stats::quantile(.data$alpha_year,  0.975, na.rm = TRUE),
      rho_med   = stats::median(.data$rho_natural, na.rm = TRUE),
      rho_lo    = stats::quantile(.data$rho_natural, 0.025, na.rm = TRUE),
      rho_hi    = stats::quantile(.data$rho_natural, 0.975, na.rm = TRUE),
      .groups   = "drop"
    )
}


#' Format a Posterior Median with 95% Credible Interval as a String
#'
#' Produces a formatted character string of the form
#' `"median (lower–upper)"`, suitable for tables and reports.
#'
#' @param med Numeric scalar or vector of posterior medians.
#' @param lo  Numeric scalar or vector of 2.5% credible-interval lower bounds.
#' @param hi  Numeric scalar or vector of 97.5% credible-interval upper bounds.
#' @param digits Integer. Number of decimal places to display. Default `2`.
#' @param sci Logical. If `TRUE`, format using scientific notation
#'   (`"e"` format). Default `FALSE`.
#'
#' @return A character vector of formatted strings, one element per input row.
#'
#' @examples
#' fmt_mci(1234.56, 890.1, 1678.9, digits = 2)
#' #> "1234.56 (890.10–1678.90)"
#'
#' fmt_mci(3.2e-9, 1.1e-10, 8.7e-8, digits = 2, sci = TRUE)
#' #> "3.20e-09 (1.10e-10–8.70e-08)"
#'
#' @export
fmt_mci <- function(med, lo, hi, digits = 2, sci = FALSE) {
  f <- function(x) {
    if (sci) {
      formatC(x, format = "e", digits = digits)
    } else {
      formatC(x, format = "f", digits = digits)
    }
  }
  sprintf("%s (%s\u2013%s)", f(med), f(lo), f(hi))
}


#' Format a Summary Table of Kinetic Parameters for Display
#'
#' Converts a summarised parameter tibble (from [summarise_pop_params()]) to a
#' display-ready format with `"Biomarker"` and formatted parameter columns.
#'
#' @param param_summary A tibble as returned by [summarise_pop_params()].
#'
#' @return A tibble with columns `Biomarker`, `y0`, `y1`, `t1`, `alpha`,
#'   `rho` (all character, formatted via [fmt_mci()]).
#'
#' @examples
#' \dontrun{
#' draws_nat <- prep_pop_params(overall_IpaB_pop_6, "IpaB") |>
#'   transform_pop_params()
#' tbl <- summarise_pop_params(draws_nat) |> format_param_table()
#' tbl
#' }
#'
#' @importFrom dplyr mutate select
#' @export
format_param_table <- function(param_summary) {
  param_summary |>
    dplyr::mutate(
      Biomarker = paste0(
        as.character(.data$antigen), "\u2013",
        as.character(.data$Iso_type)
      ),
      y0    = fmt_mci(.data$y0_med,    .data$y0_lo,    .data$y0_hi,    digits = 2),
      y1    = fmt_mci(.data$y1_med,    .data$y1_lo,    .data$y1_hi,    digits = 2),
      t1    = fmt_mci(.data$t1_med,    .data$t1_lo,    .data$t1_hi,    digits = 1),
      alpha = fmt_mci(.data$alpha_med, .data$alpha_lo, .data$alpha_hi, digits = 8),
      rho   = fmt_mci(.data$rho_med,   .data$rho_lo,   .data$rho_hi,   digits = 2)
    ) |>
    dplyr::select("Biomarker", "y0", "y1", "t1", "alpha", "rho")
}
