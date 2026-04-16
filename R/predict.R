#' Predict Antibody Levels at Specified Times from Posterior Draws
#'
#' For one or more individuals and a single antigen-isotype combination,
#' evaluates the two-phase antibody kinetic model at each requested time
#' point across all MCMC posterior draws.
#'
#' @param model A fitted individual-level serodynamics model object, with
#'   columns `Subject`, `Iso_type`, `Chain`, `Iteration`, `Parameter`,
#'   and `value`.
#' @param ids Character vector of subject IDs to predict for.
#' @param antigen_iso A character scalar identifying the isotype to filter
#'   (e.g., `"IgG"`, `"IgA"`).
#' @param times Numeric vector of time points (in days since symptom onset)
#'   at which to evaluate the antibody trajectory. Default `seq(0, 200, 1)`.
#'
#' @return A long tibble with one row per draw x time point combination,
#'   including columns `id`, `t`, `Chain`, `Iteration`, and `res`
#'   (predicted antibody level on the natural MFI scale).
#'
#' @examples
#' \dontrun{
#' preds <- predict_posterior_at_times(
#'   model       = overall_IpaB_pop_6,
#'   ids         = "SOSAR-22008",
#'   antigen_iso = "IgG",
#'   times       = seq(0, 200, by = 5)
#' )
#' head(preds)
#' }
#'
#' @importFrom dplyr filter select mutate summarise arrange
#' @importFrom tidyr pivot_wider crossing pivot_longer
#' @export
predict_posterior_at_times <- function(model,
                                       ids,
                                       antigen_iso,
                                       times = seq(0, 200, by = 1)) {

  sr_model_sub <- model |>
    dplyr::filter(
      .data$Subject %in% ids,
      .data$Iso_type == antigen_iso
    )

  param_medians_wide <- sr_model_sub |>
    dplyr::select(
      "Chain", "Iteration", "Iso_type", "Parameter", "value", "Subject"
    ) |>
    tidyr::pivot_wider(
      names_from  = "Parameter",
      values_from = "value"
    ) |>
    dplyr::arrange(.data$Chain, .data$Iteration) |>
    dplyr::mutate(
      antigen_iso = factor(.data$Iso_type),
      id          = as.factor(.data$Subject),
      r           = .data$shape
    ) |>
    dplyr::select(-c("Iso_type", "Subject"))

  if (!"sample_id" %in% names(param_medians_wide)) {
    param_medians_wide <- param_medians_wide |>
      dplyr::mutate(sample_id = dplyr::row_number())
  }

  dt1 <- data.frame(t = times) |>
    dplyr::mutate(idx = dplyr::row_number()) |>
    tidyr::pivot_wider(
      names_from  = "idx",
      values_from = "t",
      names_prefix = "time"
    ) |>
    dplyr::slice(
      rep(seq_len(dplyr::n()), each = nrow(param_medians_wide))
    )

  cbind(param_medians_wide, dt1) |>
    tidyr::pivot_longer(
      cols      = dplyr::starts_with("time"),
      values_to = "t"
    ) |>
    dplyr::select(-"name") |>
    dplyr::mutate(
      res = .ab(
        .data$t,
        .data$y0,
        .data$y1,
        .data$t1,
        .data$alpha,
        .data$shape
      )
    )
}


#' Summarise Posterior Predictions for a Single Subject and Isotype
#'
#' A lightweight wrapper around [predict_posterior_at_times()] that returns
#' the posterior median and 95% credible interval of predicted antibody
#' levels at each time point.
#'
#' @param model_output A fitted model object (same format as `model` in
#'   [predict_posterior_at_times()]).
#' @param sid A single character subject ID.
#' @param iso A character scalar isotype (`"IgG"` or `"IgA"`).
#' @param times Numeric vector of prediction times (days). Default
#'   `seq(0, 200, by = 1)`.
#'
#' @return A tibble with columns `t`, `med`, `lo`, `hi` (posterior median,
#'   2.5%, and 97.5% quantiles of predicted MFI at each time point).
#'   Returns an empty tibble with a warning if no data are found.
#'
#' @examples
#' \dontrun{
#' pred <- get_prediction_summary(overall_IpaB_pop_6, "SOSAR-22008", "IgG")
#' }
#'
#' @importFrom dplyr filter select mutate
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats median quantile
#' @export
get_prediction_summary <- function(model_output,
                                   sid,
                                   iso,
                                   times = seq(0, 200, by = 1)) {

  sub <- model_output |>
    dplyr::filter(.data$Subject == sid, .data$Iso_type == iso)

  if (nrow(sub) == 0) {
    warning(paste("No data for subject", sid, "isotype", iso))
    return(tibble::tibble(
      t = numeric(), med = numeric(), lo = numeric(), hi = numeric()
    ))
  }

  wide <- sub |>
    dplyr::select("Iteration", "Chain", "Parameter", "value") |>
    tidyr::pivot_wider(names_from = "Parameter", values_from = "value") |>
    dplyr::mutate(sample_id = dplyr::row_number())

  y0_vec    <- wide$y0
  y1_vec    <- wide$y1
  t1_vec    <- wide$t1
  alpha_vec <- wide$alpha
  shape_vec <- wide$shape
  n_samples <- length(y0_vec)

  pred_mat <- matrix(NA_real_, nrow = n_samples, ncol = length(times))
  for (j in seq_along(times)) {
    pred_mat[, j] <- .ab(
      times[j], y0_vec, y1_vec, t1_vec, alpha_vec, shape_vec
    )
  }

  tibble::tibble(
    t   = times,
    med = apply(pred_mat, 2, stats::median, na.rm = TRUE),
    lo  = apply(pred_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    hi  = apply(pred_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )
}
