# ── Internal helpers ────────────────────────────────

# Observed (id, t, obs, antigen_iso) for the requested ids + isotype; aborts if
# the filter leaves nothing.
#' @keywords internal
#' @noRd
.observed_long <- function(dataset, ids, antigen_iso) {
  time_var  <- dataset |> serodynamics:::get_timeindays_var()
  value_var <- dataset |> serocalculator::get_values_var()

  observed_data <- dataset |>
    dplyr::rename(t = {{ time_var }}, obs = {{ value_var }}) |>
    dplyr::select(dplyr::all_of(c("id", "t", "obs", "antigen_iso"))) |>
    dplyr::mutate(id = as.character(.data$id)) |>
    dplyr::filter(.data$id %in% .env$ids, .data$antigen_iso == .env$antigen_iso)

  if (nrow(observed_data) == 0) {
    cli::cli_abort(c(
      "No observed data found for the specified IDs and antigen_iso.",
      "i" = "IDs: {.val {ids}}",
      "i" = "antigen_iso: {.val {antigen_iso}}"
    ))
  }
  observed_data
}

# Posterior-median prediction at each observed time, joined to the observations,
# with residual / abs / squared columns on the chosen scale.
#' @keywords internal
#' @noRd
.pointwise_residuals <- function(model, observed_data, ids, antigen_iso, scale) { # nolint: line_length_linter.
  obs_times <- sort(unique(observed_data$t))

  pred_summary <- predict_posterior_at_times(
    model = model, ids = ids, antigen_iso = antigen_iso, times = obs_times
  ) |>
    dplyr::summarise(
      .by = dplyr::all_of(c("id", "t")),
      pred_med   = stats::median(.data$res, na.rm = TRUE),
      pred_lower = stats::quantile(.data$res, 0.025, na.rm = TRUE),
      pred_upper = stats::quantile(.data$res, 0.975, na.rm = TRUE)
    ) |>
    dplyr::mutate(id = as.character(.data$id))

  residual_data <- observed_data |>
    dplyr::inner_join(pred_summary, by = c("id", "t"))

  if (scale == "log") {
    n_nonpos <- sum(residual_data$obs <= 0) + sum(residual_data$pred_med <= 0)
    if (n_nonpos > 0) {
      cli::cli_warn(
        "Dropping {n_nonpos} non-positive observation{?s} before log-scale residuals." # nolint: line_length_linter.
      )
      residual_data <- dplyr::filter(residual_data, .data$obs > 0, .data$pred_med > 0) # nolint: line_length_linter.
    }
    residual_data |>
      dplyr::mutate(
        residual     = log(.data$obs) - log(.data$pred_med),
        abs_residual = abs(.data$residual),
        sq_residual  = .data$residual^2
      )
  } else {
    residual_data |>
      dplyr::mutate(
        residual     = .data$obs - .data$pred_med,
        abs_residual = abs(.data$residual),
        sq_residual  = .data$residual^2
      )
  }
}

# ── Main function ─────────────────────────────────

#' Per-observation / per-individual residual metrics on the log scale
#'
#' Builds observed values (`.observed_long()`), joins posterior-median
#' predictions and computes residuals (`.pointwise_residuals()`), then aggregates
#' according to `summary_level`.
#'
#' @param model Fitted `sr_model`.
#' @param dataset Case data the model was fit to.
#' @param ids Subject ids to include.
#' @param antigen_iso Isotype (`"IgG"`/`"IgA"`).
#' @param scale `"original"` or `"log"` (default `"log"`, as in the manuscript).
#' @param summary_level `"id_antigen"` (default), `"pointwise"`, `"antigen"`,
#'   or `"overall"`.
#' @return Tibble of residual metrics at the requested level (MAE, RMSE, ...).
#' @export
compute_residual_metrics <- function(model, dataset, ids, antigen_iso,
                                      scale = c("original", "log"),
                                      summary_level = c("id_antigen", "pointwise", # nolint: line_length_linter.
                                                        "antigen", "overall")) {
  scale <- match.arg(scale)
  summary_level <- match.arg(summary_level)

  observed_data <- .observed_long(dataset, ids, antigen_iso)
  residual_data <- .pointwise_residuals(model, observed_data, ids, antigen_iso, scale) # nolint: line_length_linter.

  if (summary_level == "pointwise") {
    return(dplyr::select(residual_data, dplyr::all_of(c(
      "id", "antigen_iso", "t", "obs", "pred_med",
      "pred_lower", "pred_upper", "residual", "abs_residual", "sq_residual"
    ))))
  }

  group <- switch(summary_level,
                  id_antigen = c("id", "antigen_iso"),
                  antigen    = "antigen_iso",
                  overall    = character(0))

  residual_data |>
    dplyr::summarise(
      .by = dplyr::all_of(group),
      MAE   = mean(.data$abs_residual, na.rm = TRUE),
      RMSE  = sqrt(mean(.data$sq_residual, na.rm = TRUE)),
      SSE   = sum(.data$sq_residual, na.rm = TRUE),
      n_obs = dplyr::n()
    )
}
