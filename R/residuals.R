#' Extract Observed Antibody Data for a Single Subject and Isotype
#'
#' Retrieves the longitudinal observed MFI values for one participant and
#' one antigen-isotype combination from a dataset that may use varying
#' column name conventions.
#'
#' @param dataset A data frame or case_data object.
#' @param sid A single character subject ID.
#' @param iso A character scalar isotype filter (e.g., `"IgG"`).
#'
#' @return A tibble with columns `t` (time in days) and `value` (observed
#'   MFI), sorted by `t`, with missing values removed.
#'
#' @examples
#' \dontrun{
#' obs <- get_observed(dL_clean_Ipab, "SOSAR-22008", "IgG")
#' head(obs)
#' }
#'
#' @importFrom dplyr filter transmute arrange
#' @importFrom tibble tibble
#' @export
get_observed <- function(dataset, sid, iso) {


  id_col <- if ("id" %in% names(dataset)) "id" else
    if ("sid" %in% names(dataset)) "sid" else
      stop("Dataset has no 'id' or 'sid' column. Available: ",
           paste(names(dataset), collapse = ", "))

  iso_col <- if ("antigen_iso"   %in% names(dataset)) "antigen_iso"   else
    if ("isotype_name" %in% names(dataset)) "isotype_name" else
      if ("Iso_type"    %in% names(dataset)) "Iso_type"    else
        stop("Dataset has no isotype column. Available: ",
             paste(names(dataset), collapse = ", "))

  time_col <- if ("timeindays"  %in% names(dataset)) "timeindays"  else
    if ("Actual day"  %in% names(dataset)) "Actual day"  else
      if ("timepoint"   %in% names(dataset)) "timepoint"   else
        stop("Dataset has no time column.")

  val_col <- if ("result" %in% names(dataset)) {
    "result"
  } else {
    mfi_cols <- names(dataset)[grepl("_MFI$", names(dataset))]
    if (length(mfi_cols) == 1) mfi_cols else
      stop("No 'result' column and ambiguous MFI columns: ",
           paste(mfi_cols, collapse = ", "))
  }

  out <- dataset |>
    dplyr::filter(
      .data[[id_col]]  == sid,
      .data[[iso_col]] == iso
    ) |>
    dplyr::transmute(
      t     = .data[[time_col]],
      value = .data[[val_col]]
    ) |>
    dplyr::arrange(.data$t) |>
    dplyr::filter(!is.na(.data$value))

  if (nrow(out) == 0) {
    warning(sprintf(
      "No observed data for sid='%s', iso='%s'. Available IDs (first 5): %s",
      sid, iso,
      paste(head(unique(dataset[[id_col]]), 5), collapse = ", ")
    ))
  }

  out
}


#' Compute Residual-Based Prediction Accuracy Metrics
#'
#' For a fitted hierarchical model and a corresponding case dataset, computes
#' absolute and squared residuals between observed antibody measurements and
#' posterior median predictions.
#'
#' @param model A fitted model object (individual-level draws).
#' @param dataset A case_data object matching `model`.
#' @param ids Character vector of subject IDs to include.
#' @param antigen_iso A character scalar for the isotype to evaluate.
#' @param scale One of `"original"` or `"log"` (default `"log"`).
#' @param summary_level One of `"id_antigen"`, `"pointwise"`, `"antigen"`,
#'   or `"overall"`.
#'
#' @return A tibble whose rows depend on `summary_level`.
#'
#' @examples
#' \dontrun{
#' ids <- unique(as.character(dL_clean_Ipab$id))
#' metrics <- compute_residual_metrics(
#'   model         = overall_IpaB_pop_6,
#'   dataset       = dL_clean_Ipab,
#'   ids           = ids,
#'   antigen_iso   = "IgG",
#'   scale         = "log",
#'   summary_level = "id_antigen"
#' )
#' }
#'
#' @importFrom dplyr filter mutate inner_join summarise select rename
#' @export
compute_residual_metrics <- function(model,
                                     dataset,
                                     ids,
                                     antigen_iso,
                                     scale         = c("original", "log"),
                                     summary_level = c("id_antigen",
                                                       "pointwise",
                                                       "antigen",
                                                       "overall")) {

  scale         <- match.arg(scale)
  summary_level <- match.arg(summary_level)

  time_var  <- .get_timeindays_var(dataset)
  value_var <- .get_values_var(dataset)

  observed_data <- dataset |>
    dplyr::rename(t = !!time_var, obs = !!value_var) |>
    dplyr::select("id", "t", "obs", "antigen_iso") |>
    dplyr::mutate(id = as.character(.data$id)) |>
    dplyr::filter(
      .data$id          %in% ids,
      .data$antigen_iso == antigen_iso
    )

  if (nrow(observed_data) == 0) {
    stop(
      "No observed data for the specified IDs and antigen_iso.\n",
      "  IDs: ", paste(head(ids, 3), collapse = ", "), "\n",
      "  antigen_iso: ", antigen_iso
    )
  }

  obs_times <- sort(unique(observed_data$t))

  predictions_all <- predict_posterior_at_times(
    model       = model,
    ids         = ids,
    antigen_iso = antigen_iso,
    times       = obs_times
  )

  pred_summary <- predictions_all |>
    dplyr::summarise(
      .by        = c("id", "t"),
      pred_med   = stats::median(.data$res, na.rm = TRUE),
      pred_lower = stats::quantile(.data$res, probs = 0.025, na.rm = TRUE),
      pred_upper = stats::quantile(.data$res, probs = 0.975, na.rm = TRUE)
    ) |>
    dplyr::mutate(id = as.character(.data$id))

  residual_data <- observed_data |>
    dplyr::inner_join(pred_summary, by = c("id", "t"))

  if (scale == "log") {
    n_nonpos_obs  <- sum(residual_data$obs      <= 0)
    n_nonpos_pred <- sum(residual_data$pred_med <= 0)
    if (n_nonpos_obs > 0 || n_nonpos_pred > 0) {
      warning(
        "Removing ", n_nonpos_obs + n_nonpos_pred,
        " non-positive values for log-scale residuals."
      )
      residual_data <- residual_data |>
        dplyr::filter(.data$obs > 0, .data$pred_med > 0)
    }
    residual_data <- residual_data |>
      dplyr::mutate(
        obs_log      = log(.data$obs),
        pred_med_log = log(.data$pred_med),
        residual     = .data$obs_log - .data$pred_med_log,
        abs_residual = abs(.data$residual),
        sq_residual  = .data$residual^2
      ) |>
      dplyr::select(-"obs_log", -"pred_med_log")
  } else {
    residual_data <- residual_data |>
      dplyr::mutate(
        residual     = .data$obs - .data$pred_med,
        abs_residual = abs(.data$residual),
        sq_residual  = .data$residual^2
      )
  }

  if (summary_level == "pointwise") {
    return(residual_data |>
      dplyr::select(
        "id", "antigen_iso", "t", "obs", "pred_med",
        "pred_lower", "pred_upper",
        "residual", "abs_residual", "sq_residual"
      ))
  }

  if (summary_level == "id_antigen") {
    residual_data |>
      dplyr::summarise(
        .by    = c("id", "antigen_iso"),
        MAE    = mean(.data$abs_residual, na.rm = TRUE),
        RMSE   = sqrt(mean(.data$sq_residual, na.rm = TRUE)),
        SSE    = sum(.data$sq_residual, na.rm = TRUE),
        n_obs  = dplyr::n()
      )
  } else if (summary_level == "antigen") {
    residual_data |>
      dplyr::summarise(
        .by   = "antigen_iso",
        MAE   = mean(.data$abs_residual, na.rm = TRUE),
        RMSE  = sqrt(mean(.data$sq_residual, na.rm = TRUE)),
        SSE   = sum(.data$sq_residual, na.rm = TRUE),
        n_obs = dplyr::n()
      )
  } else {
    residual_data |>
      dplyr::summarise(
        MAE   = mean(.data$abs_residual, na.rm = TRUE),
        RMSE  = sqrt(mean(.data$sq_residual, na.rm = TRUE)),
        SSE   = sum(.data$sq_residual, na.rm = TRUE),
        n_obs = dplyr::n()
      )
  }
}


#' Compute Per-Individual MAE for One Model-Antigen-Isotype Combination
#'
#' A thin wrapper around [compute_residual_metrics()] that automatically
#' extracts all IDs from the dataset and returns a tidy tibble.
#'
#' @param model A fitted model object.
#' @param dataset A case_data object.
#' @param antigen_label Character scalar, e.g. `"IpaB"`.
#' @param iso Character scalar isotype, e.g. `"IgG"`.
#' @param scale Character scalar. Default `"log"`.
#'
#' @return A tibble with columns `sid`, `antigen`, `Iso_type`, `mae`.
#'
#' @examples
#' \dontrun{
#' get_mae(overall_IpaB_pop_6, dL_clean_Ipab, "IpaB", "IgG")
#' }
#'
#' @importFrom dplyr mutate select
#' @importFrom tibble tibble
#' @export
get_mae <- function(model, dataset, antigen_label, iso, scale = "log") {

  all_ids <- unique(as.character(dataset$id))

  tryCatch({
    res <- compute_residual_metrics(
      model         = model,
      dataset       = dataset,
      ids           = all_ids,
      antigen_iso   = iso,
      scale         = scale,
      summary_level = "id_antigen"
    )

    res |>
      dplyr::mutate(
        sid      = as.character(.data$id),
        antigen  = antigen_label,
        Iso_type = iso,
        mae      = .data$MAE
      ) |>
      dplyr::select("sid", "antigen", "Iso_type", "mae")

  }, error = function(e) {
    message(paste("Skipped:", antigen_label, iso, "\u2014", e$message))
    tibble::tibble(
      sid      = character(),
      antigen  = character(),
      Iso_type = character(),
      mae      = numeric()
    )
  })
}
