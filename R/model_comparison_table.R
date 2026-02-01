#' Compare serotype-specific vs overall models using residual metrics
#'
#' Computes per-ID residual metrics for two models on the intersection of IDs
#' present in both datasets, then reports absolute and percent differences.
#'
#' @param model_serospecific Posterior draws (long format) for the serotype-specific model.
#' @param data_serospecific Case dataset used for the serotype-specific model.
#' @param model_overall Posterior draws (long format) for the overall model.
#' @param data_overall Case dataset used for the overall model.
#' @param antigen_iso Character scalar antigen/isotype label.
#' @param scale "original" or "log".
#' @param tie_tol Numeric tolerance to declare a tie.
#'
#' @return A tibble with per-ID MAE/RMSE for each model and deltas/winner labels.
#'
#' @export
make_model_comparison_table <- function(model_serospecific, data_serospecific,
                                        model_overall,     data_overall,
                                        antigen_iso,
                                        scale = c("original", "log"),
                                        tie_tol = 1e-8) {
  scale <- match.arg(scale)

  ids_common <- intersect(unique(data_serospecific$id), unique(data_overall$id))

  m_sero <- compute_residual_metrics(
    model = model_serospecific,
    dataset = data_serospecific,
    ids = ids_common,
    antigen_iso = antigen_iso,
    scale = scale,
    summary_level = "id_antigen"
  ) |>
    dplyr::select(.data$id, .data$antigen_iso, .data$MAE, .data$RMSE, .data$n_obs) |>
    dplyr::rename(
      MAE_serospecific  = .data$MAE,
      RMSE_serospecific = .data$RMSE,
      n_obs_serospecific = .data$n_obs
    )

  m_over <- compute_residual_metrics(
    model = model_overall,
    dataset = data_overall,
    ids = ids_common,
    antigen_iso = antigen_iso,
    scale = scale,
    summary_level = "id_antigen"
  ) |>
    dplyr::select(.data$id, .data$antigen_iso, .data$MAE, .data$RMSE, .data$n_obs) |>
    dplyr::rename(
      MAE_overall  = .data$MAE,
      RMSE_overall = .data$RMSE,
      n_obs_overall = .data$n_obs
    )

  dplyr::full_join(m_sero, m_over, by = c("id", "antigen_iso")) |>
    dplyr::mutate(
      delta_MAE  = .data$MAE_overall  - .data$MAE_serospecific,
      delta_RMSE = .data$RMSE_overall - .data$RMSE_serospecific,
      pct_improve_MAE = dplyr::case_when(
        is.na(.data$MAE_overall) ~ NA_real_,
        abs(.data$MAE_overall) <= .Machine$double.eps ~ NA_real_,
        TRUE ~ 100 * .data$delta_MAE / .data$MAE_overall
      ),
      pct_improve_RMSE = dplyr::case_when(
        is.na(.data$RMSE_overall) ~ NA_real_,
        abs(.data$RMSE_overall) <= .Machine$double.eps ~ NA_real_,
        TRUE ~ 100 * .data$delta_RMSE / .data$RMSE_overall
      ),
      best_MAE = dplyr::case_when(
        is.na(.data$MAE_overall) | is.na(.data$MAE_serospecific) ~ NA_character_,
        abs(.data$delta_MAE) <= tie_tol ~ "tie",
        .data$delta_MAE > 0 ~ "serospecific",
        TRUE ~ "overall"
      ),
      best_RMSE = dplyr::case_when(
        is.na(.data$RMSE_overall) | is.na(.data$RMSE_serospecific) ~ NA_character_,
        abs(.data$delta_RMSE) <= tie_tol ~ "tie",
        .data$delta_RMSE > 0 ~ "serospecific",
        TRUE ~ "overall"
      ),
      best_overall = dplyr::case_when(
        .data$best_MAE == "serospecific" & .data$best_RMSE == "serospecific" ~ "serospecific",
        .data$best_MAE == "overall"      & .data$best_RMSE == "overall"      ~ "overall",
        .data$best_MAE == "tie"          & .data$best_RMSE == "tie"          ~ "tie",
        TRUE ~ "mixed"
      )
    ) |>
    dplyr::arrange(.data$id)
}
