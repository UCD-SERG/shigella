#' Per-individual MAE for one model x antigen x isotype (safe wrapper)
#'
#' Thin wrapper around [compute_residual_metrics()] returning a uniform
#' `sid / antigen / Iso_type / mae` tibble, with errors downgraded to an empty
#' tibble (so a missing antigen-isotype does not abort a whole comparison).
#'
#' @param model Fitted `sr_model`.
#' @param dataset Case data the model was fit to.
#' @param antigen_label Antigen label for the output.
#' @param iso Isotype (`"IgG"`/`"IgA"`).
#' @param scale `"log"` (default) or `"original"`.
#' @return Tibble `sid, antigen, Iso_type, mae`.
#' @export
get_mae <- function(model, dataset, antigen_label, iso, scale = "log") {
  all_ids <- unique(as.character(dataset$id))

  tryCatch({
    compute_residual_metrics(
      model = model, dataset = dataset, ids = all_ids,
      antigen_iso = iso, scale = scale, summary_level = "id_antigen"
    ) |>
      dplyr::mutate(sid = as.character(.data$id),
                    antigen = antigen_label, Iso_type = iso, mae = .data$MAE) |>
      dplyr::select("sid", "antigen", "Iso_type", "mae")
  }, error = function(e) {
    message("Skipped: ", antigen_label, " ", iso, " - ", e$message)
    tibble::tibble(sid = character(), antigen = character(),
                   Iso_type = character(), mae = numeric())
  })
}
