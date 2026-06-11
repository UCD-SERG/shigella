#' Build the prior-sensitivity MAE results across all sensitivity fits
#'
#' For each `sensitivity_*.rda` file, load `fit_obj`, look up its case data, and
#' compute the mean per-individual MAE (via [get_mae()]). Used by
#' [table_s3_sensitivity()].
#'
#' @param files Character vector of `sensitivity_*.rda` paths.
#' @param datasets Named list of case-data objects (see
#'   [get_dataset_for_sensitivity()]).
#' @param scale `"log"` (default) or `"original"`, passed to [get_mae()].
#' @return Tibble `biomarker, prior, model, mae` (ordered factors).
#' @export
build_sensitivity_results <- function(files, datasets, scale = "log") {
  file_info <- purrr::map_dfr(files, parse_sensitivity_filename)

  results <- purrr::pmap_dfr(
    list(file_info$file, file_info$antigen, file_info$iso,
         file_info$prior_key, file_info$model),
    function(file, antigen, iso, prior_key, model) {

      dataset <- get_dataset_for_sensitivity(antigen, model, datasets)
      if (is.null(dataset)) {
        cli::cli_inform("Skipped (no dataset mapping): {.file {basename(file)}}")
        return(tibble::tibble(biomarker = character(), prior = character(),
                              model = character(), mae = numeric()))
      }

      fit_obj <- load_fit_obj(file)
      mae_tbl <- get_mae(model = fit_obj, dataset = dataset,
                         antigen_label = antigen, iso = iso, scale = scale)
      if (nrow(mae_tbl) == 0) {
        cli::cli_inform("Skipped (empty MAE): {.file {basename(file)}}")
        return(tibble::tibble(biomarker = character(), prior = character(),
                              model = character(), mae = numeric()))
      }

      tibble::tibble(
        biomarker = paste(antigen, iso),
        prior     = recode_prior_label(prior_key),
        model     = model,
        mae       = mean(mae_tbl$mae, na.rm = TRUE)
      )
    }
  )

  results |>
    dplyr::mutate(
      biomarker = factor(.data$biomarker,
                         levels = c("Sonnei IgG", "Sonnei IgA",
                                    "Sf3a IgG", "Sf3a IgA")),
      prior = factor(.data$prior, levels = c("Primary", "Diffuse", "Informative")),
      model = factor(.data$model, levels = c("overall", "serotype"))
    ) |>
    dplyr::arrange(.data$biomarker, .data$prior, .data$model)
}
