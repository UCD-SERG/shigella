#' Parse a Sensitivity Analysis RDA Filename into Metadata Components
#'
#' Extracts antigen, isotype, prior configuration, and model type from a
#' standardised sensitivity-analysis RDA filename.
#'
#' @param path Character scalar. Full or relative path to one RDA file.
#'
#' @return A one-row tibble with columns `file`, `antigen`, `iso`,
#'   `prior_key`, `model`.
#'
#' @examples
#' parse_sensitivity_filename(
#'   "sensitivity_Sonnei_IgG_primary_overall.rda"
#' )
#'
#' @importFrom tibble tibble
#' @importFrom stringr str_remove
#' @export
parse_sensitivity_filename <- function(path) {
  nm <- basename(path) |>
    stringr::str_remove("^sensitivity_") |>
    stringr::str_remove("\\.rda$")

  parts <- strsplit(nm, "_")[[1]]

  if (length(parts) != 4L) {
    stop(
      "Unexpected file name format (expected 4 underscore-delimited tokens): ",
      basename(path)
    )
  }

  tibble::tibble(
    file      = path,
    antigen   = parts[1],
    iso       = parts[2],
    prior_key = parts[3],
    model     = parts[4]
  )
}


#' Recode Sensitivity Prior Keys to Display Labels
#'
#' @param x Character vector of prior keys.
#' @return A character vector with recoded values.
#'
#' @examples
#' recode_prior_label(c("primary", "diffuse", "informative", "unknown"))
#'
#' @importFrom dplyr case_when
#' @export
recode_prior_label <- function(x) {
  dplyr::case_when(
    x == "primary"     ~ "Primary",
    x == "diffuse"     ~ "Diffuse",
    x == "informative" ~ "Informative",
    TRUE               ~ x
  )
}


#' Retrieve the Case-Data Object for a Sensitivity Analysis
#'
#' Returns the appropriate case_data object for a given antigen and model
#' type by looking it up from the package's exported data.
#'
#' @param antigen Character scalar (e.g., `"Sonnei"`, `"Sf3a"`).
#' @param model Character scalar (`"overall"` or `"serotype"`).
#'
#' @return The relevant case_data object, or `NULL` if not found.
#'
#' @examples
#' \dontrun{
#' get_dataset_for_sensitivity("Sonnei", "overall")
#' }
#'
#' @export
get_dataset_for_sensitivity <- function(antigen, model) {
  antigen_lower <- tolower(antigen)

  obj_name <- if (model == "overall") {
    paste0("dL_clean_", antigen_lower, "_new")
  } else if (model == "serotype") {
    paste0("dL_serotype_", antigen_lower)
  } else {
    return(NULL)
  }

  # Try to get from package data namespace, then from global env
  tryCatch(
    utils::getFromNamespace(obj_name, "shigella"),
    error = function(e1) {
      tryCatch(
        get(obj_name, envir = .GlobalEnv),
        error = function(e2) {
          warning("Dataset not found: ", obj_name)
          NULL
        }
      )
    }
  )
}


#' Load a Fitted Model Object from an RDA File
#'
#' @param path Character scalar. Full path to the RDA file.
#' @return The R object stored as `fit_obj` in the RDA file.
#'
#' @examples
#' \dontrun{
#' model <- load_fit_obj("path/to/sensitivity_Sonnei_IgG_primary_overall.rda")
#' }
#'
#' @export
load_fit_obj <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)

  if (!exists("fit_obj", envir = e, inherits = FALSE)) {
    stop("'fit_obj' not found in: ", basename(path))
  }
  get("fit_obj", envir = e, inherits = FALSE)
}


#' Compute Per-Individual MAE for All Sensitivity Analysis Models
#'
#' Iterates over sensitivity-analysis RDA files, loads each fitted model,
#' retrieves the corresponding case-data, and computes per-individual MAE.
#'
#' @param files Character vector of full paths to sensitivity RDA files.
#' @param scale Character scalar passed to [get_mae()]. Default `"log"`.
#'
#' @return A tibble with columns `biomarker`, `prior`, `model`, `mae`.
#'
#' @examples
#' \dontrun{
#' files <- list.files("path/to/data", pattern = "^sensitivity_.*\\.rda$",
#'                     full.names = TRUE)
#' results <- build_sensitivity_results(files)
#' }
#'
#' @importFrom purrr map_dfr pmap_dfr
#' @importFrom dplyr mutate arrange
#' @importFrom tibble tibble
#' @export
build_sensitivity_results <- function(files, scale = "log") {

  file_info <- purrr::map_dfr(files, parse_sensitivity_filename)

  results <- purrr::pmap_dfr(
    list(
      file_info$file,
      file_info$antigen,
      file_info$iso,
      file_info$prior_key,
      file_info$model
    ),
    function(file, antigen, iso, prior_key, model) {

      dataset <- get_dataset_for_sensitivity(antigen = antigen, model = model)

      if (is.null(dataset)) {
        message("Skipped (no dataset mapping): ", basename(file))
        return(tibble::tibble(
          biomarker = character(), prior = character(),
          model = character(), mae = numeric()
        ))
      }

      fit_obj <- load_fit_obj(file)

      mae_tbl <- get_mae(
        model         = fit_obj,
        dataset       = dataset,
        antigen_label = antigen,
        iso           = iso,
        scale         = scale
      )

      if (nrow(mae_tbl) == 0) {
        message("Skipped (empty MAE): ", basename(file))
        return(tibble::tibble(
          biomarker = character(), prior = character(),
          model = character(), mae = numeric()
        ))
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
      biomarker = factor(
        .data$biomarker,
        levels = c("Sonnei IgG", "Sonnei IgA", "Sf3a IgG", "Sf3a IgA")
      ),
      prior = factor(.data$prior,
                     levels = c("Primary", "Diffuse", "Informative")),
      model = factor(.data$model, levels = c("overall", "serotype"))
    ) |>
    dplyr::arrange(.data$biomarker, .data$prior, .data$model)
}


#' Build a gt Sensitivity Analysis Table (S3 Table)
#'
#' @param sensitivity_results Tibble from [build_sensitivity_results()].
#' @return A `gt_tbl` object.
#'
#' @examples
#' \dontrun{
#' results <- build_sensitivity_results(files)
#' gt_s3 <- build_table_s3_gt(results)
#' }
#'
#' @importFrom dplyr mutate select case_when
#' @importFrom tidyr pivot_wider
#' @importFrom gt gt fmt_number cols_label cols_width tab_header
#'   tab_source_note tab_options md
#' @export
build_table_s3_gt <- function(sensitivity_results) {

  tbl <- sensitivity_results |>
    tidyr::pivot_wider(
      names_from  = "model",
      values_from = "mae",
      names_prefix = "mae_"
    ) |>
    dplyr::mutate(
      delta_mae = .data$mae_serotype - .data$mae_overall,
      conclusion = dplyr::case_when(
        .data$delta_mae < -0.1 ~ "Serotype better",
        .data$delta_mae >  0.1 ~ "Overall better",
        TRUE                   ~ "Similar"
      ),
      biomarker_group = dplyr::case_when(
        .data$biomarker == "Sonnei IgG" ~ "*S. sonnei* IgG",
        .data$biomarker == "Sonnei IgA" ~ "*S. sonnei* IgA",
        .data$biomarker == "Sf3a IgG"   ~ "*S. flexneri* 3a IgG",
        .data$biomarker == "Sf3a IgA"   ~ "*S. flexneri* 3a IgA",
        TRUE                            ~ as.character(.data$biomarker)
      )
    ) |>
    dplyr::select(
      "biomarker_group",
      Biomarker        = "biomarker",
      Prior            = "prior",
      `MAE (Overall)`  = "mae_overall",
      `MAE (Serotype)` = "mae_serotype",
      `ΔMAE`           = "delta_mae",
      Conclusion       = "conclusion"
    )

  tbl |>
    gt::gt(groupname_col = "biomarker_group") |>
    gt::fmt_number(
      columns = c("MAE (Overall)", "MAE (Serotype)", "\u0394MAE"),
      decimals = 2
    ) |>
    gt::cols_label(
      Biomarker        = "Biomarker",
      Prior            = "Prior",
      `MAE (Overall)`  = "MAE (Overall)",
      `MAE (Serotype)` = "MAE (Serotype)",
      `ΔMAE`           = gt::md("&Delta;MAE"),
      Conclusion       = "Conclusion"
    ) |>
    gt::tab_header(
      title = gt::md("**Sensitivity analysis: prior specification robustness**")
    ) |>
    gt::tab_source_note(gt::md(
      "&Delta;MAE = MAE(serotype) &minus; MAE(overall); negative = serotype advantage."
    )) |>
    gt::tab_options(
      table.font.size        = gt::px(10),
      source_notes.font.size = gt::px(10)
    )
}
