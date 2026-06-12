#' Convenience driver: locate sensitivity files and build the S3 table
#'
#' Mirrors the driver in the supplement .qmd: lists `sensitivity_*.rda` in
#' `dir`, keeps the Sonnei/Sf3a x IgG/IgA x prior x model files, and renders the
#' S3 table.
#'
#' @param dir Directory of the `.rda` files (`manuscript_data_dir`).
#' @param datasets Named list of case-data objects (see
#'   [get_dataset_for_sensitivity()]).
#' @param scale `"log"` (default).
#' @return A `gt_tbl`.
#' @export
build_table_s3 <- function(dir, datasets, scale = "log") {
  files <- list.files(
    path = fs::path_expand(dir),
    pattern = "^sensitivity_.*\\.rda$", full.names = TRUE
  )
  files <- files[stringr::str_detect(
    basename(files),
    "sensitivity_(Sonnei|Sf3a)_(IgG|IgA)_(primary|diffuse|informative)_(overall|serotype)\\.rda$" # nolint: line_length_linter.
  )]
  res <- build_sensitivity_results(files, datasets = datasets, scale = scale)
  table_s3_sensitivity(res)
}
