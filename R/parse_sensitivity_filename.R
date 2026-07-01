#' Parse a `sensitivity_*.rda` file name into its components
#'
#' Expects `sensitivity_<Antigen>_<Iso>_<prior>_<model>.rda`.
#'
#' @param path File path.
#' @return One-row tibble: `file, antigen, iso, prior_key, model`.
#' @export
parse_sensitivity_filename <- function(path) {
  nm <- basename(path)
  pattern <- "^sensitivity_(.+)_(IgG|IgA)_(primary|diffuse|informative)_(overall|serotype)\\.rda$" # nolint: line_length_linter.
  m <- regmatches(nm, regexec(pattern, nm))[[1]]
  if (length(m) != 5) {
    cli::cli_abort("Unexpected sensitivity file name format: {.file {nm}}")
  }
  tibble::tibble(
    file = path, antigen = m[2], iso = m[3],
    prior_key = m[4], model = m[5]
  )
}
