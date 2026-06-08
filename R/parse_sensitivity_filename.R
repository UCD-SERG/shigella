#' Parse a `sensitivity_*.rda` file name into its components
#'
#' Expects `sensitivity_<Antigen>_<Iso>_<prior>_<model>.rda`.
#'
#' @param path File path.
#' @return One-row tibble: `file, antigen, iso, prior_key, model`.
#' @export
parse_sensitivity_filename <- function(path) {
  nm <- basename(path) |>
    stringr::str_remove("^sensitivity_") |>
    stringr::str_remove("\\.rda$")

  parts <- strsplit(nm, "_")[[1]]
  if (length(parts) != 4) {
    stop("Unexpected file name format: ", basename(path))
  }
  tibble::tibble(
    file = path, antigen = parts[1], iso = parts[2],
    prior_key = parts[3], model = parts[4]
  )
}
