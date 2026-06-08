#' Load the `fit_obj` object from a single `sensitivity_*.rda`
#' @param path File path.
#' @return The fitted `sr_model`.
#' @export
load_fit_obj <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("fit_obj", envir = e, inherits = FALSE)) {
    stop("fit_obj not found in: ", basename(path))
  }
  get("fit_obj", envir = e, inherits = FALSE)
}
