# Helper: locate the Stan source file for the given model.
#' @keywords internal
#' @noRd
.locate_stan_file <- function(model, stan_dir) {
  stan_basename <- paste0(model, ".stan")

  if (is.null(stan_dir)) {
    stan_file <- system.file(
      "stan",
      stan_basename,
      package = "shigella",
      mustWork = FALSE
    )

    # Fallback for interactive development before the package is installed.
    if (identical(stan_file, "") || !file.exists(stan_file)) {
      stan_file <- file.path("inst", "stan", stan_basename)
    }
  } else {
    stan_file <- file.path(stan_dir, stan_basename)
  }

  if (!file.exists(stan_file)) {
    cli::cli_abort(c(
      "Cannot locate Stan file: {.file {stan_file}}",
      "i" = "Working directory is: {.path {getwd()}}",
      "i" = "If running interactively, check that
      {.file inst/stan/{stan_basename}} exists.",
      "i" = "If running from an installed package, use system.file('stan',
      '{stan_basename}', package = 'shigella')."
    ))
  }

  stan_file
}
