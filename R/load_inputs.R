#' @title Load input data objects from `.rda` files
#' @description Loads a set of named objects from `.rda` files into the caller's
#'   environment. Each `<stem>.rda` must exist in `dir` and contain an object
#'   named `<stem>`, as produced by `data-raw/00_build_case_data.R`.
#' @param stems Character vector of file stems.
#' @param dir Character. Directory containing the `.rda` files. Required; no
#'   default to avoid dependency on globals.
#' @param envir Environment into which objects are loaded. Defaults to the
#'   caller's environment.
#' @return `stems`, invisibly.
#' @export
load_inputs <- function(stems, dir, envir = parent.frame()) {
  for (s in stems) {
    f <- file.path(dir, paste0(s, ".rda"))
    if (!file.exists(f)) {
      cli::cli_abort(c("Missing input: {.file {f}}",
        "i" = "Run {.file data-raw/00_build_case_data.R} first."
      ))
    }
    load(f, envir = envir)
  }
  invisible(stems)
}
