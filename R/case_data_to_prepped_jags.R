# Helper: convert a case_data object to prepped_jags_data via serodynamics.
#' @keywords internal
#' @noRd
.case_data_to_prepped_jags <- function(data) {
  if (!requireNamespace("serodynamics", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg serodynamics} is required.",
      "i" = "Install it before using {.fn prep_data_stan}."
    ))
  }
  serodynamics::prep_data(
    data,
    add_newperson = FALSE
  )
}
