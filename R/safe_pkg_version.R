# Returns the installed version of a package as a string, or "not installed"
# if the package is absent. Used for robust version-logging in diagnostic
# functions where suggested packages may not be present.
#' @keywords internal
#' @noRd
.safe_pkg_version <- function(pkg) {
  tryCatch(
    as.character(utils::packageVersion(pkg)),
    error = function(e) "not installed"
  )
}
