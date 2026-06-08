# -----------------------------------------------------------------------------
# Internal: re-attach serocalculator/serodynamics case-data attributes + class
# after a dplyr verb (dplyr drops custom attributes on filter/mutate).
# Copies every attribute from `template` except the base data-frame ones.
# -----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
.reattach_case_attrs <- function(new, template) {
  src <- attributes(template)
  drop <- c("names", "row.names", "class")
  for (nm in setdiff(names(src), drop)) {
    attr(new, nm) <- src[[nm]]
  }
  class(new) <- class(template)
  new
}
