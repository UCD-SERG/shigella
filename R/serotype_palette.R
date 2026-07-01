#' Infecting-serotype colour palette (manuscript v2)
#' @export
serotype_palette <- function() {
  c(
    "Sf2a"   = "#D81B60", # magenta-red (homologous, dominant)
    "Sonnei" = "#1E88E5", # blue (homologous, second dominant)
    "Sf3a"   = "#E65100", # deep orange
    "Sf6"    = "#2E7D32", # forest green
    "Other"  = "grey45"
  )
}

#' @rdname serotype_palette
#' @export
serotype_line_alpha <- function() {
  c(
    "Sf2a" = 0.45, "Sonnei" = 0.40, "Sf3a" = 0.50, "Sf6" = 0.45, "Other" = 0.20
  )
}

#' @rdname serotype_palette
#' @export
serotype_linewidth <- function() {
  c(
    "Sf2a" = 0.5, "Sonnei" = 0.4, "Sf3a" = 0.45, "Sf6" = 0.40, "Other" = 0.2
  )
}
