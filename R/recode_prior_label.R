#' Recode the prior key to a display label
#' @param x Character vector of prior keys.
#' @return Character vector of labels.
#' @export
recode_prior_label <- function(x) {
  dplyr::case_when(
    x == "primary" ~ "Primary",
    x == "diffuse" ~ "Diffuse",
    x == "informative" ~ "Informative",
    TRUE ~ x
  )
}
