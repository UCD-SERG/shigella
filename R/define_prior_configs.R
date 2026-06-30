#' @title Define prior configurations for sensitivity analysis
#' @description Build a named list of primary, diffuse, and informative prior
#'   configurations from a base prior specification. Used by the sensitivity
#'   analysis in `data-raw/05_sensitivity_models.R`.
#' @param base Named list of prior hyperparameters (see `prior_settings` in
#'   `data-raw/_config.R`). Required; no default to avoid dependency on globals.
#' @return Named list with elements `primary`, `diffuse`, and `informative`,
#'   each a modified copy of `base`.
#' @export
define_prior_configs <- function(base) {
  list(
    primary = utils::modifyList(base, list(label = "Primary")),
    diffuse = utils::modifyList(base, list(
      label = "Diffuse",
      prec_hyp_param = base$prec_hyp_param / 4
    )),
    informative = utils::modifyList(base, list(
      label = "Informative",
      prec_hyp_param = base$prec_hyp_param * 4
    ))
  )
}
