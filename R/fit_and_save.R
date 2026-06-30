#' @title Fit a population model and save the result
#' @description Runs `serodynamics::run_mod_pop()` with the given data and MCMC
#'   settings, then saves the fitted object to an `.rda` file. The object inside
#'   the file is named `object_name` so that downstream `load()` calls find the
#'   expected variable name.
#' @param data Data frame of longitudinal immunological observations.
#' @param name Character. File stem; output is written as `<name>.rda` in `dir`.
#' @param object_name Character. Variable name stored inside the `.rda`.
#'   Defaults to `name`. Script 05 passes `"fit_obj"` to match the S3-table
#'   loader contract in `R/tables_supp.R`.
#' @param settings Named list of MCMC settings with elements `nchain`, `nadapt`,
#'   `nburn`, `nmc`, `niter`. Required; no default to avoid dependency on 
#'   globals.
#' @param priors Named list of prior hyperparameters with elements
#'   `mu_hyp_param`, `prec_hyp_param`, `omega_param`, `wishdf_param`,
#'   `prec_logy_hyp_param`. Required; no default to avoid dependency on globals.
#' @param dir Character. Output directory for the `.rda` file. Required; no
#'   default to avoid dependency on globals.
#' @param with_post Logical. Whether to compute posterior summaries. Default
#'   `TRUE`.
#' @return The fitted model object, invisibly.
#' @export
fit_and_save <- function(data,
                         name,
                         object_name = name,
                         settings,
                         priors,
                         dir,
                         with_post = TRUE) {
  obj <- serodynamics::run_mod_pop(
    data = data,
    file_mod = serodynamics::serodynamics_example("model.jags"),
    nchain = settings$nchain,
    nadapt = settings$nadapt,
    nburn = settings$nburn,
    nmc = settings$nmc,
    niter = settings$niter,
    mu_hyp_param = priors$mu_hyp_param,
    prec_hyp_param = priors$prec_hyp_param,
    omega_param = priors$omega_param,
    wishdf_param = priors$wishdf_param,
    prec_logy_hyp_param = priors$prec_logy_hyp_param,
    with_post = with_post
  )

  assign(object_name, obj)
  save(
    list = object_name,
    file = file.path(dir, paste0(name, ".rda")),
    compress = "xz",
    envir = environment()
  )
  cli::cli_inform("saved: {.file {file.path(dir, paste0(name, '.rda'))}}")
  invisible(obj)
}
