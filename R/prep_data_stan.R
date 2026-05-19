#' @title Prepare data for the Stan backend
#' @description
#' Converts case data into the structured list that the Stan models
#' (`model_1.stan`, `model_2.stan`) expect. Accepts either a raw
#' `case_data` object (the typical entry point) or a
#' `prepped_jags_data` object already produced by
#' [serodynamics::prep_data()] (the JAGS-side prep step).
#'
#' When given a `case_data` object, this function internally calls
#' [serodynamics::prep_data()] with `add_newperson = FALSE` (Stan
#' handles posterior prediction in the `generated quantities` block
#' rather than via a dummy missing-data subject).
#'
#' The Stan models expect:
#' - `N`: number of subjects
#' - `K`: number of antigen-isotype biomarkers
#' - `P`: number of kinetic parameters (always 5)
#' - `max_obs`: max number of observations per subject
#' - `n_obs[N]`: actual number of observations per subject
#' - `time_obs[N, max_obs]`: observation times (NA -> 0, ignored by
#'   the likelihood via the `n_obs[i]` guard)
#' - `log_y[N, max_obs, K]`: log-transformed antibody observations
#'
#' @param data either a `case_data` object (output of
#'   [sim_correlated_case_data()] or [serodynamics::as_case_data()])
#'   or a `prepped_jags_data` object (output of
#'   [serodynamics::prep_data()]).
#' @param drop_newperson [logical] whether to drop the JAGS dummy
#'   "newperson" row if it is present in a `prepped_jags_data`
#'   input. Default `TRUE`. Has no effect when `data` is a
#'   `case_data` object because the internal `prep_data()` call uses
#'   `add_newperson = FALSE`.
#'
#' @returns a named [list] with attributes `ids` and `antigens`,
#'   ready to pass to a compiled Stan model via
#'   `cmdstanr::cmdstan_model()$sample()`.
#' @export
#' @example inst/examples/prep_data_stan-examples.R
prep_data_stan <- function(data,
                           drop_newperson = TRUE) {

  # Route prepped_jags_data directly; convert case_data via helper
  if (inherits(data, "prepped_jags_data")) {
    prepped_jags_data <- data
  } else if (inherits(data, "case_data")) {
    prepped_jags_data <- .case_data_to_prepped_jags(data)
  } else {
    cli::cli_abort(c(
      "{.arg data} must be a {.cls case_data} or {.cls prepped_jags_data}
      object", "i" = "Got an object of class {.cls {class(data)}}."
    ))
  }

  # Extract arrays
  smpl_t  <- prepped_jags_data$smpl.t    # [nsubj, max_visits]
  logy    <- prepped_jags_data$logy      # [nsubj, max_visits, K]
  nsmpl   <- as.integer(prepped_jags_data$nsmpl)
  K       <- prepped_jags_data$n_antigen_isos
  N_full  <- prepped_jags_data$nsubj
  ids_all <- attr(prepped_jags_data, "ids")

  # Drop the "newperson" dummy row if present
  if (drop_newperson && "newperson" %in% ids_all) {
    keep_idx <- which(ids_all != "newperson")
    smpl_t   <- smpl_t[keep_idx, , drop = FALSE]
    logy     <- logy[keep_idx, , , drop = FALSE]
    nsmpl    <- nsmpl[keep_idx]
    ids_kept <- ids_all[keep_idx]
    N        <- length(keep_idx)
  } else {
    ids_kept <- ids_all
    N        <- N_full
  }

  max_obs <- ncol(smpl_t)
  P       <- 5L

  # Replace NA with 0; Stan ignores these via the n_obs[i] guard in
  # the likelihood loop (for (t_idx in 1:n_obs[i])).
  time_obs <- smpl_t
  time_obs[is.na(time_obs)] <- 0
  log_y <- logy
  log_y[is.na(log_y)] <- 0

  # Sanity checks
  .validate_stan_arrays(nsmpl, max_obs)

  antigens <- attr(prepped_jags_data, "antigens")
  stan_data <- list(
    N        = N,
    K        = as.integer(K),
    P        = P,
    max_obs  = as.integer(max_obs),
    n_obs    = nsmpl,
    time_obs = time_obs,
    log_y    = log_y
  )

  # Attach metadata for postprocessing
  attr(stan_data, "ids")      <- ids_kept
  attr(stan_data, "antigens") <- antigens
  return(stan_data)
}
