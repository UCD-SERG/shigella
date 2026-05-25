# Helper: run prep + sample + postprocess for one stratum.
# Accepts the full dataset plus strat/stratum identifiers and slices
# internally, so the caller loop body only needs one function call.
# Returns a list with sr_tibble, cov_summaries, stan_fit, and priors.
#' @keywords internal
#' @noRd
.run_single_stratum <- function(stratum, data, strat,
                                mod, model, chains, parallel_chains,
                                iter_warmup, iter_sampling, seed,
                                adapt_delta, max_treedepth, init,
                                refresh, show_messages, ...) {
  dl_sub <- if (is.na(strat)) {
    data
  } else {
    sub <- data[data[[strat]] == stratum, , drop = FALSE]
    # Restore all non-structural attributes (class, id_var, biomarker_var,
    # time_in_days, value_var, etc.) dropped by [.data.frame subsetting.
    standard_attrs <- c("names", "class", "row.names")
    for (a in setdiff(names(attributes(data)), standard_attrs)) {
      attr(sub, a) <- attr(data, a)
    }
    class(sub) <- class(data)
    sub
  }

  prepped   <- serodynamics::prep_data(dl_sub, add_newperson = FALSE)
  stan_data <- prep_data_stan(prepped)
  priors    <- prep_priors_stan(model = model, ...)
  full_data <- c(stan_data, priors)

  cli::cli_inform(c("i" = "Sampling {.strong {model}} with {chains} chains..."))
  fit <- mod$sample(
    data            = full_data,
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    seed            = seed,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    init            = init,
    refresh         = refresh,
    show_messages   = show_messages
  )

  processed <- postprocess_stan_output(
    stan_fit       = fit,
    ids            = attr(stan_data, "ids"),
    antigens       = attr(stan_data, "antigens"),
    model          = model,
    stratification = stratum
  )

  list(
    sr_tibble     = processed$sr_tibble,
    cov_summaries = processed$cov_summaries,
    stan_fit      = fit,
    priors        = priors
  )
}
