# Helper: run prep + sample + postprocess for one stratum.
# Returns a list with sr_tibble, cov_summaries, stan_fit, and priors.
#' @keywords internal
#' @noRd
.run_single_stratum <- function(mod, dl_sub, model, chains, parallel_chains,
                                iter_warmup, iter_sampling, seed,
                                adapt_delta, max_treedepth, init,
                                refresh, show_messages, stratum, ...) {
  prepped   <- serodynamics::prep_data(dl_sub)
  stan_data <- shigella::prep_data_stan(prepped)
  priors    <- shigella::prep_priors_stan(model = model, ...)
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

  processed <- shigella::postprocess_stan_output(
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
