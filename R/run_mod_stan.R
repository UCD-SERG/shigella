#' @title Run Stan model — cmdstanr version (Shiva-compatible)
#' @description
#' Fits the two-phase antibody kinetics model using **cmdstanr** instead of
#' rstan. This is required on HPC systems like Shiva where:
#'   - rstan can have toolchain conflicts with conda R
#'   - cmdstanr's `dir` argument lets us write compiled binaries to a
#'     writable/executable location (e.g., /tmp), bypassing /home noexec
#'
#' Output: an `sr_model` tibble with the same column schema as `run_mod()`,
#' so all existing plot / summary functions work unchanged. Stan-specific
#' attributes are also attached:
#'   - `Omega_eps`, `Sigma_eps`: residual covariance (all models)
#'   - `Omega_B`, `Sigma_B`:     biomarker covariance (Model 2 only)
#'   - `Omega_P`, `Sigma_P`:     parameter covariance
#'   - `stan_fit`:               raw CmdStanMCMC object (when with_post = TRUE)
#'
#' @param data case_data object (from sim_correlated_case_data() or 
#' as_case_data())
#' @param model character: "model_1", "model_2"
#' @param chains, iter_sampling, iter_warmup, adapt_delta, max_treedepth, seed, 
#' parallel_chains
#'        standard cmdstanr arguments
#' @param strat optional stratification variable (default NA)
#' @param with_post return raw CmdStanMCMC object as attribute (default FALSE)
#' @param stan_dir Optional directory containing `model_*.stan` files. 
#' If `NULL`, the function first looks for Stan files installed with the 
#' package using `system.file("stan", ..., package = "shigella")`, then falls 
#' back to `inst/stan` for interactive development.
#' @param compile_dir directory where cmdstanr writes compiled binaries.
#'                    Default uses STAN_COMPILE_DIR env var, or 
#'                    /tmp/<user>/cmdstan_bin if /home is noexec.
#' @param init initial value strategy. Numeric value scales down random init
#'             (default 0.1 to avoid -inf in multi_normal_cholesky_lpdf)
#' @param ... additional priors passed to prep_priors_stan()
#' @param iter_sampling Number of post-warmup iterations per chain.
#' @param iter_warmup Number of warmup iterations per chain.
#' @param adapt_delta Target average acceptance probability for Stan sampling.
#' @param max_treedepth Maximum tree depth for Stan NUTS sampling.
#' @param seed Random seed passed to Stan.
#' @param parallel_chains Number of chains to run in parallel.
#' @param refresh Stan progress refresh interval.
#' @param show_messages Logical; whether to show CmdStan messages.
#'
#' @returns sr_model tibble
#' @example inst/examples/run_mod_stan-examples.R
#' @export
run_mod_stan <- function(data,
                         model           = c("model_2", "model_1"),
                         chains          = 4,
                         iter_sampling   = 1000,
                         iter_warmup     = 1000,
                         adapt_delta     = 0.95,
                         max_treedepth   = 12,
                         seed            = sample.int(.Machine$integer.max, 1),
                         strat           = NA,
                         parallel_chains = chains,
                         with_post       = FALSE,
                         stan_dir        = NULL,
                         compile_dir     = NULL,
                         init            = 0.1,
                         refresh         = 200,
                         show_messages   = TRUE,
                         ...) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg cmdstanr} is required.",
      "i" = "Install it with: install.packages('cmdstanr', 
      repos = 'https://mc-stan.org/r-packages/')"
    ))
  }
  if (!requireNamespace("serodynamics", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg serodynamics} is required for prep_data().")
  }

  model <- match.arg(model)

  # ---- Locate Stan source file ----
  stan_basename <- paste0(model, ".stan")
  
  if (is.null(stan_dir)) {
    stan_file <- system.file(
      "stan",
      stan_basename,
      package = "shigella",
      mustWork = FALSE
    )
    
    # Fallback for interactive development before the package is installed.
    if (identical(stan_file, "") || !file.exists(stan_file)) {
      stan_file <- file.path("inst", "stan", stan_basename)
    }
  } else {
    stan_file <- file.path(stan_dir, stan_basename)
  }
  
  if (!file.exists(stan_file)) {
    cli::cli_abort(c(
      "Cannot locate Stan file: {.file {stan_file}}",
      "i" = "Working directory is: {.path {getwd()}}",
      "i" = "If running interactively, check that 
      {.file inst/stan/{stan_basename}} exists.",
      "i" = "If running from an installed package, use system.file('stan', 
      '{stan_basename}', package = 'shigella')."
    ))
  }
  
  cli::cli_inform(c("i" = "Using Stan file: {.file {stan_file}}"))

  # ---- Determine compile output directory ----
  # Priority: argument > environment variable > /tmp fallback
  if (is.null(compile_dir)) {
    compile_dir <- Sys.getenv("STAN_COMPILE_DIR", unset = "")
    if (compile_dir == "") {
      user <- Sys.getenv("USER", unset = "default")
      compile_dir <- file.path("/tmp", user, "cmdstan_bin")
    }
  }
  if (!dir.exists(compile_dir)) {
    dir.create(compile_dir, recursive = TRUE, mode = "0755")
  }
  cli::cli_inform(c("i" = "Compile output directory: {.path {compile_dir}}"))

  # ---- Stratification ----
  if (is.na(strat)) {
    strat_list <- "None"
  } else {
    strat_list <- unique(data[[strat]])
  }

  combined_out <- list()
  stanfit_list <- list()
  cov_list     <- list()

  for (i in strat_list) {
    if (is.na(strat)) {
      dl_sub <- data
    } else {
      dl_sub <- data[data[[strat]] == i, , drop = FALSE]
    }

    # ---- Prep data + priors ----
    prepped   <- serodynamics::prep_data(dl_sub)
    stan_data <- shigella::prep_data_stan(prepped)
    priors <- shigella::prep_priors_stan(model = model, ...)
    full_data <- c(stan_data, priors)

    # ---- Compile model ----
    cli::cli_inform(c("i" = "Compiling {.strong {model}} (or using cache)..."))
    mod <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      dir       = compile_dir,    
      compile   = TRUE
    )

    # ---- Sample ----
    cli::cli_inform(c("i" = "Sampling {.strong {model}} with {chains} 
                      chains..."))
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

    # ---- Postprocess ----
    processed <- shigella::postprocess_stan_output(
      stan_fit       = fit,
      ids            = attr(stan_data, "ids"),
      antigens       = attr(stan_data, "antigens"),
      model          = model,
      stratification = i
    )

    combined_out[[i]] <- processed$sr_tibble
    cov_list[[i]]     <- processed$cov_summaries
    stanfit_list[[i]] <- fit
  }

  sr_out <- dplyr::bind_rows(combined_out)

  sr_out <- sr_out |>
    structure(
      nChain      = chains,
      nParameters = 5L,
      nIterations = iter_sampling + iter_warmup,
      nWarmup     = iter_warmup,
      model_type  = model,
      priors      = attr(priors, "used_stan_priors")
    )

  if (length(cov_list) == 1) {
    for (nm in names(cov_list[[1]])) {
      attr(sr_out, nm) <- cov_list[[1]][[nm]]
    }
  } else {
    attr(sr_out, "cov_by_stratum") <- cov_list
  }

  # Calculate fitted/residuals
  if (requireNamespace("serodynamics", quietly = TRUE)) {
    fit_res <- tryCatch(
      serodynamics:::calc_fit_mod(modeled_dat = sr_out, original_data = dl_sub),
      error = function(e) {
        cli::cli_warn("calc_fit_mod failed: {e$message}")
        NULL
      }
    )
    if (!is.null(fit_res)) {
      attr(sr_out, "fitted_residuals") <- fit_res
    }
  }

  if (with_post) {
    attr(sr_out, "stan_fit") <- stanfit_list
  }

  class(sr_out) <- union("sr_model", class(sr_out))
  return(sr_out)
}
