# ==========================================================================
# phase0_no_slurm_reproducibility.R
# ==========================================================================

# ----- 0. Setup + paths -----
setwd("~/shigella/chapter2")

cat("\n", strrep("=", 70), "\n", sep = "")
cat(" PHASE 0: NO-SLURM REPRODUCIBILITY TEST\n")
cat(" Purpose: Run identical fit OUTSIDE Slurm to isolate environment\n")
cat(" Goal: 'Does this only happen on Slurm?'\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat(sprintf("Started at:     %s\n", format(Sys.time())))
cat(sprintf("Host:           %s\n", Sys.info()[["nodename"]]))
cat(sprintf("R version:      %s\n", R.version.string))
cat(sprintf("Working dir:    %s\n\n", getwd()))

dir.create("logs/phase0",   recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/phase0", recursive = TRUE, showWarnings = FALSE)

# Status tracker — written incrementally so we know where we crashed even
# if the script dies mid-way. Mirrors Ezra's request for crash-resilient logs.
status_file <- "outputs/phase0/PHASE0_STATUS.txt"
unlink(status_file)
write_status(status_file, "INIT", "Phase 0 started")

# ----- 1. Package loading with version capture -----
cat(strrep("-", 70), "\n", sep = "")
cat("STEP 1: Load packages + capture versions\n")
cat(strrep("-", 70), "\n", sep = "")
write_status(status_file,"LOAD_PACKAGES", "loading")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(cmdstanr)
  library(posterior)
  library(tibble)
  library(shigella)
})

pkg_versions <- c(
  R            = R.version.string,
  cmdstanr     = as.character(packageVersion("cmdstanr")),
  posterior    = as.character(packageVersion("posterior")),
  shigella     = as.character(packageVersion("shigella")),
  serodynamics = as.character(packageVersion("serodynamics"))
)
for (n in names(pkg_versions)) {
  cat(sprintf("  %-15s %s\n", n, pkg_versions[n]))
}

# cmdstan version itself
cmdstan_ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) "UNKNOWN")
cat(sprintf("  %-15s %s\n", "cmdstan", cmdstan_ver))

# Save for later comparison with Phase 1 logs
saveRDS(c(pkg_versions, cmdstan = cmdstan_ver),
        "outputs/phase0/env_versions.rds")

write_status(status_file,"LOAD_PACKAGES", "OK")
cat("\n")

# ----- 2. Verify Stan files present -----
cat(strrep("-", 70), "\n", sep = "")
cat("STEP 2: Verify Stan files\n")
cat(strrep("-", 70), "\n", sep = "")
write_status(status_file,"VERIFY_STAN", "checking")

stan_files <- c(
  m1 = system.file("stan", "model_1.stan", package = "shigella"),
  m2 = system.file("stan", "model_2.stan", package = "shigella")
)
for (n in names(stan_files)) {
  f <- stan_files[n]
  if (file.exists(f)) {
    cat(sprintf("  [OK]   %s -> %s\n", n, f))
  } else {
    cat(sprintf("  [MISS] %s — falling back to inst/stan/%s.stan\n", n, n))
    # Fall back to local inst/stan/ if installed package missing the file
    fallback <- file.path("inst/stan", paste0(n, ".stan"))
    if (file.exists(fallback)) {
      stan_files[n] <- normalizePath(fallback)
      cat(sprintf("         using fallback: %s\n", stan_files[n]))
    }
  }
}

write_status(status_file,"VERIFY_STAN", "OK")
cat("\n")

# ----- 3. compile_dir setup -----
cat(strrep("-", 70), "\n", sep = "")
cat("STEP 3: Set up compile directory\n")
cat(strrep("-", 70), "\n", sep = "")
write_status(status_file,"COMPILE_DIR", "setting up")

user <- Sys.getenv("USER", unset = "unknown")
compile_dir <- file.path("/tmp", user, "cmdstan_bin_phase0")
if (!dir.exists(compile_dir)) {
  dir.create(compile_dir, recursive = TRUE, mode = "0755")
}
cat(sprintf("  compile_dir: %s\n", compile_dir))
cat(sprintf("  existing files: %d\n", length(list.files(compile_dir))))


testfile <- file.path(compile_dir, "_test_write")
writeLines("test", testfile)
if (file.exists(testfile)) {
  cat("  [OK] write test passed\n")
  unlink(testfile)
} else {
  stop("compile_dir is not writable — cannot proceed")
}

# Executable test 
shellscript <- file.path(compile_dir, "_test_exec.sh")
writeLines(c("#!/bin/bash", "echo executable"), shellscript)
Sys.chmod(shellscript, "0755")
exec_out <- tryCatch(
  system(shellscript, intern = TRUE),
  warning = function(w) NULL,
  error   = function(e) NULL
)
if (length(exec_out) > 0 && exec_out == "executable") {
  cat("  [OK] exec test passed (no noexec issue)\n")
} else {
  cat("  [WARN] exec test failed — Stan may not be able to run binaries here\n")
}
unlink(shellscript)

write_status(status_file,"COMPILE_DIR", "OK")
cat("\n")

# ----- 4. Simulate small data (n=5, fixed seed) -----
cat(strrep("-", 70), "\n", sep = "")
cat("STEP 4: Simulate small synthetic data\n")
cat(strrep("-", 70), "\n", sep = "")
write_status(status_file,"SIMULATE", "running")

set.seed(20260513)   # matches meeting date, so reproducible

true_rho_B <- 0.6
omega_B_true <- matrix(c(1, true_rho_B, true_rho_B, 1), 2, 2)

sim_data <- sim_correlated_case_data(
  n                 = 5,                            
  omega_B           = omega_B_true,
  antigen_isos      = c("IgG", "IgA"),
  n_obs_per_subject = 5L
)

cat(sprintf("  n_subjects:  %d\n", length(unique(sim_data$id))))
cat(sprintf("  total rows:  %d\n", nrow(sim_data)))
cat(sprintf("  isotypes:    %s\n", paste(unique(sim_data$antigen_iso), collapse = ", ")))
cat(sprintf("  true rho_B:  %.3f\n", true_rho_B))

saveRDS(sim_data, "outputs/phase0/sim_data_n5.rds")
cat("  saved -> outputs/phase0/sim_data_n5.rds\n")

write_status(status_file,"SIMULATE", "OK")
cat("\n")

# ----- 5. Run single fit -----
cat(strrep("-", 70), "\n", sep = "")
cat("STEP 5: Fit model_2 (Kronecker), light settings, single chain\n")
cat(strrep("-", 70), "\n", sep = "")
write_status(status_file,"FIT", "running")

t_start <- Sys.time()

saveRDS(
  list(scenario = "phase0_no_slurm",
       status   = "FIT_STARTED",
       started_at = format(t_start),
       true_rho_B = true_rho_B,
       n = 5),
  "outputs/phase0/one_fit_n5.rds"
)

fit <- tryCatch({
  run_mod_stan(
    data            = sim_data,
    model           = "model_2",
    chains          = 2,
    iter_warmup     = 500,
    iter_sampling   = 500,
    parallel_chains = 2,
    adapt_delta     = 0.95,
    max_treedepth   = 12,
    init            = 0.1,
    with_post       = TRUE,
    stan_dir        = "inst/stan",
    compile_dir     = compile_dir,
    refresh         = 100,
    show_messages   = TRUE
  )
}, error = function(e) {
  cat("\n  [FIT ERROR]:", conditionMessage(e), "\n")
  cat("  [STACK TRACE]:\n")
  print(sys.calls())
  write_status(status_file,"FIT", paste("CRASHED:", conditionMessage(e)))
  saveRDS(
    list(scenario = "phase0_no_slurm",
         status   = "FIT_FAILED",
         error    = conditionMessage(e),
         crashed_at = format(Sys.time()),
         true_rho_B = true_rho_B,
         n = 5),
    "outputs/phase0/one_fit_n5.rds"
  )
  NULL
})

elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
cat(sprintf("\n  Fit elapsed: %.2f min\n", elapsed))

if (is.null(fit)) {
  cat("\n", strrep("!", 70), "\n", sep = "")
  cat(" PHASE 0 RESULT: FIT CRASHED OUTSIDE SLURM\n")
  cat(" This is a HIGH-SIGNAL finding for Ezra.\n")
  cat(" → Confirms problem is NOT Slurm-specific.\n")
  cat(" → Skip Phase 1-3, jump to Phase 4 over-parameterization diagnosis.\n")
  cat(strrep("!", 70), "\n", sep = "")
  write_status(status_file,"DONE", "Phase 0 FAILED — see error above")
  quit(status = 1)
}

write_status(status_file,"FIT", "OK")
cat("\n")

# ----- 6. Extract diagnostics from cmdstanr fit -----
cat(strrep("-", 70), "\n", sep = "")
cat("STEP 6: Extract diagnostics\n")
cat(strrep("-", 70), "\n", sep = "")
write_status(status_file,"DIAG", "extracting")

sf <- attr(fit, "stan_fit")[[1]]

diag <- sf$diagnostic_summary(diagnostics = c("divergences",
                                              "treedepth",
                                              "ebfmi"))
n_total_draws <- sum(diag$num_divergent) + sum(diag$num_max_treedepth)
total_iters   <- 2 * 500   # chains * iter_sampling

cat(sprintf("  divergent transitions:   %d / %d (%.2f%%)\n",
            sum(diag$num_divergent), total_iters,
            100 * sum(diag$num_divergent) / total_iters))
cat(sprintf("  max-treedepth hits:      %d / %d (%.2f%%)\n",
            sum(diag$num_max_treedepth), total_iters,
            100 * sum(diag$num_max_treedepth) / total_iters))
cat(sprintf("  E-BFMI by chain:         %s\n",
            paste(sprintf("%.3f", diag$ebfmi), collapse = ", ")))

# ESS + R-hat for key parameter Omega_B[1,2]
draws_summary <- tryCatch({
  posterior::summarise_draws(
    sf$draws(variables = "Omega_B[1,2]"),
    "median", "mean", "sd",
    ~quantile(.x, c(0.025, 0.975), na.rm = TRUE),
    "ess_bulk", "rhat"
  )
}, error = function(e) NULL)

if (!is.null(draws_summary)) {
  cat("\n  Omega_B[1,2] posterior summary:\n")
  print(draws_summary)
}

write_status(status_file,"DIAG", "OK")
cat("\n")

# ----- 7. Save full diagnostic bundle -----
cat(strrep("-", 70), "\n", sep = "")
cat("STEP 7: Save diagnostic bundle\n")
cat(strrep("-", 70), "\n", sep = "")
write_status(status_file,"SAVE", "writing rds")

result_bundle <- list(
  scenario            = "phase0_no_slurm",
  status              = "OK",
  elapsed_min         = elapsed,
  started_at          = format(t_start),
  completed_at        = format(Sys.time()),
  host                = Sys.info()[["nodename"]],
  pkg_versions        = pkg_versions,
  cmdstan_version     = cmdstan_ver,
  true_rho_B          = true_rho_B,
  n_subjects          = 5,
  fit_settings        = list(chains = 2, warmup = 500, sampling = 500,
                             adapt_delta = 0.95, max_treedepth = 12),
  diagnostic_summary  = diag,
  omega_B_summary     = draws_summary,
  # Pull full posterior of rho_B (small file, ~2000 doubles)
  rho_B_posterior     = as.vector(posterior::as_draws_array(
                          sf$draws("Omega_B[1,2]")))
)

saveRDS(result_bundle, "outputs/phase0/one_fit_n5.rds")
saveRDS(result_bundle$diagnostic_summary, "outputs/phase0/one_fit_n5_diag.rds")

write_status(status_file,"SAVE", "OK")
cat("  saved -> outputs/phase0/one_fit_n5.rds\n")
cat("  saved -> outputs/phase0/one_fit_n5_diag.rds\n\n")

# ----- 8. Final summary block -----
cat(strrep("=", 70), "\n", sep = "")
cat(" PHASE 0 RESULT SUMMARY\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("  Status:              OK\n"))
cat(sprintf("  Elapsed:             %.2f min\n", elapsed))
cat(sprintf("  True rho_B:          %+.3f\n", true_rho_B))
if (!is.null(draws_summary)) {
  cat(sprintf("  Recovered median:    %+.3f  [%.3f, %.3f]\n",
              draws_summary$median,
              draws_summary$`2.5%`,
              draws_summary$`97.5%`))
  cat(sprintf("  ESS_bulk:            %.0f\n", draws_summary$ess_bulk))
  cat(sprintf("  R-hat:               %.3f\n", draws_summary$rhat))
}
cat(sprintf("  Divergent:           %d / %d\n",
            sum(diag$num_divergent), total_iters))
cat(sprintf("  Max-treedepth hits:  %d / %d\n",
            sum(diag$num_max_treedepth), total_iters))
cat(strrep("=", 70), "\n\n", sep = "")

cat(" NEXT STEP:\n")
cat("   1. Inspect outputs/phase0/one_fit_n5.rds + logs/phase0/*.log\n")
cat("   2. If divergent rate ≤ 5% AND R-hat ≤ 1.01:\n")
cat("        → Proceed to Phase 1 (sbatch slurm/phase1_single.sbatch)\n")
cat("   3. If divergent rate > 10% OR R-hat > 1.02:\n")
cat("        → This is a NO-Slurm reproducible failure.\n")
cat("        → Skip Phase 1-3, jump to Phase 4 over-param diagnosis.\n\n")

write_status(status_file,"DONE", "Phase 0 completed successfully")
