# ==========================================================================
# 00_compile_test.R — Shiva-compatible
#
# Goal: Verify all Stan models compile correctly on Shiva.
# Run this FIRST before any actual fitting.
#
# KEY DIFFERENCE from Mercury version:
#   - Compiled binaries are written to /tmp/<user>/cmdstan_bin/ instead of
#     next to the .stan source file. This bypasses /home noexec restrictions
#     that cause "system error 13, Permission denied" on Shiva.
#
# ==========================================================================

# Set working directory (use Shiva absolute path)
setwd("~/chapter2")

library(cmdstanr)

# ==========================================================================
# 1. Configure compile output directory
# ==========================================================================
# We write compiled binaries to /tmp because /home may be mounted noexec on
# HPC systems. /tmp is always writable + executable for the current user.

user <- Sys.getenv("USER", unset = "default")
compile_dir <- Sys.getenv("STAN_COMPILE_DIR",
                          unset = file.path("/tmp", user, "cmdstan_bin"))

if (!dir.exists(compile_dir)) {
  dir.create(compile_dir, recursive = TRUE, mode = "0755")
}

cat("=== cmdstan setup ===\n")
cat("cmdstanr version:", as.character(packageVersion("cmdstanr")), "\n")
cat("CmdStan path:", cmdstan_path(), "\n")
cat("CmdStan version:", cmdstan_version(), "\n")
cat("Compile output dir:", compile_dir, "\n\n")

# ==========================================================================
# 2. CRITICAL: Clean stale binaries from inst/stan/
# ==========================================================================

stan_dir <- "inst/stan"
stale_files <- list.files(stan_dir, full.names = TRUE)
# Keep only files ending in .stan
binaries_to_remove <- stale_files[!grepl("\\.stan$", stale_files)]
binaries_to_remove <- binaries_to_remove[!grepl("\\.hpp$", binaries_to_remove)]

if (length(binaries_to_remove) > 0) {
  cat("=== Removing stale binaries from", stan_dir, "===\n")
  for (f in binaries_to_remove) {
    cat("  rm:", f, "\n")
    file.remove(f)
  }
  cat("\n")
} else {
  cat("=== inst/stan/ is clean (no stale binaries) ===\n\n")
}

# ==========================================================================
# 3. Compile each Stan model
# ==========================================================================
stan_files <- c(
  "model_1"           = "inst/stan/model_1.stan",
  "model_2"           = "inst/stan/model_2.stan",
  "model_1_time_est"  = "inst/stan/model_1_time_est.stan",
  "model_2_time_est"  = "inst/stan/model_2_time_est.stan"
)

compile_results <- list()

for (model_name in names(stan_files)) {
  stan_path <- stan_files[[model_name]]

  if (!file.exists(stan_path)) {
    cat(sprintf("[SKIP] %s — file not found: %s\n", model_name, stan_path))
    compile_results[[model_name]] <- "MISSING"
    next
  }

  cat(sprintf("\n=== Compiling %s ===\n", model_name))
  t0 <- Sys.time()

  result <- tryCatch({
    mod <- cmdstan_model(
      stan_file = stan_path,
      dir       = compile_dir,    # KEY: write binary to /tmp, not /home
      compile   = TRUE
    )
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf("[OK] %s compiled in %.1f seconds\n", model_name, elapsed))
    cat(sprintf("     Binary: %s\n", mod$exe_file()))
    "OK"
  }, error = function(e) {
    cat(sprintf("[ERROR] %s failed:\n", model_name))
    cat(conditionMessage(e), "\n")
    "FAILED"
  })

  compile_results[[model_name]] <- result
}

# ==========================================================================
# 4. Summary
# ==========================================================================
cat("\n=== Compilation Summary ===\n")
for (name in names(compile_results)) {
  status <- compile_results[[name]]
  symbol <- switch(status,
    "OK"      = "OK     ",
    "FAILED"  = "FAIL   ",
    "MISSING" = "MISSING")
  cat(sprintf("  [%s] %s\n", symbol, name))
}

n_failed  <- sum(unlist(compile_results) == "FAILED")
n_missing <- sum(unlist(compile_results) == "MISSING")

if (n_failed > 0) {
  stop(sprintf("\n%d models failed to compile. Fix errors before proceeding.\n",
               n_failed))
}
if (n_missing > 0) {
  warning(sprintf("\n%d Stan files are missing. Place them in inst/stan/\n",
                  n_missing))
}

cat("\nAll available models compile successfully.\n")
cat("Compiled binaries cached in:", compile_dir, "\n")
cat("Ready to proceed to sanity_check.R or 02_run_scenarios.R.\n")
