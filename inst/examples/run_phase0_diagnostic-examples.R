## Example: run_phase0_diagnostic()
##
## Runs the Phase 0 interactive SLURM reproducibility diagnostic.
## Requires a compiled cmdstan installation and is intended for use
## under an interactive SLURM allocation (salloc), not CI.

if (interactive()) {

  if (requireNamespace("cmdstanr", quietly = TRUE)) {

    result <- run_phase0_diagnostic(
      n             = 5,
      iter_warmup   = 200,
      iter_sampling = 200,
      tag           = "n5",
      output_dir    = file.path(tempdir(), "phase0"),
      chains        = 1L
    )

    names(result)  # status, elapsed_min, diagnostic_summary, omega_B_summary, ...
  }
}
