## Example: run_phase1_diagnostic()
##
## Runs the Phase 1 SLURM single-job reproducibility diagnostic.
## Intended for use inside a SLURM sbatch job. Requires a compiled
## cmdstan installation and matching Phase 0 outputs for comparison.

if (interactive()) {

  if (requireNamespace("cmdstanr", quietly = TRUE)) {

    result <- run_phase1_diagnostic(
      n             = 5,
      iter_warmup   = 200,
      iter_sampling = 200,
      tag           = "n5",
      output_dir    = file.path(tempdir(), "phase1"),
      phase0_dir    = file.path(tempdir(), "phase0"),
      chains        = 1L
    )

    names(result)  # status, elapsed_min, diagnostic_summary, omega_B_summary, ...
  }
}
