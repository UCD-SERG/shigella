# ============================================================================
# phase0_interactive_reproducibility.R  —  n = 5 pilot
#
# Execution:
#   salloc --time=04:00:00 --cpus-per-task=2 --mem=10G
#   Rscript scripts/phase0_interactive_reproducibility.R
#   exit
#
# Purpose: reproduce a Phase 1 (sbatch) fit under interactive SLURM to confirm
# determinism across allocation modes. Run the n=5 version first as a fast
# smoke test before committing to the n=48 full-cohort run.
# ============================================================================
local({
  flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(flag)) {
    root <- dirname(dirname(normalizePath(sub("^--file=", "", flag))))
    if (file.exists(file.path(root, "DESCRIPTION"))) setwd(root)
  }
})
suppressPackageStartupMessages(library(shigella))

shigella:::run_phase0_diagnostic(
  n             = 5,
  iter_warmup   = 500,
  iter_sampling = 500,
  tag           = "n5"
)
