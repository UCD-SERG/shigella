#!/bin/bash
# ==========================================================================
# diagnose_mini_array.sbatch — 6-task array test
#
# Goal: After single srun works, test if SMALL array runs OK. This catches
# the most likely culprit: concurrent-write conflict on /tmp/<user>/cmdstan_bin
# when multiple tasks try to compile/cache the same binary simultaneously.
#
# This uses the SAME 02_run_array.R logic but with R_TOTAL=2 (so 6 tasks total
# = 3 scenarios x 2 reps each), 3 concurrent.
#
# Submit:
#   sbatch slurm/diagnose_mini_array.sbatch
#
# Expected: all 6 tasks complete in ~30 min
# ==========================================================================

#SBATCH --job-name=diag_mini
#SBATCH --output=logs/diag_mini_%A_%a.out
#SBATCH --error=logs/diag_mini_%A_%a.err

# 2 reps x 3 scenarios = 6 tasks, 3 concurrent
#SBATCH --array=1-6%3
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate r_chapter2

echo "=== Job: $SLURM_JOB_NAME  Task: $SLURM_ARRAY_TASK_ID ==="
echo "=== Host: $(hostname)  Date: $(date) ==="

# CRITICAL: per-task subdir to avoid /tmp collisions
# Each task writes its compile cache to a UNIQUE directory.
# This prevents the "everyone writing to same /tmp" hang.
export STAN_COMPILE_DIR="/tmp/$USER/cmdstan_bin_task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$STAN_COMPILE_DIR"
echo "=== Stan compile dir: $STAN_COMPILE_DIR ==="

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# CRITICAL: R_TOTAL=2 so 02_run_array.R sees small per-scenario count
export R_TOTAL=2

cd "$HOME/chapter2"

mkdir -p logs outputs/02_array_test

# Override OUTPUT_DIR if 02_run_array.R supports it; otherwise use main array
# We'll use a separate output folder to keep test results separate.
export OUTPUT_DIR="outputs/02_array_test"

Rscript scripts/02_run_array_v2.R

echo "=== Finished at $(date) ==="
