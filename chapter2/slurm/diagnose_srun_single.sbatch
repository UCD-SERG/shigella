#!/bin/bash
# ==========================================================================
# diagnose_srun_single.sbatch — Single SLURM task test
#
# Goal: After diagnose_single_fit.R works on login node, verify that the
# SAME fit works inside a SLURM compute job. This isolates the
# "is SLURM environment OK?" question from "is the code OK?".
#
# Submit:
#   cd ~/chapter2
#   sbatch slurm/diagnose_srun_single.sbatch
#
# Monitor:
#   squeue -u $USER
#   tail -f logs/diag_srun_<jobid>.out
# ==========================================================================

#SBATCH --job-name=diag_single
#SBATCH --output=logs/diag_srun_%j.out
#SBATCH --error=logs/diag_srun_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate r_chapter2

echo "=== Job: $SLURM_JOB_NAME  Job ID: $SLURM_JOB_ID ==="
echo "=== Host: $(hostname)  Date: $(date) ==="
echo "=== R: $(which R) ==="
echo ""

# Per-task local /tmp
export STAN_COMPILE_DIR="/tmp/$USER/cmdstan_bin"
mkdir -p "$STAN_COMPILE_DIR"
echo "=== Stan compile dir: $STAN_COMPILE_DIR ==="
echo "=== Existing binaries: $(ls $STAN_COMPILE_DIR | wc -l) ==="
echo ""

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "$HOME/chapter2"

# Run the same diagnose script on a compute node
Rscript scripts/diagnose_single_fit.R

echo ""
echo "=== Finished at $(date) ==="
