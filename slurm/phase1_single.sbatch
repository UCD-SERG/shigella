#!/bin/bash
# ==========================================================================
# phase1_single.sbatch — SLURM single-task diagnostic run
#
# Purpose: Run the SAME n=5 fit as phase0_interactive_reproducibility.R, but
# inside a SLURM compute job. Direct comparison isolates "is the issue
# SLURM or the code itself?" 
#
# ==========================================================================

set -euo pipefail

# ==========================================================================
# 1. Environment activation
# ==========================================================================
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate r_chapter2

# ==========================================================================
# 2. Banner
# ==========================================================================
echo "========================================================================"
echo " PHASE 1 SLURM SINGLE JOB"
echo "========================================================================"
echo " Job ID:        $SLURM_JOB_ID"
echo " Job name:      $SLURM_JOB_NAME"
echo " Submit dir:    $SLURM_SUBMIT_DIR"
echo " Compute node:  $(hostname)"
echo " Started:       $(date)"
echo " R binary:      $(which R)"
echo " Rscript:       $(which Rscript)"
echo " CPUs/task:     $SLURM_CPUS_PER_TASK"
echo " Mem/node:      ${SLURM_MEM_PER_NODE:-unset}"
echo "========================================================================"

# ==========================================================================
# 3. Per-task compile dir — avoids /home noexec + array collision risk
# ==========================================================================
export STAN_COMPILE_DIR="/tmp/$USER/cmdstan_bin_phase1_${SLURM_JOB_ID}"
mkdir -p "$STAN_COMPILE_DIR"
echo "=== STAN_COMPILE_DIR=$STAN_COMPILE_DIR ==="

# ==========================================================================
# 4. Disable thread oversubscription
# ==========================================================================
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ==========================================================================
# 5. Change to submit dir and run
# ==========================================================================
cd "$SLURM_SUBMIT_DIR"   

mkdir -p logs/phase1 outputs/phase1

Rscript scripts/phase1_single_diagnostic.R

EXIT_CODE=$?

# ==========================================================================
# 6. Footer with exit code
# ==========================================================================
echo ""
echo "========================================================================"
echo " PHASE 1 FINISHED"
echo "========================================================================"
echo " Exit code:    $EXIT_CODE"
echo " Finished:     $(date)"
echo " STAN_COMPILE_DIR final size: $(du -sh $STAN_COMPILE_DIR 2>/dev/null || echo 'gone')"
echo "========================================================================"

# Cleanup compile dir only if success — keep on failure for inspection
if [ "$EXIT_CODE" -eq 0 ]; then
    rm -rf "$STAN_COMPILE_DIR"
fi

exit $EXIT_CODE
