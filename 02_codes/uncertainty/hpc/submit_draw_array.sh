#!/bin/bash
# ==============================================================================
# submit_draw_array.sh
#
# Run one random uncertainty draw per Slurm array task.
# ==============================================================================

#SBATCH --job-name=edu_draw
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
# A random draw takes about 46 min and peaks below 32G. The walltime carries
# the case where the node is shared, which can stretch a draw several times over.
#SBATCH --time=12:00:00

set -euo pipefail

: "${EDU_PARENT_LOSS_REPO:?set EDU_PARENT_LOSS_REPO to the project directory}"
: "${EDU_PARENT_LOSS_SHARD_ROOT:?set EDU_PARENT_LOSS_SHARD_ROOT to the batch root}"
: "${EDU_PARENT_LOSS_BATCH_ID:?set EDU_PARENT_LOSS_BATCH_ID}"
: "${EDU_PARENT_LOSS_M_EXPECTED:?set EDU_PARENT_LOSS_M_EXPECTED}"
: "${EDU_PARENT_LOSS_CODE_VERSION:?set EDU_PARENT_LOSS_CODE_VERSION}"
: "${SLURM_ARRAY_TASK_ID:?submit this script as a job array}"

module purge
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate "${EDU_PARENT_LOSS_CONDA_ENV:-edu_parental_death}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1
export EDU_PARENT_LOSS_BATCH_MODE=random
export EDU_PARENT_LOSS_DRAW_ID="$SLURM_ARRAY_TASK_ID"
export EDU_PARENT_LOSS_DRAW_MODE=random

cd "$EDU_PARENT_LOSS_REPO"

test -f "$EDU_PARENT_LOSS_SHARD_ROOT/batch_info.rds"
test -f "$EDU_PARENT_LOSS_SHARD_ROOT/seed_ledger.rds"

DRAW_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/edu_draw_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_XXXXXX")"
cleanup() {
  if [[ -n "${DRAW_ROOT:-}" && -d "$DRAW_ROOT" ]]; then
    rm -rf -- "$DRAW_ROOT"
  fi
}
trap cleanup EXIT
export EDU_PARENT_LOSS_DRAW_ROOT="$DRAW_ROOT"

echo "batch $EDU_PARENT_LOSS_BATCH_ID, draw $SLURM_ARRAY_TASK_ID"
echo "node $(hostname), start $(date -Is), scratch $DRAW_ROOT"

Rscript --vanilla 02_codes/uncertainty/02_run_draw.r
