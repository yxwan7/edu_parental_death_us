#!/bin/bash
# ==============================================================================
# submit_smoke_draw.sh
#
# Run one identity draw and one random draw.
# ==============================================================================

#SBATCH --job-name=edu_smoke
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
# Two draws run back to back here, so this job gets the same limit as an array
# task rather than a tighter one.
#SBATCH --time=12:00:00

set -euo pipefail

: "${EDU_PARENT_LOSS_REPO:?set EDU_PARENT_LOSS_REPO to the project directory}"
: "${EDU_PARENT_LOSS_SHARD_ROOT:?set EDU_PARENT_LOSS_SHARD_ROOT to the smoke batch root}"
: "${EDU_PARENT_LOSS_BATCH_ID:?set EDU_PARENT_LOSS_BATCH_ID}"
: "${EDU_PARENT_LOSS_M_EXPECTED:?set EDU_PARENT_LOSS_M_EXPECTED}"
: "${EDU_PARENT_LOSS_CODE_VERSION:?set EDU_PARENT_LOSS_CODE_VERSION}"

module purge
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate "${EDU_PARENT_LOSS_CONDA_ENV:-edu_parental_death}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1
export EDU_PARENT_LOSS_BATCH_MODE=smoke

cd "$EDU_PARENT_LOSS_REPO"
Rscript --vanilla 02_codes/uncertainty/hpc/record_run_environment.r

DRAW_ROOT=""
cleanup() {
  if [[ -n "${DRAW_ROOT:-}" && -d "$DRAW_ROOT" ]]; then
    rm -rf -- "$DRAW_ROOT"
  fi
}
trap cleanup EXIT

run_one_draw() {
  local draw_mode="$1"
  local draw_id="$2"
  DRAW_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/edu_smoke_${SLURM_JOB_ID}_${draw_mode}_XXXXXX")"

  EDU_PARENT_LOSS_DRAW_ROOT="$DRAW_ROOT" \
  EDU_PARENT_LOSS_DRAW_ID="$draw_id" \
  EDU_PARENT_LOSS_DRAW_MODE="$draw_mode" \
    Rscript --vanilla 02_codes/uncertainty/02_run_draw.r

  cleanup
  DRAW_ROOT=""
}

echo "smoke batch $EDU_PARENT_LOSS_BATCH_ID on $(hostname) at $(date -Is)"
run_one_draw identity identity
run_one_draw random 1

find "$EDU_PARENT_LOSS_SHARD_ROOT" -maxdepth 2 -type f -printf '%P\n' | sort
