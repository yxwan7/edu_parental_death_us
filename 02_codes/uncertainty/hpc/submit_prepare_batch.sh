#!/bin/bash
# ==============================================================================
# submit_prepare_batch.sh
#
# Create or verify the batch settings before draw jobs are released.
# ==============================================================================

#SBATCH --job-name=edu_prepare
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00

set -euo pipefail

: "${EDU_PARENT_LOSS_REPO:?set EDU_PARENT_LOSS_REPO to the project directory}"
: "${EDU_PARENT_LOSS_SHARD_ROOT:?set EDU_PARENT_LOSS_SHARD_ROOT to the batch root}"
: "${EDU_PARENT_LOSS_BATCH_ID:?set EDU_PARENT_LOSS_BATCH_ID}"
: "${EDU_PARENT_LOSS_BATCH_MODE:?set EDU_PARENT_LOSS_BATCH_MODE}"
: "${EDU_PARENT_LOSS_M_EXPECTED:?set EDU_PARENT_LOSS_M_EXPECTED}"
: "${EDU_PARENT_LOSS_CODE_VERSION:?set EDU_PARENT_LOSS_CODE_VERSION}"

module purge
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate "${EDU_PARENT_LOSS_CONDA_ENV:-edu_parental_death}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1

cd "$EDU_PARENT_LOSS_REPO"
Rscript --vanilla 02_codes/uncertainty/hpc/record_run_environment.r
