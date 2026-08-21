#!/bin/bash
# ==============================================================================
# run_batch.sh
#
# Submit, inspect and retrieve uncertainty jobs from an existing HPC project
# directory. Transfer the project files before using this helper.
# ==============================================================================

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# ---- Cluster settings ----

# Set EDU_PARENT_LOSS_SSH_HOST to the cluster login host (an ssh config alias works).
ssh_host="${EDU_PARENT_LOSS_SSH_HOST:-}"
# These are written as $HOME rather than ~ because they end up inside sbatch
# arguments such as --output=$root/logs/draw_%A_%a.out, where the tilde is not
# at the start of a word and so nothing expands it -- the logs then land in a
# directory literally named "~". A remote shell expands $HOME anywhere in the
# string. scp is the one caller that needs the other form; see scp_path.
remote_repo="${EDU_PARENT_LOSS_REMOTE_ROOT:-\$HOME/edu_parental_death}"
remote_run_root="${EDU_PARENT_LOSS_REMOTE_RUN_ROOT:-\$HOME/edu_parental_death_runs}"
conda_env="${EDU_PARENT_LOSS_CONDA_ENV:-edu_parental_death}"
default_m=1000
# sacct reports only today unless given a start time.
sacct_since="${EDU_PARENT_LOSS_SACCT_SINCE:-now-30days}"
ssh_options=(-o BatchMode=yes -o ConnectTimeout=30)

# ---- Command line ----

dry_run=0
args=()
for arg in "$@"; do
  if [[ "$arg" == "--dry-run" ]]; then dry_run=1; else args+=("$arg"); fi
done
verb="${args[0]:-help}"
verb_args=("${args[@]:1}")

if command -v sbatch >/dev/null 2>&1; then on_cluster=1; else on_cluster=0; fi

execute() {
  local shown=() argument
  for argument in "$@"; do
    if [[ "$argument" == *[[:space:]]* ]]; then
      shown+=("\"$argument\"")
    else
      shown+=("$argument")
    fi
  done
  echo "  \$ ${shown[*]}"
  if (( dry_run )); then return 0; fi
  "$@"
}

on_hpc() {
  local shell_line="$1"
  if (( on_cluster )); then
    execute bash -lc "$shell_line"
  else
    execute ssh "${ssh_options[@]}" "$ssh_host" "$shell_line"
  fi
}

batch_id() {
  local value="${verb_args[0]:-}"
  if [[ -z "$value" ]]; then
    echo "This command needs an explicit batch ID." >&2
    exit 1
  fi
  if [[ ! "$value" =~ ^[A-Za-z0-9][A-Za-z0-9_-]{2,63}$ ]]; then
    echo "Batch ID must be 3-64 letters, numbers, underscores or hyphens." >&2
    exit 1
  fi
  echo "$value"
}

code_version() {
  if [[ -n "${EDU_PARENT_LOSS_CODE_VERSION:-}" ]]; then
    echo "$EDU_PARENT_LOSS_CODE_VERSION"
    return
  fi
  local value
  value="$(git -C "$repo_root" rev-parse HEAD 2>/dev/null || true)"
  if [[ -z "$value" ]]; then
    echo "Set EDU_PARENT_LOSS_CODE_VERSION, or run from a Git checkout." >&2
    exit 1
  fi
  echo "$value"
}

run_path() { echo "$remote_run_root/$1"; }

# scp speaks SFTP, which has no shell behind it: it expands a leading ~ but
# never $HOME, the reverse of every other remote call here. An absolute path or
# one that already starts with ~ passes through untouched.
scp_path() { echo "${1/#\$HOME/\~}"; }

with_environment() {
  echo "cd $remote_repo && module purge && source \$HOME/miniforge3/etc/profile.d/conda.sh && conda activate $conda_env && $1"
}

batch_exports() {
  local id="$1" mode="$2" m="$3"
  echo "export EDU_PARENT_LOSS_REPO=$remote_repo EDU_PARENT_LOSS_SHARD_ROOT=$(run_path "$id") EDU_PARENT_LOSS_BATCH_ID=$id EDU_PARENT_LOSS_BATCH_MODE=$mode EDU_PARENT_LOSS_M_EXPECTED=$m EDU_PARENT_LOSS_CODE_VERSION=$(code_version) EDU_PARENT_LOSS_CONDA_ENV=$conda_env"
}

# ---- Submission ----

do_smoke() {
  local id root m command submit
  id="$(batch_id)"
  m="${verb_args[1]:-$default_m}"
  root="$(run_path "$id")"

  submit="sbatch --parsable --job-name=$id --output=$root/logs/smoke_%j.out --error=$root/logs/smoke_%j.err 02_codes/hpc/submit_smoke_draw.sh"

  command="mkdir -p $root/logs $root/accounting"
  command="$command && $(batch_exports "$id" smoke "$m")"
  command="$command && job=\$($submit)"
  command="$command && job=\${job%%;*}"
  command="$command && echo smoke job \$job"
  command="$command && squeue -u \$USER"

  on_hpc "$(with_environment "$command")"
}

do_draws() {
  local id root first last concurrency m array_spec command
  id="$(batch_id)"
  if (( ${#verb_args[@]} < 3 )); then
    echo "draws needs BATCH_ID START END [CONCURRENCY] [M]." >&2
    exit 1
  fi
  first="${verb_args[1]}"
  last="${verb_args[2]}"
  concurrency="${verb_args[3]:-40}"
  m="${verb_args[4]:-$default_m}"
  if (( first < 1 || first > last || last > m )); then
    echo "The draw range must satisfy 1 <= START <= END <= M." >&2
    exit 1
  fi
  if (( concurrency > 40 )); then
    echo "CONCURRENCY cannot exceed 40 concurrent array tasks." >&2
    exit 1
  fi

  root="$(run_path "$id")"
  array_spec="$first-$last%$concurrency"

  command="mkdir -p $root/logs $root/accounting"
  command="$command && $(batch_exports "$id" random "$m")"
  command="$command && prepare_job=\$(sbatch --parsable --job-name=prep_$id --output=$root/logs/prepare_%j.out --error=$root/logs/prepare_%j.err 02_codes/hpc/submit_prepare_batch.sh)"
  command="$command && prepare_job=\${prepare_job%%;*}"
  command="$command && array_job=\$(sbatch --parsable --dependency=afterok:\$prepare_job --job-name=draw_$id --array=$array_spec --output=$root/logs/draw_%A_%a.out --error=$root/logs/draw_%A_%a.err 02_codes/hpc/submit_draw_array.sh)"
  command="$command && array_job=\${array_job%%;*}"
  command="$command && echo prepare job \$prepare_job, array job \$array_job"
  command="$command && squeue -u \$USER"

  on_hpc "$(with_environment "$command")"
}

# ---- Inspection and retrieval ----

do_status() {
  local id root command
  id="$(batch_id)"
  root="$(run_path "$id")"

  # The three job names carry the batch ID, so sacct finds this batch's jobs
  # without a submission log to look them up in.
  command="squeue -u \$USER"
  command="$command; sacct --name=$id,prep_$id,draw_$id -S $sacct_since --format=JobID%24,JobName%20,State,ExitCode,Elapsed,MaxRSS,ReqMem,NodeList | tee $root/accounting/sacct_latest.txt"
  command="$command; if [ -d $root ]; then printf 'shards '; find $root/shards -maxdepth 1 -name 'draw_*.rds' | wc -l; printf 'failures '; find $root/failures -maxdepth 1 -name '*.fail' | wc -l; fi"

  on_hpc "$command"
}

do_fetch() {
  local id target
  id="$(batch_id)"
  if (( ${#verb_args[@]} < 2 )); then
    echo "fetch needs BATCH_ID and a local destination." >&2
    exit 1
  fi
  target="${verb_args[1]}"
  mkdir -p "$(dirname "$target")"
  execute scp "${ssh_options[@]}" -r "$ssh_host:$(scp_path "$(run_path "$id")")" "$target"
  echo "fetched batch $id to $target"
}

usage() {
  cat <<'USAGE'
usage: 02_codes/hpc/run_batch.sh <command> [args] [--dry-run]

Sync the project directory to HPC before submitting jobs.

  smoke BATCH_ID [M]
  draws BATCH_ID START END [CONCURRENCY] [M]
  status BATCH_ID
  fetch BATCH_ID LOCAL_DESTINATION
USAGE
}

case "$verb" in
  smoke|draws|status|fetch) ;;
  help) usage; exit 0 ;;
  *) usage; exit 1 ;;
esac

if (( on_cluster )); then
  where="cluster"
else
  : "${ssh_host:?set EDU_PARENT_LOSS_SSH_HOST to the cluster login host}"
  where="ssh $ssh_host"
fi
if (( dry_run )); then suffix=" [dry run]"; else suffix=""; fi
echo "$verb ($where)$suffix"

"do_$verb"
