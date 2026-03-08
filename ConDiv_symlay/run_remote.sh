#!/bin/bash
#SBATCH --job-name=condiv-mem
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=16G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load cmake
  module load openmpi
fi

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export PYTHONUNBUFFERED=1

source "$SCRIPT_DIR/../.venv/bin/activate"
source "$SCRIPT_DIR/../source.sh"

export CONDIV_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT/py:$PROJECT_ROOT/obj:${PYTHONPATH:-}"

PROFILE="${PROFILE:-dimer3}"
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR/test_${PROFILE}}"
RUN_STEPS="${1:-${RUN_STEPS:-20}}"
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"
LAYER_MANIFEST_JSON="${LAYER_MANIFEST_JSON:-$BASE_DIR/layer_manifest.json}"
PROGRESS_LOG_JSONL="${PROGRESS_LOG_JSONL:-$BASE_DIR/training_progress.jsonl}"
STATUS_JSON="${STATUS_JSON:-$BASE_DIR/training_status.json}"
CONDIV_CONVERGENCE_GRAD_NORM="${CONDIV_CONVERGENCE_GRAD_NORM:-1.0}"
CONDIV_CONVERGENCE_UPDATE_NORM="${CONDIV_CONVERGENCE_UPDATE_NORM:-0.05}"
CONDIV_CONVERGENCE_PATIENCE="${CONDIV_CONVERGENCE_PATIENCE:-3}"
CONDIV_AUTO_RESUBMIT="${CONDIV_AUTO_RESUBMIT:-1}"
CONDIV_AUTO_MAX_SUBMISSIONS="${CONDIV_AUTO_MAX_SUBMISSIONS:-0}"
CONDIV_RESUBMIT_COUNT="${CONDIV_RESUBMIT_COUNT:-0}"

export PROFILE BASE_DIR WORKER_LAUNCH LAYER_MANIFEST_JSON
export PROGRESS_LOG_JSONL STATUS_JSON
export CONDIV_CONVERGENCE_GRAD_NORM CONDIV_CONVERGENCE_UPDATE_NORM CONDIV_CONVERGENCE_PATIENCE
export CONDIV_AUTO_RESUBMIT CONDIV_AUTO_MAX_SUBMISSIONS CONDIV_RESUBMIT_COUNT

parse_slurm_task_slots() {
  local raw="${1:-}"
  raw="${raw%%,*}"
  raw="${raw%%(*}"
  printf '%s' "$raw"
}

SLURM_TASK_SLOTS="${SLURM_NTASKS:-}"
if [ -z "$SLURM_TASK_SLOTS" ]; then
  SLURM_TASK_SLOTS="$(parse_slurm_task_slots "${SLURM_NTASKS_PER_NODE:-}")"
fi

export CONDIV_N_THREADS="${CONDIV_N_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
if [ -n "$SLURM_TASK_SLOTS" ] && [ -z "${CONDIV_MAX_PARALLEL_WORKERS:-}" ]; then
  export CONDIV_MAX_PARALLEL_WORKERS="$SLURM_TASK_SLOTS"
fi

export CONDIV_WORKER_LAUNCH="$WORKER_LAUNCH"
export CONDIV_SYMLAY_LAYER_MANIFEST="$LAYER_MANIFEST_JSON"

RESOLVED_WORKER_LAUNCH="$CONDIV_WORKER_LAUNCH"
if [ "$RESOLVED_WORKER_LAUNCH" = "auto" ] && [ -n "${SLURM_JOB_ID:-}" ]; then
  RESOLVED_WORKER_LAUNCH="srun"
fi

find_latest_checkpoint() {
  local base_dir="$1"
  local ckpt="$base_dir/initial_checkpoint.pkl"
  local all_dirs
  all_dirs=$(ls -d "$base_dir"/epoch_*_minibatch_* 2>/dev/null | sort -r || true)
  for d in $all_dirs; do
    if [ -f "$d/checkpoint.pkl" ]; then
      ckpt="$d/checkpoint.pkl"
      break
    fi
  done
  echo "$ckpt"
}

checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
if [ ! -f "$checkpoint" ]; then
  echo "ERROR: checkpoint not found at $checkpoint"
  echo "Run ./run_init.sh first."
  exit 1
fi

echo "Running Slurm membrane ConDiv restart"
echo "  base dir: $BASE_DIR"
echo "  checkpoint: $checkpoint"
echo "  steps: $RUN_STEPS"
echo "  layer manifest: $LAYER_MANIFEST_JSON"
echo "  slurm job id: ${SLURM_JOB_ID:-N/A}"
echo "  worker launch request: $CONDIV_WORKER_LAUNCH"
echo "  worker launch resolved: $RESOLVED_WORKER_LAUNCH"
echo "  threads/worker: ${CONDIV_N_THREADS}"
echo "  max parallel workers: ${CONDIV_MAX_PARALLEL_WORKERS:-auto}"
echo "  allocated task slots: ${SLURM_TASK_SLOTS:-unknown}"
echo "  convergence grad threshold: ${CONDIV_CONVERGENCE_GRAD_NORM}"
echo "  convergence update threshold: ${CONDIV_CONVERGENCE_UPDATE_NORM}"
echo "  convergence patience: ${CONDIV_CONVERGENCE_PATIENCE}"
echo "  auto resubmit: ${CONDIV_AUTO_RESUBMIT}"
echo "  resubmit count: ${CONDIV_RESUBMIT_COUNT}"
echo "  progress log: ${PROGRESS_LOG_JSONL}"
echo "  status json: ${STATUS_JSON}"

python3 -u "$SCRIPT_DIR/ConDiv_mem.py" restart "$checkpoint" "$RUN_STEPS" 2>&1 | tee -a "$BASE_DIR/run.output"

checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
set +e
python3 "$SCRIPT_DIR/training_control.py" \
  --base-dir "$BASE_DIR" \
  --checkpoint "$checkpoint" \
  --progress-log "$PROGRESS_LOG_JSONL" \
  --status-json "$STATUS_JSON" \
  --grad-threshold "$CONDIV_CONVERGENCE_GRAD_NORM" \
  --update-threshold "$CONDIV_CONVERGENCE_UPDATE_NORM" \
  --patience "$CONDIV_CONVERGENCE_PATIENCE" \
  --run-steps "$RUN_STEPS" \
  --slurm-job-id "${SLURM_JOB_ID:-}" \
  --resubmit-count "$CONDIV_RESUBMIT_COUNT" 2>&1 | tee -a "$BASE_DIR/run.output"
control_rc=${PIPESTATUS[0]}
set -e

if [ "$control_rc" -eq 0 ]; then
  echo "Training converged; no resubmission needed."
  exit 0
fi

if [ "$control_rc" -ne 10 ]; then
  echo "ERROR: training control failed with exit code $control_rc"
  exit "$control_rc"
fi

if [ "$RUN_STEPS" -le 0 ]; then
  echo "RUN_STEPS <= 0; skipping auto-resubmission."
  exit 0
fi

if [ "$CONDIV_AUTO_RESUBMIT" != "1" ]; then
  echo "Auto-resubmission disabled; stopping before convergence."
  exit 0
fi

if [ -z "${SLURM_JOB_ID:-}" ]; then
  echo "Not running inside Slurm; skipping auto-resubmission."
  exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not available; skipping auto-resubmission."
  exit 0
fi

next_resubmit_count=$((CONDIV_RESUBMIT_COUNT + 1))
if [ "$CONDIV_AUTO_MAX_SUBMISSIONS" -gt 0 ] && [ "$next_resubmit_count" -gt "$CONDIV_AUTO_MAX_SUBMISSIONS" ]; then
  echo "Reached CONDIV_AUTO_MAX_SUBMISSIONS=${CONDIV_AUTO_MAX_SUBMISSIONS}; stopping before convergence."
  exit 0
fi

submit_output="$(sbatch --parsable --dependency=afterok:${SLURM_JOB_ID} --export=ALL,CONDIV_RESUBMIT_COUNT=${next_resubmit_count} "$SCRIPT_DIR/run_remote.sh" "$RUN_STEPS")"
echo "Submitted continuation job: $submit_output"
