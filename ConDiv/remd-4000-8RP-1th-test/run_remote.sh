#!/bin/bash
#SBATCH --job-name=condiv-train
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -euo pipefail

WORK_DIR="$(pwd)"
VENV_PYTHON="$WORK_DIR/venv/bin/python3"

PROJECT_ROOT="$WORK_DIR/../.."
UPSIDE_LIB="$PROJECT_ROOT/obj"
UPSIDE_PY="$PROJECT_ROOT/py"

if [ ! -x "$VENV_PYTHON" ]; then
    echo "ERROR: venv python not found at $VENV_PYTHON"
    echo "Run setup_venv.sh first."
    exit 1
fi

source "$WORK_DIR/venv/bin/activate"
export PYTHONPATH="$WORK_DIR:$UPSIDE_LIB:$UPSIDE_PY:$PROJECT_ROOT/src:$PYTHONPATH"

mode=restart
step="${1:-${RUN_STEPS:-152}}"

case "$step" in
    ''|*[!0-9]*)
        echo "ERROR: step must be a positive integer, got '$step'"
        exit 1
        ;;
esac

if [ "$step" -le 0 ]; then
    echo "ERROR: step must be > 0, got '$step'"
    exit 1
fi

checkpoint="$WORK_DIR/test_00/initial_checkpoint.pkl"
all_dirs=$(ls -d "$WORK_DIR/test_00"/epoch_*_minibatch_* 2>/dev/null | sort -r || true)

for d in $all_dirs; do
    if [ -f "$d/checkpoint.pkl" ]; then
        echo "--> Found Valid Resume Point: $d"
        checkpoint="$d/checkpoint.pkl"
        break
    fi
done

if [ ! -f "$checkpoint" ]; then
    echo "ERROR: No checkpoint found. Run ./run_init.sh on login node first."
    exit 1
fi

echo "Running ConDiv on Slurm with: python3 (from venv)"
echo "Python executable: $VENV_PYTHON"
echo "Loading Checkpoint: $checkpoint"
echo "Restart iterations (step): $step"
echo "Slurm job id: ${SLURM_JOB_ID:-N/A}"

"$VENV_PYTHON" -u ConDiv.py "$mode" "$checkpoint" "$step" 2>&1 | tee -a "$WORK_DIR/run.output"
