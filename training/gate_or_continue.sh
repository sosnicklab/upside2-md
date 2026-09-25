#!/bin/bash
# Decide what a training run does once it reaches its target: validate, or train one more epoch.
#
#   bash gate_or_continue.sh <run_dir> <ff_name> <max_epochs>
#
# Run by <run_dir>/after_training.sbatch, which train_chain.sbatch submits when the run reaches its
# target. convergence_gate.py judges the last full epoch:
#   exit 0 (every group at a fixed point)  -> validate_ff.sh releases <ff_name> and starts validation
#   exit 3 (a group still pulled)          -> train_chain.sbatch again, target one epoch further,
#                                             unless <max_epochs> is reached: then stop for review
#   anything else (the gate itself failed) -> stop; nothing is released and training does not go on
# Each verdict is kept as gate_step<N>.txt in the run directory.

set -eo pipefail
RUN_DIR="$(cd "${1:?usage: gate_or_continue.sh <run_dir> <ff_name> <max_epochs>}" && pwd)"
FF="${2:?usage: gate_or_continue.sh <run_dir> <ff_name> <max_epochs>}"
MAX_EPOCHS="${3:?usage: gate_or_continue.sh <run_dir> <ff_name> <max_epochs>}"
TRAIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MB_PER_EPOCH=19

cd "$RUN_DIR"
source ./env.sh
CKPT=$(find run_output -name checkpoint.pkl -path "*/epoch_*/checkpoint.pkl" | sort | tail -1)
NAME=$(basename "$(dirname "$CKPT")")
[[ "$NAME" =~ ^epoch_([0-9]+)_minibatch_([0-9]+)$ ]] || { echo "unexpected checkpoint $NAME" >&2; exit 1; }
STEP=$(( 10#${BASH_REMATCH[1]} * MB_PER_EPOCH + 10#${BASH_REMATCH[2]} + 1 ))

REPORT="gate_step$STEP.txt"
set +e
python3 "$TRAIN_DIR/convergence_gate.py" . > "$REPORT" 2>&1
RC=$?
set -e
cat "$REPORT"

case $RC in
    0)
        echo "step $STEP: converged, releasing $FF"
        bash "$TRAIN_DIR/validate_ff.sh" "$RUN_DIR" "$FF"
        ;;
    3)
        if [ "$STEP" -ge $(( MAX_EPOCHS * MB_PER_EPOCH )) ]; then
            echo "step $STEP: NOT converged after $MAX_EPOCHS epochs; stopped for review, nothing released"
        else
            NEXT=$(( STEP + MB_PER_EPOCH ))
            echo "step $STEP: not converged, training on to step $NEXT"
            sbatch "$TRAIN_DIR/train_chain.sbatch" "$RUN_DIR" "$NEXT"
        fi
        ;;
    *)
        echo "step $STEP: the convergence gate itself failed (exit $RC); stopped, nothing released" >&2
        exit 1
        ;;
esac
