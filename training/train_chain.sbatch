#!/bin/bash
#SBATCH --job-name=condiv-train
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=36:00:00
#SBATCH --no-requeue
#SBATCH --output=%x_%j.out

# Self-chaining FF2 dual-target ConDiv training.
#
#   sbatch $(cat <run_dir>/slurm.args) train_chain.sbatch <run_dir> <target_steps> [<max_epochs>]
#
# <run_dir>/slurm.args holds the cluster's flags: account, partition and any excluded nodes, e.g.
#   --account=<account> --partition=<partition> --exclude=<nodes>
# They are given on the command line of every submission, first link, successor and gate alike,
# because a site's submit filter can remap #SBATCH partition directives.
#
# --ntasks must equal the minibatch size (24) and --cpus-per-task the systems per protein (14:
# restrained native, 12 free replicas, SARW); each task is one protein worker dispatched by srun.
#
# The successor is queued BEFORE training starts, as insurance: a wall kill is SIGKILL and cannot
# submit anything, so a link that dies would otherwise end the chain silently.
#
# The link that reaches the target cancels that successor and runs the convergence gate
# (`ConDiv.py gate`, report in gate_step<N>.txt). Converged: the chain stops and lists the
# epoch-end checkpoints; the one to release is chosen by simulating them, since the fixed-point
# test does not say which iterate simulates best, and `ConDiv.py extract` writes it. Not converged:
# the chain trains one more epoch, until <max_epochs> (default: no extension), then stops for
# review. A failure of the gate itself stops the chain.
#
# The driver is the run directory's own copy of ConDiv.py, written at initialisation, so a run
# is never continued by later code.

set -eo pipefail
RUN_DIR="$(cd "${1:?usage: sbatch train_chain.sbatch <run_dir> <target_steps>}" && pwd)"
TARGET="${2:?usage: sbatch train_chain.sbatch <run_dir> <target_steps>}"
MAX_EPOCHS="${3:-0}"
TRAIN_DIR="$(cd "$RUN_DIR/.." && pwd)"     # run dirs live directly under training/
SLURM_ARGS=$(cat "$RUN_DIR/slurm.args")

MB_PER_EPOCH=19            # 456 proteins / 24
STEPS_PER_LINK=54          # a step is ~36 min at dt 0.009 (5vhg's 1242 s at dt 0.015, scaled), so ~32 h in the 36 h wall

cd "$RUN_DIR"
CKPT=$(find run_output -name checkpoint.pkl -path "*/epoch_*/checkpoint.pkl" 2>/dev/null | sort | tail -1)
if [ -z "$CKPT" ]; then
    CKPT=run_output/initial_checkpoint.pkl
    STEP=0
else
    NAME=$(basename "$(dirname "$CKPT")")
    [[ "$NAME" =~ ^epoch_([0-9]+)_minibatch_([0-9]+)$ ]] || { echo "unexpected checkpoint name $NAME" >&2; exit 1; }
    STEP=$(( 10#${BASH_REMATCH[1]} * MB_PER_EPOCH + 10#${BASH_REMATCH[2]} + 1 ))
fi

echo "step $STEP of $TARGET, resuming from $CKPT"
if [ "$STEP" -ge "$TARGET" ]; then
    echo "target already met, nothing to do"
    exit 0
fi

# A link that finishes no step leaves its successor the same resume point. That is right after a
# node failure, but a failure that repeats would resubmit these 336 CPUs without end, so three
# links in a row that start from the same step stop the chain. Delete .chain_starts to retry.
touch .chain_starts
if [ "$(tail -3 .chain_starts | grep -cx "$STEP")" -eq 3 ]; then
    echo "STOPPED: the last three links all started from step $STEP and finished none" >&2
    exit 1
fi
echo "$STEP" >> .chain_starts

# The last link runs only the steps that remain.
REMAINING=$(( TARGET - STEP ))
if [ "$REMAINING" -lt "$STEPS_PER_LINK" ]; then
    STEPS_PER_LINK=$REMAINING
fi
echo "this link runs $STEPS_PER_LINK steps"

SUCCESSOR=$(sbatch --parsable $SLURM_ARGS --dependency=afterany:$SLURM_JOB_ID \
                   "$TRAIN_DIR/train_chain.sbatch" "$RUN_DIR" "$TARGET" "$MAX_EPOCHS") \
    || echo "WARNING: successor was not queued, the chain stops after this link"

source "$RUN_DIR/env.sh"
python3 "$RUN_DIR/run_output/ConDiv.py" restart "$CKPT" "$STEPS_PER_LINK"

[ "$STEPS_PER_LINK" -eq "$REMAINING" ] || exit 0
echo "target $TARGET reached"
[ -n "$SUCCESSOR" ] && scancel "$SUCCESSOR"

set +e
python3 "$RUN_DIR/run_output/ConDiv.py" gate "$RUN_DIR" > "gate_step$TARGET.txt" 2>&1
RC=$?
set -e
cat "gate_step$TARGET.txt"
case $RC in
    0)
        echo "converged; choose among the epoch-end checkpoints by simulation:"
        find run_output -path "*/epoch_*_minibatch_$(( MB_PER_EPOCH - 1 ))/checkpoint.pkl" | sort
        echo "then: python3 run_output/ConDiv.py extract <checkpoint> \$UPSIDE_HOME/parameters/<ff_name>"
        ;;
    3)
        if [ "$TARGET" -ge $(( MAX_EPOCHS * MB_PER_EPOCH )) ]; then
            echo "not converged at step $TARGET; stopped for review"
        else
            echo "not converged; training on to step $(( TARGET + MB_PER_EPOCH ))"
            sbatch $SLURM_ARGS "$TRAIN_DIR/train_chain.sbatch" "$RUN_DIR" \
                   "$(( TARGET + MB_PER_EPOCH ))" "$MAX_EPOCHS"
        fi
        ;;
    *)
        echo "the convergence gate itself failed (exit $RC); stopped" >&2
        exit 1
        ;;
esac
