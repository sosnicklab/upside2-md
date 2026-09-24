#!/bin/bash
#SBATCH --job-name=condiv-train
#SBATCH --account=pi-trsosnic
#SBATCH --partition=broadwl
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=36:00:00
#SBATCH --no-requeue
#SBATCH --exclude=midway2-0003,midway2-[0010-0011],midway2-[0342-0345]
#SBATCH --output=%x_%j.out

# Self-chaining FF2 dual-target ConDiv training.
#
#   sbatch train_chain.sbatch <run_dir> <target_steps>
#
# --ntasks must equal the minibatch size (24) and --cpus-per-task the systems per protein (14:
# restrained native, 12 free replicas, SARW); each task is one protein worker dispatched by srun.
#
# The successor is queued BEFORE training starts, as insurance: a wall kill is SIGKILL and cannot
# submit anything, so a link that dies would otherwise end the chain silently. The link that
# reaches the target cancels that successor and submits <run_dir>/after_training.sbatch, if the
# run has one. It is submitted from here rather than from a successor because a successor asks for
# this job's 336 CPUs and waited three hours in the queue to run a single sbatch.
#
# The driver is the run directory's own copy of ConDiv.py, written at initialisation, so a run
# is never continued by later code.
#
# --exclude carries midway2-[0010-0011], an adjacent pair that produced three NODE_FAILs in one
# campaign while other nodes ran clean for 14 h.  If a NODE_FAIL happens outside that pair, say so
# rather than just adding the node: the bad-pair explanation would then be wrong.

set -eo pipefail
RUN_DIR="$(cd "${1:?usage: sbatch train_chain.sbatch <run_dir> <target_steps>}" && pwd)"
TARGET="${2:?usage: sbatch train_chain.sbatch <run_dir> <target_steps>}"
TRAIN_DIR="$(cd "$RUN_DIR/.." && pwd)"     # run dirs live directly under training/

MB_PER_EPOCH=19            # 456 proteins / 24
STEPS_PER_LINK=90          # a step is ~21 min (5vhg, 150 res, 1242 s measured), so ~32 h in the 36 h wall

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

# The last link runs only the steps that remain.
REMAINING=$(( TARGET - STEP ))
if [ "$REMAINING" -lt "$STEPS_PER_LINK" ]; then
    STEPS_PER_LINK=$REMAINING
fi
echo "this link runs $STEPS_PER_LINK steps"

SUCCESSOR=$(sbatch --parsable --dependency=afterany:$SLURM_JOB_ID \
                   "$TRAIN_DIR/train_chain.sbatch" "$RUN_DIR" "$TARGET") \
    || echo "WARNING: successor was not queued, the chain stops after this link"

source "$RUN_DIR/env.sh"
python3 "$RUN_DIR/run_output/ConDiv.py" restart "$CKPT" "$STEPS_PER_LINK"

if [ "$STEPS_PER_LINK" -eq "$REMAINING" ]; then
    echo "target $TARGET reached"
    [ -n "$SUCCESSOR" ] && scancel "$SUCCESSOR"
    if [ -f "$RUN_DIR/after_training.sbatch" ]; then
        sbatch "$RUN_DIR/after_training.sbatch" \
            || echo "WARNING: after_training.sbatch was NOT submitted; submit it by hand"
    fi
fi
