#!/bin/bash
#SBATCH --job-name=condiv-train
#SBATCH --account=pi-trsosnic
#SBATCH --partition=broadwl
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=36:00:00
#SBATCH --no-requeue
#SBATCH --exclude=midway2-0003,midway2-[0010-0011],midway2-[0342-0345]
#SBATCH --output=%x_%j.out

# Self-chaining ConDiv training.
#
#   sbatch train_chain.sbatch <run_dir> [target_steps]
#
# --ntasks must equal the minibatch size (12); each task is one protein worker dispatched by srun.
#
# The successor is queued BEFORE training starts.  A wall kill is SIGKILL and cannot submit
# anything, so a link that runs out of time would otherwise end the chain silently.
#
# --exclude carries midway2-[0010-0011], an adjacent pair that produced three NODE_FAILs in one
# campaign while other nodes ran clean for 14 h.  If a NODE_FAIL happens outside that pair, say so
# rather than just adding the node: the bad-pair explanation would then be wrong.

set -eo pipefail
RUN_DIR="${1:?usage: sbatch train_chain.sbatch <run_dir> [target_steps]}"
TARGET="${2:-500}"

RUN_DIR="$(cd "$RUN_DIR" && pwd)"
TRAIN_DIR="$(cd "$RUN_DIR/.." && pwd)"     # run dirs live directly under training/

MB_PER_EPOCH=38
STEPS_PER_LINK=70          # ~26 min/step with sheet training on, so ~30 h inside the 36 h wall

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
    echo "target met, chain ends here"
    exit 0
fi

sbatch --dependency=afterany:$SLURM_JOB_ID "$TRAIN_DIR/train_chain.sbatch" "$RUN_DIR" "$TARGET" \
    || echo "WARNING: successor was not queued, the chain stops after this link"

source "$TRAIN_DIR/env.sh"
export UPSIDE_TRAIN_NTHREADS=$SLURM_CPUS_PER_TASK
python3 "$TRAIN_DIR/ConDiv.py" restart "$CKPT" "$STEPS_PER_LINK"
