#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

source "$SCRIPT_DIR/venv/bin/activate"
export PYTHONPATH="${PYTHONPATH:-}"
export CPLUS_INCLUDE_PATH="${CPLUS_INCLUDE_PATH:-}"
export LIBRARY_PATH="${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "$PROJECT_ROOT/source.sh"

export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT/py:$PROJECT_ROOT/obj:$PROJECT_ROOT/src:${PYTHONPATH:-}"
export CONDIV_PROJECT_ROOT="$PROJECT_ROOT"

BASE_DIR="${BASE_DIR:-$SCRIPT_DIR/validate_test_00}"
INIT_PARAM_DIR="${INIT_PARAM_DIR:-$SCRIPT_DIR/init_param}"
PROTEIN_DIR="${PROTEIN_DIR:-$SCRIPT_DIR/upside_input}"
PDB_LIST="${PDB_LIST:-$SCRIPT_DIR/pdb_list}"
ROUNDS="${ROUNDS:-1}"
STEPS_PER_ROUND="${STEPS_PER_ROUND:-1}"
FD_SAMPLES="${FD_SAMPLES:-6}"
FD_EPS="${FD_EPS:-1e-3}"
FD_REL_MEDIAN_THRESHOLD="${FD_REL_MEDIAN_THRESHOLD:-2e-1}"
FD_REL_MAX_THRESHOLD="${FD_REL_MAX_THRESHOLD:-8e-1}"

# Smoke reduction knobs (checkpoint rewrite before restart)
SIM_TIME_OVERRIDE="${SIM_TIME_OVERRIDE:-200}"
MINIBATCH_OVERRIDE="${MINIBATCH_OVERRIDE:-1}"

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

prepare_smoke_checkpoint() {
  local src_ckpt="$1"
  local out_ckpt="$2"
  CHECKPOINT_IN="$src_ckpt" CHECKPOINT_OUT="$out_ckpt" CONDIV_SCRIPT="$SCRIPT_DIR/ConDiv.py" SIM_TIME_OVERRIDE="$SIM_TIME_OVERRIDE" MINIBATCH_OVERRIDE="$MINIBATCH_OVERRIDE" python3 - <<'PY'
import importlib.util
import os
import pickle
import sys

ckpt_in = os.environ['CHECKPOINT_IN']
ckpt_out = os.environ['CHECKPOINT_OUT']
condiv_script = os.environ['CONDIV_SCRIPT']
sim_time = float(os.environ.get('SIM_TIME_OVERRIDE', '200'))
minibatch_n = int(os.environ.get('MINIBATCH_OVERRIDE', '1'))

spec = importlib.util.spec_from_file_location('ConDiv', condiv_script)
if spec is None or spec.loader is None:
    raise RuntimeError(f'Unable to load ConDiv module from {condiv_script}')
mod = importlib.util.module_from_spec(spec)
sys.modules['ConDiv'] = mod
spec.loader.exec_module(mod)

class _MappedUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == '__main__' and hasattr(mod, name):
            return getattr(mod, name)
        return super().find_class(module, name)

with open(ckpt_in, 'rb') as fh:
    state = _MappedUnpickler(fh).load()

state['sim_time'] = sim_time
if minibatch_n > 0:
    all_items = [item for mb in state['minibatches'] for item in mb]
    all_items = all_items[:minibatch_n]
    state['minibatches'] = [all_items]
    state['n_prot'] = len(all_items)
    state['i_mb'] = 0
    state['epoch'] = 0

with open(ckpt_out, 'wb') as fh:
    pickle.dump(state, fh, -1)

mb = len(state['minibatches'][0]) if state.get('minibatches') else 0
print(f'Prepared smoke checkpoint: sim_time={sim_time}, minibatch={mb}')
PY
}

mkdir -p "$BASE_DIR"

if [ ! -f "$BASE_DIR/initial_checkpoint.pkl" ]; then
  echo "Initializing ConDiv checkpoint"
  python3 "$SCRIPT_DIR/ConDiv.py" initialize "$INIT_PARAM_DIR" "$PROTEIN_DIR" "$PDB_LIST" "$BASE_DIR"
fi

SMOKE_CKPT="$BASE_DIR/smoke_checkpoint.pkl"
prepare_smoke_checkpoint "$BASE_DIR/initial_checkpoint.pkl" "$SMOKE_CKPT"

echo "Running ConDiv validation rounds"
echo "  base dir: $BASE_DIR"
echo "  rounds: $ROUNDS"
echo "  steps/round: $STEPS_PER_ROUND"
echo "  sim_time override: $SIM_TIME_OVERRIDE"
echo "  minibatch override: $MINIBATCH_OVERRIDE"

checkpoint="$SMOKE_CKPT"
for round in $(seq 1 "$ROUNDS"); do
  echo "=== ROUND $round / $ROUNDS ==="
  python3 -u "$SCRIPT_DIR/ConDiv.py" restart "$checkpoint" "$STEPS_PER_ROUND" 2>&1 | tee -a "$BASE_DIR/round_${round}.log"

  checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
  report="$BASE_DIR/gradient_round_${round}.json"
  python3 "$SCRIPT_DIR/check_condiv_gradient.py" \
    --checkpoint "$checkpoint" \
    --report "$report" \
    --fd-samples "$FD_SAMPLES" \
    --fd-eps "$FD_EPS" \
    --rel-median-threshold "$FD_REL_MEDIAN_THRESHOLD" \
    --rel-max-threshold "$FD_REL_MAX_THRESHOLD"
done

echo "All validation rounds completed successfully"
