#!/bin/bash
# Release a trained force field and launch its validation.
#
#   bash validate_ff.sh <run_dir> <ff_name>
#
# Run by <run_dir>/after_training.sbatch, which train_chain.sbatch submits when the run reaches its
# target. A launcher, not a simulation: it
#   1. extracts the newest checkpoint with extract_ff.py (the run's own expand_param);
#   2. backs up and overwrites parameters/<ff_name> in this tree and in the shared /beagle3
#      deployment the benchmark reads, md5-verifying the copy;
#   3. moves the superseded benchmark runs of <ff_name> aside and submits the 32 Peng arms;
#   4. proves the glpG patch method on a pristine ff_2.1 seed, patches the four live seeds, clears
#      their superseded replicas and submits the four REMD chains.
# Every step that fails stops the launcher: nothing downstream runs on a half-released force field.

set -eo pipefail
RUN_DIR="$(cd "${1:?usage: validate_ff.sh <run_dir> <ff_name>}" && pwd)"
FF="${2:?usage: validate_ff.sh <run_dir> <ff_name>}"
UP=/project/trsosnic/yinhan/upside2-md-mdw2
DEPLOY_ROOT=/beagle3/trsosnic/yinhan/upside2-md
B=/beagle3/trsosnic/yinhan/ff3_benchmark
GM=/project/trsosnic/yinhan/popepopg_REMD_mdw2
PRISTINE=/project/trsosnic/yinhan/popepopg_REMD/seeds/glpG-RKRK-79HIS.up.bak_production_handoff
VARIANTS="glpG-RKRK-79HIS glpG-RKRK-79HIS_S115T glpG-RKRK-79ALA glpG-RKRK-79ALA_S115T"
STAMP=$(date +%Y%m%d-%H%M%S)

source "$RUN_DIR/env.sh"
CKPT=$(find "$RUN_DIR/run_output" -name checkpoint.pkl -path "*/epoch_*/checkpoint.pkl" | sort | tail -1)
[ -n "$CKPT" ] || { echo "no checkpoint under $RUN_DIR/run_output" >&2; exit 1; }

# --- 1-2. extract and release -------------------------------------------------------------------
NEW="$RUN_DIR/release_$STAMP"
python3 "$UP/training/extract_ff.py" "$CKPT" "$NEW"
for f in sidechain.h5 environment.h5 bb_env.dat hbond.h5 sheet rama.dat; do
    [ -s "$NEW/$f" ] || { echo "MISSING $NEW/$f" >&2; exit 1; }
done

release() {                      # release <parameters_root>
    local dst="$1/$FF"
    if [ -d "$dst" ]; then
        mkdir -p "$1/../backup"
        mv "$dst" "$1/../backup/${FF}_superseded_$STAMP"
        echo "backed up $dst -> $1/../backup/${FF}_superseded_$STAMP"
    fi
    mkdir -p "$dst"
    cp -p "$NEW"/* "$dst"/
    (cd "$NEW" && md5sum *) | (cd "$dst" && md5sum -c --quiet -) \
        || { echo "copy to $dst does not match $NEW" >&2; exit 1; }
    echo "released $FF to $dst, md5 verified"
}
release "$UP/parameters"
release "$DEPLOY_ROOT/parameters"

# --- 3. the Peng benchmark: 16 proteins, native and de novo, each arm self-chaining ------------
if ls -d "$B/runs/"*"_$FF" >/dev/null 2>&1; then
    mkdir -p "$B/runs_superseded/${FF}_$STAMP"
    mv "$B/runs/"*"_$FF" "$B/runs_superseded/${FF}_$STAMP/"
    echo "moved superseded $FF benchmark runs to $B/runs_superseded/${FF}_$STAMP"
fi
PROTS="alpha3d BBA BBL cspA gpW homeodomain hyp lambda NTL9 NuG2 proteinB proteinG proteinL top7 ubiquitin WWdomain"
mkdir -p "$B/logs"
n=0
for P in $PROTS; do
    for KIND in native denovo; do
        sbatch --parsable --account=pi-trsosnic --partition=broadwl --time=36:00:00 \
               --job-name="b_${FF}_${P}_${KIND}" \
               --output="$B/logs/${P}_${KIND}_${FF}_%j.out" \
               --export=ALL,PROT=$P,KIND=$KIND,FF=$FF \
               "$B/bench.sbatch" >/dev/null && n=$((n+1))
    done
done
echo "submitted $n of 32 benchmark arms with FF=$FF"
[ "$n" -eq 32 ] || { echo "not every benchmark arm was submitted" >&2; exit 1; }

# --- 4. glpG --------------------------------------------------------------------------------------
# Patched, not rebuilt: re-running the bilayer preparation would change the starting structure and
# confound the force field with the initial condition on TM4 helicity. The live seeds alone carry
# `inner_steps = 4` on /input/brownian (TM4 helix 0.893 against 0.346), an attribute that a patch of
# arrays preserves. The method is proven on a pristine ff_2.1 seed, where a round trip can succeed.
OUT="$UP/parameters/$FF"
echo "--- method gate on a pristine ff_2.1 seed"
python3 "$UP/training/patch_glpg.py" --ff "$OUT" --seed "$PRISTINE" --out "$RUN_DIR/gate_$STAMP.up" \
        --verify-roundtrip "$UP/parameters/ff_2.1"
rm -f "$RUN_DIR/gate_$STAMP.up"

for V in $VARIANTS; do
    S="$GM/seeds/$V.up"
    [ -f "$S" ] || { echo "MISSING seed $S" >&2; exit 1; }
    cp -p "$S" "$S.bak_pre_${FF}_$STAMP"
    echo "--- $V"
    python3 "$UP/training/patch_glpg.py" --ff "$OUT" --seed "$S.bak_pre_${FF}_$STAMP" --out "$S"
done
echo "patched 4 glpG seeds in $GM/seeds (backups *.bak_pre_${FF}_$STAMP)"
for V in $VARIANTS; do
    rm -rf "$GM/$V"                 # replicas of the superseded force field
    bash "$GM/submit_remd.sh" "$V" | sed "s/^/  /"
done
echo "submitted 4 glpG REMD chains on broadwl; criterion: helix stability over time, TM4 above all"
