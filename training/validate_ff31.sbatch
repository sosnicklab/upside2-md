#!/bin/bash
#SBATCH --job-name=ff31-validate
#SBATCH --account=pi-trsosnic
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=%x_%j.out

# Validation of the trained ff3.1, submitted automatically by the last link of the training
# chain (train_gly.sbatch, "target met" branch).
#
#   sbatch validate_ff31.sbatch [checkpoint.pkl]
#
# This job is a launcher, not a simulation: it extracts the force field, checks it, then submits
# the 32 benchmark arms.  Keeping extraction here rather than in the training chain means a
# failed extraction is visible as its own failed job instead of a silent no-op at the end of a
# four-day run.
#
# WHY THE EXTRACTION IS NOT extract_ff.py.  That script writes sidechain.h5 and environment.h5
# only, which sufficed when hb, sheet and the glycine map were frozen.  ff3.1 trains all five, so
# extract_ff31.py is used instead and it calls ConDiv's own expand_param, making the output
# byte-identical to what training writes each step.

set -eo pipefail
T=/project/trsosnic/yinhan/upside2-md-mdw2/training/ff31-gly
UP=/project/trsosnic/yinhan/upside2-md-mdw2
B=/beagle3/trsosnic/yinhan/ff3_benchmark
FF=ff_3.1_trained
OUT=$UP/parameters/$FF

cd "$T"
source ./env.sh

CKPT="${1:-}"
if [ -z "$CKPT" ]; then
    CKPT=$(find run_output -name checkpoint.pkl -path "*/epoch_*/checkpoint.pkl" | sort | tail -1)
fi
[ -n "$CKPT" ] || { echo "no checkpoint found" >&2; exit 1; }
echo "extracting from $CKPT"

python3 "$UP/training/extract_ff31.py" "$CKPT" "$OUT"

# bb_env.dat is not a ConDiv parameter and is never written by training, but upside_config
# requires it.  Without this copy every benchmark arm dies at config build.
cp -n "$UP/parameters/ff_2.1/bb_env.dat" "$OUT/bb_env.dat"

for f in sidechain.h5 environment.h5 hbond.h5 sheet rama.dat bb_env.dat; do
    [ -s "$OUT/$f" ] || { echo "MISSING $OUT/$f" >&2; exit 1; }
done
echo "force field complete in $OUT"

# ---------------------------------------------------------------------------
# The Peng benchmark: 16 proteins, native and de novo.  bench.sbatch self-chains
# per arm, so each of these 32 submissions is the head of its own chain.
# ---------------------------------------------------------------------------
PROTS="alpha3d BBA BBL cspA gpW homeodomain hyp lambda NTL9 NuG2 proteinB proteinG proteinL top7 ubiquitin WWdomain"
mkdir -p "$B/logs"
n=0
for P in $PROTS; do
    for KIND in native denovo; do
        sbatch --parsable --account=pi-trsosnic --partition=broadwl --time=36:00:00 \
               --job-name="b31_${P}_${KIND}" \
               --output="$B/logs/${P}_${KIND}_${FF}_%j.out" \
               --export=ALL,PROT=$P,KIND=$KIND,FF=$FF \
               "$B/bench.sbatch" >/dev/null && n=$((n+1))
    done
done
echo "submitted $n of 32 benchmark arms with FF=$FF"

# ---------------------------------------------------------------------------
# glpG: patch the trained force field into the four midway2 POPE/POPG seeds.
#
# Patched, not rebuilt.  Re-running popepopg_prep.sbatch would redo bilayer packing and
# equilibration and so change the starting structure -- confounding the force field with the
# initial condition on exactly the observable of interest, TM4 helicity.
#
# The seeds are the LIVE ones in popepopg_REMD_mdw2, not the pristine ff2.1 backups, because
# they alone carry `inner_steps = 4` on /input/brownian.  That temperature fix is worth
# TM4 helix fraction 0.893 against 0.346 (findings 3.10a), and it is an HDF5 attribute rather
# than an array, so patching arrays preserves it.  Measured: the live seed differs from the
# pristine ff2.1 one in exactly two arrays, rama_pot and interaction_param, which are ff3.0's
# and are precisely what this overwrites.  ff3.0 is superseded, so overwriting is intended.
#
# THE METHOD GATE runs on the pristine ff2.1 seed, where a round trip CAN succeed (measured
# 3.55e-15).  It cannot succeed on a live seed, because those carry ff3.0 in the two arrays
# above -- a failure there would mean nothing.  So the method is proven once, then applied.
# ---------------------------------------------------------------------------
GM=/project/trsosnic/yinhan/popepopg_REMD_mdw2
PRISTINE=/project/trsosnic/yinhan/popepopg_REMD/seeds/glpG-RKRK-79HIS.up.bak_production_handoff
VARIANTS="glpG-RKRK-79HIS glpG-RKRK-79HIS_S115T glpG-RKRK-79ALA glpG-RKRK-79ALA_S115T"
STAMP=$(date +%Y%m%d-%H%M%S)

echo "--- method gate on a pristine ff_2.1 seed"
python3 "$UP/training/patch_glpg_ff31.py" --ff "$OUT" --seed "$PRISTINE"         --out /tmp/gate_$$.up --verify-roundtrip "$UP/parameters/ff_2.1" || {
    echo "METHOD GATE FAILED - glpG seeds NOT patched" >&2; rm -f /tmp/gate_$$.up; exit 1; }
rm -f /tmp/gate_$$.up

ok=0
for V in $VARIANTS; do
    S="$GM/seeds/$V.up"
    if [ ! -f "$S" ]; then echo "MISSING seed $S"; continue; fi
    cp -p "$S" "$S.bak_pre_ff31_$STAMP"
    echo "--- $V"
    if python3 "$UP/training/patch_glpg_ff31.py" --ff "$OUT"             --seed "$S.bak_pre_ff31_$STAMP" --out "$S"; then
        ok=$((ok+1))
    else
        echo "PATCH FAILED for $V; restoring"; cp -p "$S.bak_pre_ff31_$STAMP" "$S"
    fi
done
echo "patched $ok of 4 glpG seeds in $GM/seeds (backups *.bak_pre_ff31_$STAMP)"

if [ "$ok" -eq 4 ]; then
    for V in $VARIANTS; do
        rm -rf "$GM/$V"                 # stale ff3.0 replicas and block_count
        bash "$GM/submit_remd.sh" "$V" | sed "s/^/  /"
    done
    echo "submitted 4 glpG REMD chains on broadwl (28 replicas each, REMD_T_HI=0.90)"
    echo "criterion: helix stability over time, TM4 above all"
else
    echo "NOT submitting glpG: only $ok of 4 seeds patched" >&2
fi
