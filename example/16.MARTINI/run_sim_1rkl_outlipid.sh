#!/bin/bash
#SBATCH --job-name=1rkl_outlipid
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Outside-of-bilayer 1RKL start:
# - place the protein above the upper leaflet,
# - rotate it into a laid-flat orientation,
# - enlarge the box so unfolding has more room than the embedded workflow.
export RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_outlipid}"
export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1rkl_outlipid}"
export PROTEIN_PLACEMENT_MODE="${PROTEIN_PLACEMENT_MODE:-outside-top}"
export PROTEIN_ORIENTATION_MODE="${PROTEIN_ORIENTATION_MODE:-lay-flat}"
export PROTEIN_SURFACE_GAP="${PROTEIN_SURFACE_GAP:-6.0}"
export XY_SCALE="${XY_SCALE:-1.35}"
export BOX_PADDING_XY="${BOX_PADDING_XY:-20.0}"
export BOX_PADDING_Z="${BOX_PADDING_Z:-50.0}"
export PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-5.0}"
export PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-5.0}"
export PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-10.0}"

exec "${SCRIPT_DIR}/run_sim_1rkl.sh" "$@"
