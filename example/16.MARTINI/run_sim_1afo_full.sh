#!/bin/bash
#SBATCH --job-name=1afo_hybrid_full
#SBATCH --account=pi-trsosnic
#SBATCH --partition=broadwl
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# Submit from midway2: sbatch $(cat /project/trsosnic/yinhan/slurm/midway2.args) run_sim_1afo_full.sh

set -euo pipefail

source /beagle3/trsosnic/yinhan/upside2-md/env_shared.sh
export UPSIDE_SKIP_SOURCE_SH=1
export PYTHONUNBUFFERED=1

export PDB_ID="${PDB_ID:-1afo}"
export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1afo_hybrid_full}"
export RUN_DIR="${RUN_DIR:-outputs/martini_1afo_hybrid_full}"
export PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/1AFO.pdb}"
export AUTO_CONTINUE_STAGE_70="${AUTO_CONTINUE_STAGE_70:-1}"

exec bash "${UPSIDE_HOME}/example/16.MARTINI/run_sim_hybrid.sh" "$@"
