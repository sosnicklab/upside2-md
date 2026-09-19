#!/bin/bash
# Pack the initial force field into the latent vector and build the minibatches.
#
#   bash init.sh <run_dir>
#
# <run_dir> must already contain:
#   init_param/   environment.h5, sidechain.h5, hbond, sheet   (e.g. copies of parameters/ff_2.1)
#   upside_input/ per-protein .fasta/.initial.pkl/.chi plus rama.dat and rama_reference.pkl
#   pdb_list      the training set manifest
#
# Writes <run_dir>/run_output/initial_checkpoint.pkl.
#
# Check the printed `pack_param residual`: it must be far below the 1.6e-4 gate (ff_2.1 packs to
# ~3e-30).  A large residual means the latent vector does not reproduce the input force field and
# nothing downstream is meaningful.

set -euo pipefail
RUN_DIR="${1:?usage: bash init.sh <run_dir>}"
TRAIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$RUN_DIR"
source "$TRAIN_DIR/env.sh"

python3 "$TRAIN_DIR/ConDiv.py" initialize ./init_param ./upside_input ./pdb_list ./run_output
echo "initialised: $RUN_DIR/run_output/initial_checkpoint.pkl"
