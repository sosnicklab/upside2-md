#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

pdb_id=${pdb_id:-glpG-RKRK-79HIS} # CHECKME
sim_id=${sim_id:-memb_test} # CHECKME

n_rep=${n_rep:-48} # CHECKME

work_dir=./
input_dir=$work_dir/inputs
output_dir=$work_dir/outputs
run_dir=$output_dir/$sim_id

ref=$input_dir/$pdb_id.up

if [ ! -d "results" ]; then
    mkdir results
fi

for i in `seq 0 $((n_rep-1))` 
do
    traj=$output_dir/$sim_id/$pdb_id.run.$i.up
    python $UPSIDE_HOME/py/get_protection_state.py $input_dir/$pdb_id-HDX.up  $traj results/${pdb_id}_${sim_id}_${i}_PS.npy --residue results/$pdb_id.resid
done
