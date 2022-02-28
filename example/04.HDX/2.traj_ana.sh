#!/bin/bash

source ../../source.sh

pdb_id=EHEE_rd2_0005
sim_id=REMD

n_rep=16

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
    python $UPSIDE_HOME/py/get_info_from_upside_traj.py  $traj results/${pdb_id}_${sim_id}_$i
    python $UPSIDE_HOME/py/extract_vtf.py                $traj results/${pdb_id}_${sim_id}_$i.vtf
done
