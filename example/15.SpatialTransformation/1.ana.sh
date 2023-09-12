#!/bin/bash

source ../../source.sh

pdb_id=chig
sim_id=simple_test

base_dir=./
input_dir=$base_dir/inputs
output_dir=$base_dir/outputs
run_dir=$output_dir/$sim_id

ref=$input_dir/$pdb_id.up
traj=$output_dir/$sim_id/$pdb_id.run.up

mkdir results

# calculate the rmsd
python calc_rmsd.py $ref $traj > results/${pdb_id}_$sim_id.rmsd

# convert .up/.h5 file to .vtf 
python $UPSIDE_HOME/py/extract_vtf.py $traj results/${pdb_id}_$sim_id.vtf
