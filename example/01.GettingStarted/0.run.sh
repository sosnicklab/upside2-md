#!/bin/bash

source ../../source.sh
upside=$UPSIDE_HOME/obj/upside

#----------------------------------------------------------------------
## General Settings and Path
#----------------------------------------------------------------------

pdb_id=chig
pdb_dir=./pdb
sim_id=simple_test

is_native=$true
ff=ff_2.1
T=0.8
duration=1000
frame_interval=50
base_dir=./


input_dir=$base_dir/inputs
output_dir=$base_dir/outputs
run_dir=$output_dir/$sim_id

for direc in $input_dir $output_dir $run_dir
do
    if [ ! -d $direc ]; then
      mkdir $direc
    fi
done

#----------------------------------------------------------------------
## Generate Upside readable initial structure (and fasta) from PDB 
#----------------------------------------------------------------------

echo "Initial structure gen..."
python $UPSIDE_HOME/py/PDB_to_initial_structure.py $pdb_dir/$pdb_id.pdb $input_dir/$pdb_id --record-chain-breaks
echo ""

#----------------------------------------------------------------------
## Configure
#----------------------------------------------------------------------

# parameters
param_dir_base="$UPSIDE_HOME/parameters/"
param_dir_common="$param_dir_base/common/"
param_dir_ff="$param_dir_base/$ff"

echo "Configuring..."
config_base=$input_dir/$pdb_id.up
kwargs="--fasta=$input_dir/$pdb_id.fasta"
kwargs="$kwargs --output=$config_base"
kwargs="$kwargs --rama-library=$param_dir_common/rama.dat"
kwargs="$kwargs --rama-sheet-mixing-energy=$param_dir_ff/sheet"
kwargs="$kwargs --hbond-energy=$param_dir_ff/hbond.h5"
kwargs="$kwargs --reference-state-rama=$param_dir_common/rama_reference.pkl"
kwargs="$kwargs --rotamer-placement=$param_dir_ff/sidechain.h5"
kwargs="$kwargs --rotamer-interaction=$param_dir_ff/sidechain.h5"
kwargs="$kwargs --dynamic-rotamer-1body"
kwargs="$kwargs --environment-potential=$param_dir_ff/environment.h5"
kwargs="$kwargs --bb-environment-potential=$param_dir_ff/bb_env.dat"
kwargs="$kwargs --chain-break-from-file=$input_dir/$pdb_id.chain_breaks"

if [ is_native ]; then
    kwargs="$kwargs --initial-structure=$input_dir/$pdb_id.initial.npy"
fi

echo "Config commandline options:"
echo $UPSIDE_HOME/py/upside_config.py $kwargs

$UPSIDE_HOME/py/upside_config.py $kwargs

#----------------------------------------------------------------------
## Run Settings
#----------------------------------------------------------------------

#seed=`echo $RANDOM`
seed=1

h5_file=$run_dir/$pdb_id.run.up
log_file=$run_dir/$pdb_id.run.log

cp $config_base $h5_file

echo "Running..."
# run simulation
$upside --frame-interval $frame_interval --temperature $T --duration $duration --seed $seed $h5_file | tee $log_file
