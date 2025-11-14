#!/bin/bash

# Source environment variables and virtual environment
source ../../source.sh
source ../../.venv/bin/activate

# Path to upside executable
upside=$UPSIDE_HOME/obj/upside

#----------------------------------------------------------------------
## General Settings and Path for 1UBQ
#----------------------------------------------------------------------

pdb_id=1UBQ
pdb_dir=./pdb
sim_id=1UBQ_local_test

is_native=true
ff=ff_2.1
T=0.8
duration=1000  # Short duration for local testing (1000 steps)
frame_interval=100  # Output frame every 100 steps
base_dir=./

# Create directories if they don't exist
input_dir=$base_dir/inputs
output_dir=$base_dir/outputs
run_dir=$output_dir/$sim_id

for direc in $input_dir $output_dir $run_dir
do
    if [ ! -d $direc ]; then
        mkdir -p $direc
    fi
done

#----------------------------------------------------------------------
## Generate Upside readable initial structure (and fasta) from PDB
#----------------------------------------------------------------------

echo "Generating initial structure from PDB..."
python $UPSIDE_HOME/py/PDB_to_initial_structure.py $pdb_dir/$pdb_id.pdb $input_dir/$pdb_id --record-chain-breaks
echo ""

#----------------------------------------------------------------------
## Configure simulation
#----------------------------------------------------------------------

echo "Configuring simulation..."
config_base=$input_dir/$pdb_id.up

# Common parameters directory
param_dir_base="$UPSIDE_HOME/parameters/"
param_dir_common="$param_dir_base/common/"
param_dir_ff="$param_dir_base/$ff"

# Build configuration command
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

# Use native structure as initial conformation if requested
if [ "$is_native" = true ]; then
    kwargs="$kwargs --initial-structure=$input_dir/$pdb_id.initial.npy"
fi

# Run configuration
$UPSIDE_HOME/py/upside_config.py $kwargs
echo ""

#----------------------------------------------------------------------
## Run simulation locally (no slurm)
#----------------------------------------------------------------------

echo "Running simulation locally..."
h5_file=$run_dir/$pdb_id.run.up
log_file=$run_dir/$pdb_id.run.log

# Copy configuration to run directory
cp $config_base $h5_file

# Set a random seed
seed=$RANDOM

# Run upside locally (no slurm submission)
$upside --frame-interval $frame_interval --temperature $T --duration $duration --seed $seed $h5_file | tee $log_file

echo ""
echo "Converting .up output to VTF format..."
vtf_file=$run_dir/$pdb_id.run.vtf
python3 $UPSIDE_HOME/py/extract_vtf.py $h5_file $vtf_file

echo ""
echo "Simulation complete!"
echo "Output files can be found in: $run_dir"
echo "VTF trajectory: $vtf_file"