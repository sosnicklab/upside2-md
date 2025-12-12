#!/bin/bash

# Define where you are working
X=$(pwd)
Y=test_00

# Initialization
# If you haven't initialized yet, uncomment the lines below:
# mode=initialize
# initial_param=init_param       # Folder with initial parameters
# upside_input_dir=upside_input  # Folder with protein inputs
# pdb_list=pdb_list              # File with list of PDBs
# python3 ConDiv.py $mode $initial_param $upside_input_dir $pdb_list $X/$Y | tee $X/$Y_init.output

# Training Loop
# Runs locally on the M1 (srun removed)
mode=restart
checkpoint=$X/$Y/initial_checkpoint.pkl
step=40

# Run
python3 ConDiv.py $mode $checkpoint $step | tee -a $X/$Y.output
