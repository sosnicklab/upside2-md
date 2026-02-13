#!/bin/bash
#SBATCH --job-name=bilayer_diff_scan
#SBATCH --output=bilayer_scan_%A_%a.out
#SBATCH --error=bilayer_scan_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --array=1-24


TASK_LIST="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/lipid_diffusion/tasks.txt"
SCRIPT_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")

if [ -f "$SCRIPT_PATH" ]; then
    echo "Running task $SLURM_ARRAY_TASK_ID: $SCRIPT_PATH"
    "$SCRIPT_PATH"
else
    echo "Error: Script not found at $SCRIPT_PATH"
    exit 1
fi
