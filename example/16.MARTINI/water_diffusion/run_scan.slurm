#!/bin/bash
#SBATCH --job-name=water_diff_scan
#SBATCH --output=water_scan_%A_%a.out
#SBATCH --error=water_scan_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --array=1-651

TASK_LIST="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/water_diffusion/tasks.txt"
SCRIPT_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")

if [ -f "$SCRIPT_PATH" ]; then
    echo "Running task $SLURM_ARRAY_TASK_ID: $SCRIPT_PATH"
    "$SCRIPT_PATH"
else
    echo "Error: Script not found at $SCRIPT_PATH"
    exit 1
fi
