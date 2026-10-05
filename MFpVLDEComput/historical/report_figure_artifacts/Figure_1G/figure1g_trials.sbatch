#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --job-name=fig1g_scale
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
module purge
module load matlab/2025b
task=${SLURM_ARRAY_TASK_ID}
repeat=$(( (task-1)/7 + 1 ))
sizeIndex=$(( (task-1)%7 ))
sizes=(4 6 8 10 20 30 40)
field=${sizes[$sizeIndex]}
matlab -batch "addpath('$RUN_ROOT/code'); run_figure1g_field_trial('$RUN_ROOT',$field,$repeat)"
