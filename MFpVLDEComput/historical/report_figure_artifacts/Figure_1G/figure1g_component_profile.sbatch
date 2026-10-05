#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --array=1-35%6
#SBATCH --job-name=fig1g_prof
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
module purge
module load matlab/2025b
task=$SLURM_ARRAY_TASK_ID
repeat=$(( (task-1)/7 + 1 ))
index=$(( (task-1)%7 ))
sizes=(4 6 8 10 20 30 40)
field=${sizes[$index]}
matlab -batch "addpath('$SUITE_CODE'); run_figure1g_component_profile_full('$SOURCE_ROOT','$SUITE_ROOT',$field,$repeat)"
