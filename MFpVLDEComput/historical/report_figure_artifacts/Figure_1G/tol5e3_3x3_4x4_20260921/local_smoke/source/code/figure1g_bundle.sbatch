#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=1-7%2
#SBATCH --job-name=fig1g_bundle
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
module purge
module load matlab/2025b
sizes=(4 6 8 10 20 30 40)
field=${sizes[$((SLURM_ARRAY_TASK_ID-1))]}
matlab -batch "addpath('$RUN_ROOT/code'); build_figure1g_field_bundle('$RUN_ROOT',$field)"
