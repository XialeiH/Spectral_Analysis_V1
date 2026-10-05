#!/bin/bash
#SBATCH --job-name=f1g_snn3_trial
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g_snn3_trial('$OUTPUT_ROOT/snn_repo','$OUTPUT_ROOT/snn3_bundle.mat','$OUTPUT_ROOT/results','$SLURM_ARRAY_TASK_ID')"
