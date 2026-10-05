#!/bin/bash
#SBATCH --job-name=f1b_dnn_validate
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --array=1-100%16
set -euo pipefail
module purge
module load matlab/2025b
NN_THREADS=$(cat "$OUTPUT_ROOT/selected_threads.txt")
matlab -batch "addpath('$OUTPUT_ROOT/code'); validate_dnn_thread_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$SLURM_ARRAY_TASK_ID','$NN_THREADS')"
