#!/bin/bash
#SBATCH --job-name=f1b_dnn_threads
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --array=1-50%16
set -euo pipefail
module purge
module load matlab/2025b
THREADS=(1 2 4 8 16)
INDEX=$((SLURM_ARRAY_TASK_ID-1))
# Interleave thread settings so time-of-day drift is not tied to thread count.
NN_THREADS=${THREADS[INDEX%5]}
REPEAT=$((INDEX/5+1))
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_dnn_thread_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','$NN_THREADS','probe')"
