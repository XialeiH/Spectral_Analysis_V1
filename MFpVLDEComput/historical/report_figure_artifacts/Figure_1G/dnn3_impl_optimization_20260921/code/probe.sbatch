#!/bin/bash
#SBATCH --job-name=f1b_dnn_impl
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --array=1-80%16
set -euo pipefail
module purge
module load matlab/2025b
VARIANTS=(reference compact exact_inputs paged transposed paged_transposed selected_rows selected_paged_transposed)
INDEX=$((SLURM_ARRAY_TASK_ID-1))
REPEAT=$((INDEX/8+1))
# Rotate variant order across repeats to reduce submission-order confounding.
VARIANT=${VARIANTS[(INDEX%8+REPEAT-1)%8]}
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','$VARIANT','probe')"
