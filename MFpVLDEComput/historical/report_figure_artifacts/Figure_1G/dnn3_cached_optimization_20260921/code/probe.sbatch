#!/bin/bash
#SBATCH --job-name=f1b_dnn_impl
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
VARIANTS=(paged paged_quad paged_padding paged_aggregate paged_combined)
INDEX=$((SLURM_ARRAY_TASK_ID-1))
REPEAT=$((INDEX/5+1))
# Rotate variant order across repeats to reduce submission-order confounding.
VARIANT=${VARIANTS[(INDEX%5+REPEAT-1)%5]}
matlab -batch "addpath('$OUTPUT_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','$VARIANT','probe')"
