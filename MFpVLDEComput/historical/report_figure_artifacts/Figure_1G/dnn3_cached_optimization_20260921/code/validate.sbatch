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
VARIANT=$(cat "$OUTPUT_ROOT/selected_variant.txt")
if [[ "$VARIANT" == none ]]; then
    printf 'No validated probe improvement; no validation timing fabricated.\n'
    exit 0
fi
REPEAT=$SLURM_ARRAY_TASK_ID
# Both calls start a new MATLAB process on the same allocation and node.
# Reference timings are contemporaneous controls, never replacements for CG.
if (( REPEAT % 2 )); then
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','$VARIANT','validation')"
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','paged','control')"
else
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','paged','control')"
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','$VARIANT','validation')"
fi
