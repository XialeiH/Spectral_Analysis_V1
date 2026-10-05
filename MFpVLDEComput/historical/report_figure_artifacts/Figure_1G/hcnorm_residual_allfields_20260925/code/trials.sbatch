#!/bin/bash
#SBATCH --job-name=f1b_hcnorm_trials
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-610%16
#SBATCH --output=/scratch/xh2906/librarySCI_runs/figure1b_hcnorm_residual_allfields_20260925/logs/%x_%A_%a.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/figure1b_hcnorm_residual_allfields_20260925/logs/%x_%A_%a.err
set -euo pipefail
module purge
module load matlab/2025b

INDEX=$((SLURM_ARRAY_TASK_ID-1))
if (( INDEX < 290 )); then
    METHOD=CG
    if (( INDEX < 200 )); then
        FIELD=$((3+INDEX/100)); REPEAT=$((INDEX%100+1))
    elif (( INDEX < 280 )); then
        FIELDS=(6 8 10 20); J=$((INDEX-200))
        FIELD=${FIELDS[J/20]}; REPEAT=$((J%20+1))
    else
        FIELDS=(30 40); J=$((INDEX-280))
        FIELD=${FIELDS[J/5]}; REPEAT=$((J%5+1))
    fi
else
    METHOD=DNN
    J=$((INDEX-290))
    if (( J < 200 )); then
        FIELD=$((3+J/100)); REPEAT=$((J%100+1))
    else
        FIELDS=(6 8 10 20 30 40); J=$((J-200))
        FIELD=${FIELDS[J/20]}; REPEAT=$((J%20+1))
    fi
fi

if [[ "$METHOD" == CG ]]; then
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_figure1g2_tol5e3('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','CG')"
elif (( FIELD == 3 )); then
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$OUTPUT_ROOT','$REPEAT','paged_combined','results')"
else
    matlab -batch "addpath('$OUTPUT_ROOT/code'); run_allfields('$SOURCE_ROOT','$OUTPUT_ROOT','$FIELD','$REPEAT','paged_combined','results')"
fi
