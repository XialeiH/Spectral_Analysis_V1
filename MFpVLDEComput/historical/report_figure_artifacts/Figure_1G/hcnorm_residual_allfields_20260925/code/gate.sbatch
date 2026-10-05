#!/bin/bash
#SBATCH --job-name=f1b_hcnorm_gate
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=1-4%4
#SBATCH --output=/scratch/xh2906/librarySCI_runs/figure1b_hcnorm_residual_allfields_20260925/logs/%x_%A_%a.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/figure1b_hcnorm_residual_allfields_20260925/logs/%x_%A_%a.err
set -euo pipefail
module purge
module load matlab/2025b

case "$SLURM_ARRAY_TASK_ID" in
    1) FIELD=3; METHOD=CG ;;
    2) FIELD=3; METHOD=DNN ;;
    3) FIELD=4; METHOD=CG ;;
    4) FIELD=4; METHOD=DNN ;;
esac

GATE_ROOT="$OUTPUT_ROOT/gate"
if [[ "$METHOD" == CG ]]; then
    matlab -batch "addpath('$GATE_ROOT/code'); run_figure1g2_tol5e3('$SOURCE_ROOT','$GATE_ROOT','$FIELD','1','CG')"
elif (( FIELD == 3 )); then
    matlab -batch "addpath('$GATE_ROOT/code'); run_impl_probe('$SOURCE_ROOT','$GATE_ROOT','1','paged_combined','results')"
else
    matlab -batch "addpath('$GATE_ROOT/code'); run_allfields('$SOURCE_ROOT','$GATE_ROOT','$FIELD','1','paged_combined','results')"
fi
