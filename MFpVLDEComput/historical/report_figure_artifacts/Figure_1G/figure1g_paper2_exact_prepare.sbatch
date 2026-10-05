#!/bin/bash
#SBATCH --job-name=f1g_p2_exact_prep
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail
module purge
module load matlab/2025b
mkdir -p "$EXACT_ROOT/inputs" "$EXACT_ROOT/results" "$EXACT_ROOT/logs"
matlab -batch "addpath('$CODE_ROOT'); prepare_figure1g_paper2_exact_inputs('$BASE_ROOT/snn3_bundle.mat','$BASE_ROOT/snn_repo/Data/Paper2_NetworkTuning/Fig1V4/AllMFPixPara_Paper2TuneFig1V4D2.mat','$EXACT_ROOT/inputs')"
