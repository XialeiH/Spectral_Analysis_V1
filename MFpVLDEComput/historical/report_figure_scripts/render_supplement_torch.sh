#!/bin/bash
#SBATCH --job-name=supplement-render
#SBATCH --account=torch_pr_482_general
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
set -euo pipefail
module load matlab/2025b
export SUPPLEMENT_PROJECT_ROOT=/scratch/xh2906/librarySCI_runs/supplement_S1_S2_20260909
export FONTCONFIG_FILE="$SUPPLEMENT_PROJECT_ROOT/supplement_fonts.conf"
cd "$SUPPLEMENT_PROJECT_ROOT"
figure="${1:-S2}"
matlab -batch "addpath(pwd); plot_supplement_S1_S2('$figure');"
