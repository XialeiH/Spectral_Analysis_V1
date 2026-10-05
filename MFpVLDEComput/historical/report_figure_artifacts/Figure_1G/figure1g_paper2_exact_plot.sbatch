#!/bin/bash
#SBATCH --job-name=f1g_p2_exact_plot
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
module purge
module load matlab/2025b
cp "$EXACT_ROOT"/results/snn_paper2_exact_03HC_repeat_*.tsv "$BASE_ROOT/results/"
matlab -batch "addpath('$BASE_ROOT/code'); plot_figure1g2_dense_cold_snn('$BASE_ROOT','$BASE_ROOT/Figure 1G.2.pdf')"
