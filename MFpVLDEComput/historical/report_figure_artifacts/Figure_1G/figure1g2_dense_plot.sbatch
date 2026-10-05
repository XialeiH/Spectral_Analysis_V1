#!/bin/bash
#SBATCH --job-name=f1g2_plot
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
module purge; module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); plot_figure1g2_dense_cold('$OUTPUT_ROOT','$OUTPUT_ROOT/Figure 1G.2.pdf')"
