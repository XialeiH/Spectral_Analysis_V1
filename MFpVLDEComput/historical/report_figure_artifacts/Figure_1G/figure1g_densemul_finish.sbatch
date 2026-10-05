#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --job-name=fig1g_dplot
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('$RUN_ROOT/code'); aggregate_figure1g_densemul('$RUN_ROOT'); plot_figure1g_densemul('$RUN_ROOT','$RUN_ROOT/Figure 1G.2.pdf')"
