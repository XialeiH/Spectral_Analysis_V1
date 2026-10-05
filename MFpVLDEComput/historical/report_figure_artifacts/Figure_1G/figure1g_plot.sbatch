#!/bin/bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --job-name=fig1g_plot
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
module purge
module load matlab/2025b
OLD_ROOT=/scratch/xh2906/librarySCI_runs/figure1g_cg_vs_h96dnn_100trials_cpu8_20260823
matlab -batch "addpath('$OLD_ROOT/code'); aggregate_figure1g_results('$OLD_ROOT'); addpath('$RUN_ROOT/code'); aggregate_figure1g_scaling('$RUN_ROOT'); plot_figure1g_scaling(fullfile('$OLD_ROOT','figure1g_timings.tsv'),fullfile('$RUN_ROOT','field_scaling_timings.tsv'),fullfile('$RUN_ROOT','Figure 1G.pdf'))"
