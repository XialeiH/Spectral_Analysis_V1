#!/bin/bash
#SBATCH --job-name=f1g_snn3_eqplot
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); plot_figure1g_snn3_equilibrium('$OUTPUT_ROOT/snn_repo','$OUTPUT_ROOT/SNN_3x3_equilibrium_angle0.pdf','$OUTPUT_ROOT/SNN_3x3_equilibrium_angle0.png')"
