#!/bin/bash
#SBATCH --job-name=f1g_snn3_bundle
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); build_figure1g_snn3_bundle('$OUTPUT_ROOT/snn_repo','$OUTPUT_ROOT/snn3_bundle.mat')"
