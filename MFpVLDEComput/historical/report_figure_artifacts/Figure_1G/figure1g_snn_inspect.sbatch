#!/bin/bash
#SBATCH --job-name=f1g_snn_inspect
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); inspect_figure1g_snn_source('$OUTPUT_ROOT/snn_repo')"
