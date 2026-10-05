#!/bin/bash
#SBATCH --job-name=f1b_dnn_select
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
module purge
module load matlab/2025b
matlab -batch "addpath('$OUTPUT_ROOT/code'); analyze_thread_probe('$OUTPUT_ROOT')"
