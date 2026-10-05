#!/bin/bash
#SBATCH --job-name=f1b_dnn_select
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:10:00
set -euo pipefail
python3 "$OUTPUT_ROOT/code/analyze_impl.py" "$OUTPUT_ROOT" --select
