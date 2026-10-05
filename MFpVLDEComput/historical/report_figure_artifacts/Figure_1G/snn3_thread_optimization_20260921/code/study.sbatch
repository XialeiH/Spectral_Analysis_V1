#!/bin/bash
#SBATCH --job-name=snn3_threads
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
module purge
module load matlab/2025b anaconda3/2025.06
source /share/apps/anaconda3/2025.06/etc/profile.d/conda.sh
conda activate py310
python "$OUTPUT_ROOT/code/run_study.py" "$1"
