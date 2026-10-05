#!/bin/bash
#SBATCH --job-name=l6_ei_grid_next
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M
#SBATCH --time=00:05:00

set -euo pipefail
"${L6EI_CODE_ROOT}/submit_l6_ei_grid_chunk.sh" \
  "${L6EI_RUN_ROOT}" "${L6EI_NEXT_START}"
