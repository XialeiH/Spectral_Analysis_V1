#!/bin/bash
#SBATCH --job-name=l6i_grid_next
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M
#SBATCH --time=00:05:00

set -euo pipefail
"${FULL_GRID_CODE_ROOT}/submit_full_grid_chunk.sh" \
  "${FULL_GRID_RUN_ROOT}" "${FULL_GRID_NEXT_START}"
