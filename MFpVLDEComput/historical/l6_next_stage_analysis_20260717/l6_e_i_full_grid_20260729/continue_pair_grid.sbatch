#!/bin/bash
#SBATCH --job-name=pair_grid_next
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M
#SBATCH --time=00:05:00

set -euo pipefail
"${PAIR_GRID_CODE_ROOT}/submit_pair_grid_chunk.sh" \
  "${PAIR_GRID_RUN_ROOT}" "${PAIR_GRID_PAIR}" \
  "${PAIR_GRID_SETUP_FILE}" "${PAIR_GRID_NEXT_START}"
