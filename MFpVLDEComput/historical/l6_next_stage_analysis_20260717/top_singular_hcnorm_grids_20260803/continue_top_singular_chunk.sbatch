#!/bin/bash
#SBATCH --job-name=tshc_next
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M
#SBATCH --time=00:05:00

set -euo pipefail
"${TSHC_CODE_ROOT}/submit_top_singular_chunk.sh" \
  "${TSHC_RUN_ROOT}" "${TSHC_PAIR}" "${TSHC_NEXT_START}"
