#!/bin/bash
#SBATCH --job-name=fig60_next
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M
#SBATCH --time=00:05:00

set -euo pipefail
"${FIG60_CODE_ROOT}/submit_figure6_0_chunk.sh" "${FIG60_RUN_ROOT}" "${FIG60_NEXT_START}"
