#!/bin/bash
#SBATCH --job-name=pair_hc_ode_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
points_root="${PAIR_GRID_ODE_POINTS_ROOT}"
combined="${PAIR_GRID_ODE_COMBINED}"
expected=$(wc -l < "${PAIR_GRID_ODE_TASK_LIST}")
count=$(find "$points_root" -maxdepth 1 -type f -name 'ode_point_*.tsv' | wc -l)
if [[ "$count" -ne "$expected" ]]; then
  echo "Expected $expected ODE point files, found $count." >&2
  exit 1
fi
printf 'taskId\tbetaXIndex\tbetaYIndex\tbeta6\tbetaE\tconverged\tfixedPointResidual\todeIntrinsicTime\todeAcceptedSteps\tHCnorm\tminimumRateHz\tmaximumRateHz\ttermination\n' > "$combined"
find "$points_root" -maxdepth 1 -type f -name 'ode_point_*.tsv' -print0 | sort -z | xargs -0 cat >> "$combined"
module purge
module load matlab/2025b
matlab -batch "addpath('${PAIR_GRID_CODE_ROOT}','-begin'); merge_pair_grid_ode_fallback"
