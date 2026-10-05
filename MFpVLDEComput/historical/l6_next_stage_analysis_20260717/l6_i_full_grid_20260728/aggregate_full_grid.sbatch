#!/bin/bash
#SBATCH --job-name=l6i_hc_grid_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
points_root="${FULL_GRID_OUTPUT_ROOT}/points"
combined="${FULL_GRID_OUTPUT_ROOT}/all_beta6_betaI_grid_points.tsv"
count=$(find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' | wc -l)
if [[ "$count" -ne 40401 ]]; then
  echo "Expected 40401 point files, found $count." >&2
  exit 1
fi
printf 'taskId\tbeta6Index\tbetaIIndex\tbeta6\tbetaI\tconverged\tfixedPointResidual\titerations\tHCnorm\tminimumRateHz\tmaximumRateHz\n' > "$combined"
find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' -print0 | sort -z | xargs -0 cat >> "$combined"
module purge
module load matlab/2025b
matlab -batch "addpath('${FULL_GRID_CODE_ROOT}','-begin'); aggregate_full_grid"
