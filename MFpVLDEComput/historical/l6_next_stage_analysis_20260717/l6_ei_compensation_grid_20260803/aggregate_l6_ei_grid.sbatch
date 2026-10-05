#!/bin/bash
#SBATCH --job-name=l6_ei_hc_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
points_root="${L6EI_OUTPUT_ROOT}/points"
combined="${L6EI_OUTPUT_ROOT}/all_beta6_betaEI_grid_points.tsv"
expected_total="${L6EI_TOTAL_GRID_COUNT:-2911}"
count=$(find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' | wc -l)
if [[ "$count" -ne "$expected_total" ]]; then
  echo "Expected $expected_total point files, found $count." >&2
  exit 1
fi
printf 'taskId\tbeta6Index\tbetaEIIndex\tbeta6\tbetaEI\tgain6\tgainEI\tstrictConverged\taccepted\tfixedPointResidual\titerations\tHCnorm\tminimumRateHz\tmaximumRateHz\n' > "$combined"
find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' -print0 | sort -z | xargs -0 cat >> "$combined"
module purge
module load matlab/2025b
matlab -batch "addpath('${L6EI_CODE_ROOT}','-begin'); aggregate_l6_ei_grid"
