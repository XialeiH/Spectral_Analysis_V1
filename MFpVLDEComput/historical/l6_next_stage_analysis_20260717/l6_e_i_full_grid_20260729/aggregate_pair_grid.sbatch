#!/bin/bash
#SBATCH --job-name=pair_hc_grid_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
points_root="${PAIR_GRID_OUTPUT_ROOT}/points"
combined="${PAIR_GRID_OUTPUT_ROOT}/all_pair_grid_points.tsv"
count=$(find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' | wc -l)
if [[ "$count" -ne 6561 ]]; then
  echo "Expected 6561 point files, found $count." >&2
  exit 1
fi
printf 'taskId\tbetaXIndex\tbetaYIndex\tbetaX\tbetaY\tbeta6\tbetaE\tbetaI\tconverged\tfixedPointResidual\titerations\tHCnorm\tminimumRateHz\tmaximumRateHz\n' > "$combined"
find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' -print0 | sort -z | xargs -0 cat >> "$combined"
module purge
module load matlab/2025b
matlab -batch "addpath('${PAIR_GRID_CODE_ROOT}','-begin'); aggregate_pair_grid"

# Optionally launch a second pair scan after this grid has fully aggregated.
# Sequencing is required because Torch counts array elements against the
# per-user submitted-job limit.
run_root=$(dirname "${PAIR_GRID_OUTPUT_ROOT}")
next_config="${run_root}/next_pair_run.conf"
if [[ -f "$next_config" ]]; then
  # shellcheck disable=SC1090
  source "$next_config"
  "${NEXT_PAIR_RUN_ROOT}/code/submit_pair_grid.sh" \
    "${NEXT_PAIR_RUN_ROOT}" "${NEXT_PAIR_TYPE}" "${NEXT_PAIR_SETUP_FILE}"
fi
