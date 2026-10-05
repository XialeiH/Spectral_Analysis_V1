#!/bin/bash
#SBATCH --job-name=tshc_agg
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail
points_root="${TSHC_OUTPUT_ROOT}/points"
combined="${TSHC_OUTPUT_ROOT}/all_top_singular_grid_points.tsv"
count=$(find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' | wc -l)
if [[ "$count" -ne 2601 ]]; then
  echo "Expected 2601 point files, found $count." >&2
  exit 1
fi
printf 'taskId\tbetaXIndex\tbetaYIndex\tbetaX\tbetaY\tbeta6\tbetaE\tbetaI\tfixedPointAccepted\tfixedPointResidual\titerations\tmodeHCnorm\tsigma1\tsigma2\tsingularGapRelative\tmodeOverlap\tsingularResidual\tsvdsSeconds\tminimumRateHz\tmaximumRateHz\n' > "$combined"
find "$points_root" -maxdepth 1 -type f -name 'point_*.tsv' -print0 | sort -z | xargs -0 cat >> "$combined"
module purge
module load matlab/2025b
matlab -batch "addpath('${TSHC_CODE_ROOT}','-begin'); aggregate_top_singular_pair"

case "${TSHC_PAIR}" in
  l6_inhibition) next_pair=l6_excitation ;;
  l6_excitation) next_pair=excitation_inhibition ;;
  excitation_inhibition) next_pair= ;;
esac
if [[ -n "$next_pair" ]]; then
  "${TSHC_CODE_ROOT}/submit_top_singular_pair.sh" "${TSHC_RUN_ROOT}" "$next_pair"
else
  sbatch --parsable \
    --output="${TSHC_RUN_ROOT}/final_plot_%j.out" \
    --error="${TSHC_RUN_ROOT}/final_plot_%j.err" \
    --export="ALL,TSHC_RUN_ROOT=${TSHC_RUN_ROOT},TSHC_CODE_ROOT=${TSHC_CODE_ROOT}" \
    "${TSHC_CODE_ROOT}/plot_three_top_singular_hcnorm.sbatch" | \
    awk '{print "final_plot\t" $1}' >> "${TSHC_RUN_ROOT}/submission_manifest.tsv"
fi
