#!/bin/bash
#SBATCH --job-name=fig60_final
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=00:10:00

set -euo pipefail
aggregate_jobs=()
for contrast in 19 42 66 100; do
  points_root="${FIG60_OUTPUT_ROOT}/contrast_$(printf '%03d' "$contrast")/points"
  count=$(find "$points_root" -maxdepth 1 -type f -name 'point_*.mat' | wc -l)
  if [[ "$count" -ne 2601 ]]; then
    echo "Contrast $contrast expected 2601 point files, found $count." >&2
    exit 1
  fi
  job=$(sbatch --parsable \
    --output="${FIG60_OUTPUT_ROOT}/logs/aggregate_c${contrast}_%j.out" \
    --error="${FIG60_OUTPUT_ROOT}/logs/aggregate_c${contrast}_%j.err" \
    --export="ALL,FIG60_CODE_ROOT=${FIG60_CODE_ROOT},FIG60_SETUP_ROOT=${FIG60_SETUP_ROOT},FIG60_OUTPUT_ROOT=${FIG60_OUTPUT_ROOT},FIG60_CONTRAST=${contrast}" \
    "${FIG60_CODE_ROOT}/aggregate_figure6_0_contrast.sbatch")
  aggregate_jobs+=("$job")
  printf 'aggregate\t%d\t%s\n' "$contrast" "$job" >>"${FIG60_RUN_ROOT}/submission_manifest.tsv"
done
dependency=$(IFS=:; echo "${aggregate_jobs[*]}")
plot_job=$(sbatch --parsable --dependency="afterok:${dependency}" \
  --output="${FIG60_OUTPUT_ROOT}/logs/plot_%j.out" \
  --error="${FIG60_OUTPUT_ROOT}/logs/plot_%j.err" \
  --export="ALL,FIG60_CODE_ROOT=${FIG60_CODE_ROOT},FIG60_OUTPUT_ROOT=${FIG60_OUTPUT_ROOT},FIG60_FIGURE_ROOT=${FIG60_OUTPUT_ROOT}/figures" \
  "${FIG60_CODE_ROOT}/plot_figure6_0.sbatch")
printf 'plot\t%s\n' "$plot_job" >>"${FIG60_RUN_ROOT}/submission_manifest.tsv"
printf 'aggregate_jobs=%s plot_job=%s\n' "$dependency" "$plot_job"
