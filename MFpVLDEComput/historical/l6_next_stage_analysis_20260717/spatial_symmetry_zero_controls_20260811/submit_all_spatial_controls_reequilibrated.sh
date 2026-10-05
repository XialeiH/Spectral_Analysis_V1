#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT=${1:?run root required}
SETUP=${2:-"${RUN_ROOT}/input/global_bifurcation_setup.mat"}
CODE_ROOT="${RUN_ROOT}/code"
OUTPUT_ROOT="${RUN_ROOT}/reequilibrated_all_controls_20260812"
mkdir -p "${OUTPUT_ROOT}" "${RUN_ROOT}/logs_reequilibrated"

COMMON=(--account=torch_pr_482_general --partition=cpu_short --qos=cpu48 \
  --cpus-per-task=8 --mem=48G --time=06:00:00)
ARRAY_JOB=$(sbatch --parsable "${COMMON[@]}" --array=1-10%10 \
  --job-name=spatial_reeq --output="${RUN_ROOT}/logs_reequilibrated/array_%A_%a.log" \
  --export=ALL,SPATIAL_CONTROL_SETUP="${SETUP}",SPATIAL_CONTROL_OUTPUT="${OUTPUT_ROOT}" \
  --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath('${CODE_ROOT}'); addpath('${CODE_ROOT}/Utils'); run_all_spatial_controls_reequilibrated_task;\"")
PLOT_JOB=$(sbatch --parsable --account=torch_pr_482_general --partition=cpu_short \
  --qos=cpu48 --cpus-per-task=2 --mem=12G --time=01:00:00 \
  --dependency="afterok:${ARRAY_JOB}" --job-name=spatial_reeq_plot \
  --output="${RUN_ROOT}/logs_reequilibrated/plot_%j.log" \
  --export=ALL,SPATIAL_CONTROL_OUTPUT="${OUTPUT_ROOT}" \
  --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath('${CODE_ROOT}'); aggregate_all_spatial_controls_reequilibrated;\"")

MANIFEST="${OUTPUT_ROOT}/submission_manifest.tsv"
printf 'array_job\tplot_job\tsetup\toutput\n%s\t%s\t%s\t%s\n' \
  "${ARRAY_JOB}" "${PLOT_JOB}" "${SETUP}" "${OUTPUT_ROOT}" > "${MANIFEST}"
printf 'array_job=%s\nplot_job=%s\nmanifest=%s\n' "${ARRAY_JOB}" "${PLOT_JOB}" "${MANIFEST}"
