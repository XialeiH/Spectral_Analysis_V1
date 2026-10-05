#!/bin/bash
set -euo pipefail

BASE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610
RUN_ROOT=${BASE}/perturbation_7_20260807
CODE_ROOT=${RUN_ROOT}/code
UTIL_ROOT=${BASE}/code/Utils
SETUP_FILE=${CODE_ROOT}/global_bifurcation_setup.mat
TRAJECTORY_ROOT=${RUN_ROOT}/dense_finite_time_trajectories
TASK_ROOT=${RUN_ROOT}/dense_finite_time_alignment_tasks
CURVE_ROOT=${RUN_ROOT}/dense_finite_time_alignment_curves
OUTPUT_ROOT=${RUN_ROOT}/all_figures_dense_finite_time_baseline_alignment_pdfs
LOG_ROOT=${RUN_ROOT}/logs_dense_finite_time

rm -rf "${TRAJECTORY_ROOT}" "${TASK_ROOT}" "${CURVE_ROOT}" "${OUTPUT_ROOT}" "${LOG_ROOT}"
mkdir -p "${TRAJECTORY_ROOT}" "${TASK_ROOT}" "${CURVE_ROOT}" "${OUTPUT_ROOT}" "${LOG_ROOT}"

trajectory_job=$(sbatch --parsable \
    --job-name=fig7denseprep --account=torch_pr_482_general \
    --partition=cpu_short --qos=cpu48 --array=1-12%12 \
    --cpus-per-task=2 --mem=12G --time=01:00:00 \
    --output="${LOG_ROOT}/prep_%A_%a.out" \
    --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); prepare_dense_finite_time_trajectory_task(str2double(getenv('SLURM_ARRAY_TASK_ID')),'${SETUP_FILE}','${RUN_ROOT}','${TRAJECTORY_ROOT}');\"")

dense_job=$(sbatch --parsable --dependency=afterok:${trajectory_job} \
    --job-name=fig7denseft --account=torch_pr_482_general \
    --partition=cpu_short --qos=cpu48 --array=1-1004%120 \
    --cpus-per-task=2 --mem=16G --time=01:00:00 \
    --output="${LOG_ROOT}/dense_%A_%a.out" \
    --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); compute_dense_finite_time_alignment_block_task(str2double(getenv('SLURM_ARRAY_TASK_ID')),'${SETUP_FILE}','${TRAJECTORY_ROOT}','${TASK_ROOT}');\"")

aggregate_job=$(sbatch --parsable --dependency=afterok:${dense_job} \
    --job-name=fig7denseagg --account=torch_pr_482_general \
    --partition=cpu_short --qos=cpu48 --cpus-per-task=1 --mem=4G --time=00:20:00 \
    --output="${LOG_ROOT}/aggregate_%j.out" \
    --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath('${CODE_ROOT}','-begin'); aggregate_dense_finite_time_alignment('${TASK_ROOT}','${TRAJECTORY_ROOT}','${CURVE_ROOT}');\"")

plot_job=$(sbatch --parsable --dependency=afterok:${aggregate_job} \
    --job-name=fig7denseplot --account=torch_pr_482_general \
    --partition=cpu_short --qos=cpu48 --cpus-per-task=4 --mem=16G --time=01:00:00 \
    --output="${LOG_ROOT}/plot_%j.out" \
    --wrap="module purge; module load matlab/2025b; matlab -batch \"addpath('${CODE_ROOT}','-begin'); addpath('${UTIL_ROOT}','-begin'); replot_all_perturbation_figures_with_finite_time('${RUN_ROOT}','${RUN_ROOT}/all_figures_finite_time_rows','${OUTPUT_ROOT}','${CURVE_ROOT}');\"")

printf 'trajectory_job\t%s\n' "${trajectory_job}"
printf 'dense_job\t%s\n' "${dense_job}"
printf 'aggregate_job\t%s\n' "${aggregate_job}"
printf 'plot_job\t%s\n' "${plot_job}"
