#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT=${1:?run root required}
SETUP=${2:?setup file required}
PREFLIGHT_JOB=${3:-}
CODE_ROOT="${RUN_ROOT}/code"
SUPPORT_ROOT="${RUN_ROOT}/support"
RUNTIME_ROOT="/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle_extended_corrected_h96_ikp7_l6p5_matchedlib_20260702_010052"
OUTPUT_ROOT="${RUN_ROOT}/results"
mkdir -p "${OUTPUT_ROOT}" "${RUN_ROOT}/logs"

MATLAB_PREFIX="addpath('${RUNTIME_ROOT}','-begin'); addpath('${SUPPORT_ROOT}','-begin'); addpath('${CODE_ROOT}','-begin');"
MATLAB_MODULE="matlab/2025b"
COMMON=(--account=torch_pr_482_general --partition=cpu_short --qos=cpu48 \
  --cpus-per-task=8 --mem=64G --time=06:00:00)

DEPENDENCY_ARGS=()
if [[ -n "${PREFLIGHT_JOB}" ]]; then
  DEPENDENCY_ARGS=(--dependency="afterok:${PREFLIGHT_JOB}")
fi

EQUILIBRIUM_JOB=$(sbatch --parsable --account=torch_pr_482_general \
  --partition=cpu_short --qos=cpu48 --cpus-per-task=4 --mem=24G --time=02:00:00 \
  "${DEPENDENCY_ARGS[@]}" --array=1-20%20 --job-name=gke_equilibrium \
  --output="${RUN_ROOT}/logs/equilibrium_%A_%a.log" \
  --export=ALL,GKE_SETUP="${SETUP}",GKE_OUTPUT="${OUTPUT_ROOT}" \
  --wrap="module purge; module load ${MATLAB_MODULE}; matlab -batch \"${MATLAB_PREFIX} gke_run_equilibrium_condition;\"")
EQUILIBRIUM_PLOT_JOB=$(sbatch --parsable --account=torch_pr_482_general \
  --partition=cpu_short --qos=cpu48 --cpus-per-task=2 --mem=16G --time=01:00:00 \
  --dependency="afterok:${EQUILIBRIUM_JOB}" --job-name=gke_equilibrium_plots \
  --output="${RUN_ROOT}/logs/equilibrium_plots_%j.log" \
  --export=ALL,GKE_OUTPUT="${OUTPUT_ROOT}" \
  --wrap="module purge; module load ${MATLAB_MODULE}; matlab -batch \"${MATLAB_PREFIX} gke_plot_equilibrium_maps('${OUTPUT_ROOT}');\"")

BASELINE_JOB=$(sbatch --parsable "${COMMON[@]}" \
  --dependency="afterok:${EQUILIBRIUM_JOB}" --job-name=gke_baseline_full \
  --output="${RUN_ROOT}/logs/baseline_full_%j.log" \
  --wrap="module purge; module load ${MATLAB_MODULE}; matlab -batch \"${MATLAB_PREFIX} gke_run_condition(1,'${SETUP}','${OUTPUT_ROOT}');\"")
ARRAY_JOB=$(sbatch --parsable "${COMMON[@]}" --array=2-20%19 \
  --dependency="afterok:${BASELINE_JOB}" \
  --job-name=gke_controls --output="${RUN_ROOT}/logs/condition_%A_%a.log" \
  --export=ALL,GKE_SETUP="${SETUP}",GKE_OUTPUT="${OUTPUT_ROOT}" \
  --wrap="module purge; module load ${MATLAB_MODULE}; matlab -batch \"${MATLAB_PREFIX} gke_run_condition;\"")
PLOT_JOB=$(sbatch --parsable --account=torch_pr_482_general --partition=cpu_short \
  --qos=cpu48 --cpus-per-task=2 --mem=32G --time=03:00:00 \
  --dependency="afterok:${BASELINE_JOB}:${ARRAY_JOB}" --job-name=gke_aggregate \
  --output="${RUN_ROOT}/logs/aggregate_%j.log" \
  --export=ALL,GKE_OUTPUT="${OUTPUT_ROOT}" \
  --wrap="module purge; module load ${MATLAB_MODULE}; matlab -batch \"${MATLAB_PREFIX} gke_aggregate;\"")

MANIFEST="${OUTPUT_ROOT}/submission_manifest.tsv"
printf 'equilibrium_job\tequilibrium_plot_job\tbaseline_job\tarray_job\taggregate_job\tsetup\toutput\truntime\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${EQUILIBRIUM_JOB}" "${EQUILIBRIUM_PLOT_JOB}" "${BASELINE_JOB}" \
  "${ARRAY_JOB}" "${PLOT_JOB}" "${SETUP}" "${OUTPUT_ROOT}" "${RUNTIME_ROOT}" > "${MANIFEST}"
printf 'equilibrium_job=%s\nequilibrium_plot_job=%s\nbaseline_job=%s\narray_job=%s\naggregate_job=%s\nmanifest=%s\n' \
  "${EQUILIBRIUM_JOB}" "${EQUILIBRIUM_PLOT_JOB}" "${BASELINE_JOB}" \
  "${ARRAY_JOB}" "${PLOT_JOB}" "${MANIFEST}"
