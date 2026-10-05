#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis/matlab-inserting_into_CG_model"
CODE_ROOT="${PROJECT_ROOT}/l6_next_stage_analysis_20260717/gaussian_kernel_experiments_20260813"
CONTROL_ROOT="${PROJECT_ROOT}/l6_next_stage_analysis_20260717/spatial_symmetry_zero_controls_20260811"
RUNTIME_ROOT="${PROJECT_ROOT}/pixel_size_scaling_20260813/runtime_override"
UTILS_ROOT="${PROJECT_ROOT}/Complete_Code_for_Paper3/NYU-Vision-2Drive-main/Utils"
SETUP=${1:?setup MAT file required}
OUTPUT=${2:?output root required}
PARALLEL=${3:-4}
MATLAB="/Applications/MATLAB_R2025b.app/bin/matlab"
LOG_ROOT="${OUTPUT}/equilibrium_only/logs"
mkdir -p "${LOG_ROOT}"

for ((batchStart=1; batchStart<=20; batchStart+=PARALLEL)); do
  pids=()
  for ((offset=0; offset<PARALLEL; offset++)); do
    task=$((batchStart+offset))
    if ((task>20)); then break; fi
    "${MATLAB}" -batch "addpath('${UTILS_ROOT}','-begin'); addpath('${RUNTIME_ROOT}','-begin'); addpath('${CONTROL_ROOT}','-begin'); addpath('${CODE_ROOT}','-begin'); gke_run_equilibrium_condition(${task},'${SETUP}','${OUTPUT}');" \
      > "${LOG_ROOT}/condition_${task}.log" 2>&1 &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do wait "${pid}"; done
done
"${MATLAB}" -batch "addpath('${CODE_ROOT}','-begin'); gke_plot_equilibrium_maps('${OUTPUT}');"
