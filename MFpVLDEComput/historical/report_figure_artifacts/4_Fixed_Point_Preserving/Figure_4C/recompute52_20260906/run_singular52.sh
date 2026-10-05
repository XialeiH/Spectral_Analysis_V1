#!/bin/bash
set -euo pipefail
module purge
module load matlab/2025b
export GEOM_RUN_ROOT=/scratch/xh2906/librarySCI_runs/figure4c_singular52_20260906
export GEOM_ROOT=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/code
export L6_EQU_BLEND_WEIGHT=0
export MATLAB_PREFDIR=$GEOM_RUN_ROOT/prefs_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$MATLAB_PREFDIR"
matlab -batch "addpath('$GEOM_ROOT'); run('$GEOM_RUN_ROOT/compute_singular52.m')"
