#!/usr/bin/env bash
set -euo pipefail

module purge
module load matlab/2025b

export FIG1F1_REPO_ROOT=/scratch/xh2906/NYU-Vision-2Drive-main
export FIG1F1_HELPER_ROOT=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823/helpers
export FIG1F1_OUTPUT_FILE=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823/results/cg_3dlibrary_angle0_contrast100.mat

matlab -batch "run('/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823/code/run_current_parameter_3d_library_angle0_contrast100.m');"
