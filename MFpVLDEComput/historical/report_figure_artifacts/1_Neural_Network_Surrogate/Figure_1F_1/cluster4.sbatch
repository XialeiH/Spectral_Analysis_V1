#!/usr/bin/env bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-2
#SBATCH --output=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823/logs/cluster4_%A_%a.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823/logs/cluster4_%A_%a.err

set -euo pipefail
module purge
module load matlab/2025b

RUN_ROOT=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823
if [[ "${SLURM_ARRAY_TASK_ID}" == "1" ]]; then
    SOURCE=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/l6_mechanism_derivative_clamp_20260628_033413/geometry_mat/geometry_sections4_5_h96baseline_L6eqW0p00_contr100_angle_0.00.mat
    matlab -batch "addpath('${RUN_ROOT}/code'); compute_dnn_cluster4('${SOURCE}','${RUN_ROOT}/results/dnn_cluster4.mat');"
else
    export FIG1F1_REPO_ROOT=/scratch/xh2906/NYU-Vision-2Drive-main
    export FIG1F1_HELPER_ROOT=${RUN_ROOT}/helpers
    export FIG1F1_OUTPUT_FILE=${RUN_ROOT}/results/cg_cluster4.mat
    export FIG1F1_CLUSTER4_ONLY=1
    matlab -batch "run('${RUN_ROOT}/code/run_current_parameter_3d_library_angle0_contrast100.m');"
fi
