#!/usr/bin/env bash
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823/logs/cg_cluster_compare_%j.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823/logs/cg_cluster_compare_%j.err

set -euo pipefail
module purge
module load matlab/2025b

RUN_ROOT=/scratch/xh2906/librarySCI_runs/figure_1f1_current_dnn_vs_3dlibrary_20260823
export FIG1F1_REPO_ROOT=/scratch/xh2906/NYU-Vision-2Drive-main
export FIG1F1_HELPER_ROOT=${RUN_ROOT}/helpers
export FIG1F1_OUTPUT_FILE=${RUN_ROOT}/results/cg_cluster_compare_6_8_10.mat
export FIG1F1_CLUSTER_COMPARE=1

matlab -batch "run('${RUN_ROOT}/code/run_current_parameter_3d_library_angle0_contrast100.m');"
