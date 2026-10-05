#!/bin/bash
#SBATCH --job-name=norm_e_lr
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-16%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/normalized_e_heatmap_left_right_20260807/logs/%A_%a.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/normalized_e_heatmap_left_right_20260807/logs/%A_%a.err

set -euo pipefail
module purge
module load matlab/2025b

RUN_ROOT=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/normalized_e_heatmap_left_right_20260807
SOURCE_ROOT=/scratch/xh2906/librarySCI_runs/h96_default_full_geometry_dedupe_array_localhelpers_20260607_084610/l6_mechanism_derivative_clamp_20260628_033413

matlab -batch "addpath('${RUN_ROOT}/code'); compute_normalized_e_heatmap_left_right_task(${SLURM_ARRAY_TASK_ID}, '${SOURCE_ROOT}', '${RUN_ROOT}/code/left_right_mode_similarity_all_conditions.tsv', '${RUN_ROOT}/results');"
