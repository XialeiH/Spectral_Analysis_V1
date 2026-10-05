#!/bin/bash
#SBATCH --job-name=fig60123_tuning
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs,cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00

set -euo pipefail
module purge
module load matlab/2025b

export FIG60_SELECTED_INDEX="${SLURM_ARRAY_TASK_ID}"
main_root=/scratch/xh2906/NYU-Vision-2Drive-main
h96_runtime=/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle_extended_corrected_h96_ikp7_l6p5_matchedlib_20260702_010052

matlab -batch "maxNumCompThreads(4); addpath('${main_root}'); addpath('${main_root}/Utils'); addpath('${h96_runtime}','-begin'); addpath('${FIG60_RUN_ROOT}/code','-begin'); recompute_figure6_0_selected_tuning"
