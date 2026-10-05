#!/bin/bash
#SBATCH --job-name=boxed_true_stab
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
module purge
module load matlab/2025b

main_root=/scratch/xh2906/NYU-Vision-2Drive-main
h96_runtime=/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle
code_root=/scratch/xh2906/librarySCI_runs/l6_i_offset_canonicalHC_full_20260727_175000/code

matlab -batch "addpath('${main_root}'); addpath('${main_root}/Utils'); addpath('${h96_runtime}','-begin'); addpath('${code_root}','-begin'); evaluate_boxed_true_stability"
