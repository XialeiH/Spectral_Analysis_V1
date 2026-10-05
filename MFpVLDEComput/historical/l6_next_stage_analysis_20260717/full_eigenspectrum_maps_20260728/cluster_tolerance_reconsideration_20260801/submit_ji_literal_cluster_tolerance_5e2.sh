#!/bin/bash
#SBATCH --job-name=ji-cluster-5e2
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --array=1-2%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_full_eigenspectrum_maps_20260728_160249/ji_cluster_tol5e2_20260805/logs/compute_%A_%a.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/h96_full_eigenspectrum_maps_20260728_160249/ji_cluster_tol5e2_20260805/logs/compute_%A_%a.err

set -euo pipefail

RUN_ROOT=/scratch/xh2906/librarySCI_runs/h96_full_eigenspectrum_maps_20260728_160249
OUTPUT_ROOT="$RUN_ROOT/ji_cluster_tol5e2_20260805"

module purge
module load matlab/2025b

matlab -batch "addpath('$RUN_ROOT/code'); run_ji_literal_cluster_tolerance_task(str2double(getenv('SLURM_ARRAY_TASK_ID')), '$RUN_ROOT/results/data', '$OUTPUT_ROOT/data', 5e-2);"
