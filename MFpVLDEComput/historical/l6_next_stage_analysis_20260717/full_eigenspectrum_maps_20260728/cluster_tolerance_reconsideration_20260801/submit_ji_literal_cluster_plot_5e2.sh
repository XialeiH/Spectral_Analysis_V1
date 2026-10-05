#!/bin/bash
#SBATCH --job-name=ji-cluster-plot
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/xh2906/librarySCI_runs/h96_full_eigenspectrum_maps_20260728_160249/ji_cluster_tol5e2_20260805/logs/plot_%j.out
#SBATCH --error=/scratch/xh2906/librarySCI_runs/h96_full_eigenspectrum_maps_20260728_160249/ji_cluster_tol5e2_20260805/logs/plot_%j.err

set -euo pipefail

RUN_ROOT=/scratch/xh2906/librarySCI_runs/h96_full_eigenspectrum_maps_20260728_160249
OUTPUT_ROOT="$RUN_ROOT/ji_cluster_tol5e2_20260805"

module purge
module load matlab/2025b

matlab -batch "addpath('$RUN_ROOT/code'); plot_ji_literal_cluster_tolerance_5e2('$OUTPUT_ROOT/data', '$OUTPUT_ROOT/figures');"
