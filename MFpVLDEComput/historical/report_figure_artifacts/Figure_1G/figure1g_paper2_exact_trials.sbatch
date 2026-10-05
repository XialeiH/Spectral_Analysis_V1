#!/bin/bash
#SBATCH --job-name=f1g_p2_exact
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cpu_short
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-5%5
set -euo pipefail
module purge
module load matlab/2025b

TRIAL=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
TRIAL_ROOT="$EXACT_ROOT/original_trial_$TRIAL"
DATA_ROOT="$TRIAL_ROOT/Data/Paper2_NetworkTuning/Fig1V4"
SAVE_ROOT="$DATA_ROOT/Paper2NWSimulationData"
mkdir -p "$DATA_ROOT/Paper2PlotingData" "$DATA_ROOT/GlobConv" "$SAVE_ROOT" "$TRIAL_ROOT/Paper2Figs"
ln -s "$BASE_ROOT/snn_repo/Utils" "$TRIAL_ROOT/Utils"
ln -s "$BASE_ROOT/snn_repo/Data/Paper2_NetworkTuning/Fig1V4/AllMFPixPara_Paper2TuneFig1V4D2.mat" "$DATA_ROOT/AllMFPixPara_Paper2TuneFig1V4D2.mat"
ln -s "$EXACT_ROOT/inputs/LargeNWFixIni.mat" "$TRIAL_ROOT/LargeNWFixIni.mat"
ln -s "$EXACT_ROOT/inputs/Paper2DriveNW_Conn.mat" "$SAVE_ROOT/Paper2DriveNW_Conn.mat"
cp "$BASE_ROOT/snn_repo/Paper2_Fig7Comp_NW_LDE.m" "$TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m"

cd "$TRIAL_ROOT"
START=$(date +%s.%N)
matlab -batch "run('$TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m')"
END=$(date +%s.%N)
ELAPSED=$(awk -v start="$START" -v end="$END" 'BEGIN {printf "%.6f",end-start}')
OUT="$EXACT_ROOT/results/snn_paper2_exact_03HC_repeat_${TRIAL}.tsv"
printf 'FieldHC\tRepeat\tMethod\tSeconds\tSimulatedMs\tDtMs\tSteps\tCPUs\tPipeline\n' > "$OUT"
printf '3\t%d\tSNN\t%s\t10000\t0.1\t100000\t16\tPaper2OriginalExact\n' "$SLURM_ARRAY_TASK_ID" "$ELAPSED" >> "$OUT"
