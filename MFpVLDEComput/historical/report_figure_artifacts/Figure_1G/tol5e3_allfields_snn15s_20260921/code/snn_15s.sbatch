#!/bin/bash
#SBATCH --job-name=f1b_snn15s
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-5%5
set -euo pipefail
module purge
module load matlab/2025b
TRIAL=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
TRIAL_ROOT="$OUTPUT_ROOT/snn_trial_$TRIAL"
DATA_ROOT="$TRIAL_ROOT/Data/Paper2_NetworkTuning/Fig1V4"
SAVE_ROOT="$DATA_ROOT/Paper2NWSimulationData"
mkdir -p "$DATA_ROOT/Paper2PlotingData" "$DATA_ROOT/GlobConv" "$SAVE_ROOT" "$TRIAL_ROOT/Paper2Figs"
ln -s "$BASE_ROOT/snn_repo/Utils" "$TRIAL_ROOT/Utils"
ln -s "$EXACT_ROOT/inputs/AllMFPixPara_Paper2TuneFig1V4D2_torch_staged.mat" "$DATA_ROOT/AllMFPixPara_Paper2TuneFig1V4D2.mat"
ln -s "$EXACT_ROOT/inputs/LargeNWFixIni.mat" "$TRIAL_ROOT/LargeNWFixIni.mat"
ln -s "$EXACT_ROOT/inputs/Paper2DriveNW_Conn.mat" "$SAVE_ROOT/Paper2DriveNW_Conn.mat"
cp "$OUTPUT_ROOT/code/Paper2_Fig7Comp_NW_LDE.m" "$TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m"
echo "$DRIVER_SHA  $TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
cd "$TRIAL_ROOT"
START=$(date +%s.%N)
matlab -batch "run('$TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m'); assert(T==15000 && dt==0.1 && TPar==15 && WinSize==9000); assert(~BlowUp); assert(TSec==15 && TimeN==10000);"
END=$(date +%s.%N)
ELAPSED=$(awk -v start="$START" -v end="$END" 'BEGIN {printf "%.9f", end-start}')
SEGMENT_COUNT=$(find "$SAVE_ROOT" -maxdepth 1 -type f -name 'DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_*s_NewSmear.mat' | wc -l)
test "$SEGMENT_COUNT" -eq 15
test -s "$SAVE_ROOT/NWSimulationPix_0.0deg_NewSmear.mat"
echo "$DRIVER_SHA  $TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
OUT="$OUTPUT_ROOT/results/snn_paper2_15s_03HC_repeat_${TRIAL}.tsv"
printf 'FieldHC\tRepeat\tMethod\tSeconds\tSimulationMs\tdtMs\tUpdates\tCPUs\tProtocol\n' > "$OUT"
printf '3\t%d\tSNN\t%s\t15000\t0.1\t150000\t16\tPaper2OriginalDuration15s\n' "$SLURM_ARRAY_TASK_ID" "$ELAPSED" >> "$OUT"
