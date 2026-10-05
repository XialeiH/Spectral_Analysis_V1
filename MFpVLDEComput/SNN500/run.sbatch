#!/bin/bash
#SBATCH --job-name=snn_bounded
#SBATCH --account=torch_pr_482_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --cpus-per-task=16
set -euo pipefail
module purge
module load matlab/2025b
ROOT=/scratch/xh2906/librarySCI_runs/snn3_duration400s_500s_bounded_20260923
BASE=/scratch/xh2906/librarySCI_runs/figure1g2_dense_cold_equilibrium_20260829
INPUT="$BASE/snn_paper2_original_exact_20260901_160003/inputs"
case "$DURATION" in 400|500) ;; *) exit 2 ;; esac
test -s "$ROOT/validation/30s/SUCCESS"
test -s "$ROOT/validation/400s/SUCCESS"
test -s "$ROOT/validation/500s/SUCCESS"
ID=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
SEED=$((2026092500 + DURATION * 10 + SLURM_ARRAY_TASK_ID))
TRIAL="$ROOT/${DURATION}s/snn_trial_$ID"
DATA="$TRIAL/Data/Paper2_NetworkTuning/Fig1V4"
SAVE="$DATA/Paper2NWSimulationData"
mkdir -p "$DATA/Paper2PlotingData" "$DATA/GlobConv" "$SAVE" "$TRIAL/Paper2Figs"
ln -s "$BASE/snn_repo/Utils" "$TRIAL/Utils"
ln -s "$INPUT/AllMFPixPara_Paper2TuneFig1V4D2_torch_staged.mat" "$DATA/AllMFPixPara_Paper2TuneFig1V4D2.mat"
ln -s "$INPUT/LargeNWFixIni.mat" "$TRIAL/LargeNWFixIni.mat"
ln -s "$INPUT/Paper2DriveNW_Conn.mat" "$SAVE/Paper2DriveNW_Conn.mat"
cp "$ROOT/code/Paper2_Fig7Comp_NW_LDE_${DURATION}_bounded.m" "$TRIAL/Paper2_Fig7Comp_NW_LDE.m"
SHA=$(python3 -c 'import json,os; print(json.load(open(os.environ["HASH_FILE"]))[os.environ["DURATION"]])')
echo "$SHA  $TRIAL/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
cd "$TRIAL"
START=$(date +%s.%N)
matlab -batch "rng($SEED,'twister'); run('$TRIAL/Paper2_Fig7Comp_NW_LDE.m'); assert(T==$DURATION*1000 && dt==.1 && TPar==$DURATION && WinSize==9000); assert(~BlowUp && TSec==$DURATION && TimeN==10000);"
END=$(date +%s.%N)
ELAPSED=$(awk -v a="$START" -v b="$END" 'BEGIN{printf "%.9f",b-a}')
SIM_SECONDS=$(cat "$TRIAL/simulation_elapsed_seconds.txt")
for i in $(seq 1 "$DURATION"); do test -s "$SAVE/DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_${i}s_NewSmear.mat"; done
test -s "$SAVE/NWSimulationPix_0.0deg_NewSmear.mat"
echo "$SHA  $TRIAL/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
printf 'Repeat\tSeed\tSeconds\tSimulationSeconds\tSimulationMs\tdtMs\tUpdates\tCPUs\tMemoryGB\n%d\t%d\t%s\t%s\t%d000\t0.1\t%d0000\t16\t%d\n' "$SLURM_ARRAY_TASK_ID" "$SEED" "$ELAPSED" "$SIM_SECONDS" "$DURATION" "$DURATION" "$MEMORY_GB" > "$ROOT/${DURATION}s/results/trial_$ID.tsv"
