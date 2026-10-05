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

EXPECTED_DRIVER_SHA=11981ea1ae5e8927e171930770ff7858997b4c78d01755fc0f294b7f21fccd2b
TRIAL=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
TRIAL_ROOT="$EXACT_ROOT/sanitized_trial_$TRIAL"
DATA_ROOT="$TRIAL_ROOT/Data/Paper2_NetworkTuning/Fig1V4"
SAVE_ROOT="$DATA_ROOT/Paper2NWSimulationData"
mkdir -p "$DATA_ROOT/Paper2PlotingData" "$DATA_ROOT/GlobConv" "$SAVE_ROOT" "$TRIAL_ROOT/Paper2Figs"
ln -s "$BASE_ROOT/snn_repo/Utils" "$TRIAL_ROOT/Utils"
ln -s "$EXACT_ROOT/inputs/AllMFPixPara_Paper2TuneFig1V4D2_torch_staged.mat" "$DATA_ROOT/AllMFPixPara_Paper2TuneFig1V4D2.mat"
ln -s "$EXACT_ROOT/inputs/LargeNWFixIni.mat" "$TRIAL_ROOT/LargeNWFixIni.mat"
ln -s "$EXACT_ROOT/inputs/Paper2DriveNW_Conn.mat" "$SAVE_ROOT/Paper2DriveNW_Conn.mat"
cp "$BASE_ROOT/snn_repo/Paper2_Fig7Comp_NW_LDE.m" "$TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m"
echo "$EXPECTED_DRIVER_SHA  $TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -

cd "$TRIAL_ROOT"
START=$(date +%s.%N)
matlab -batch "run('$TRIAL_ROOT/Paper2_Fig7Comp_NW_LDE.m')"
END=$(date +%s.%N)
ELAPSED=$(awk -v start="$START" -v end="$END" 'BEGIN {printf "%.6f",end-start}')

SEGMENT_COUNT=$(find "$SAVE_ROOT" -maxdepth 1 -type f -name 'DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_*s_NewSmear.mat' | wc -l | tr -d ' ')
if [[ "$SEGMENT_COUNT" -ne 10 ]]; then
    echo "Expected 10 DriveWkSp segment files, found $SEGMENT_COUNT" >&2
    exit 1
fi
test -f "$SAVE_ROOT/NWSimulationPix_0.0deg_NewSmear.mat"

OUT="$EXACT_ROOT/results/snn_paper2_exact_03HC_repeat_${TRIAL}.tsv"
printf 'FieldHC\tRepeat\tMethod\tSeconds\tSimulatedMs\tDtMs\tSteps\tCPUs\tPipeline\n' > "$OUT"
printf '3\t%d\tSNN\t%s\t10000\t0.1\t100000\t16\tPaper2OriginalExact\n' "$SLURM_ARRAY_TASK_ID" "$ELAPSED" >> "$OUT"
