#!/bin/bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 RUN_ROOT SOURCE_TWO_PATHWAY_MAT" >&2
  exit 2
fi

RUN_ROOT="$1"
SOURCE_MAT="$2"
CODE_DIR="$RUN_ROOT/code"
OUT_ROOT="$RUN_ROOT/results"
LOG_DIR="$RUN_ROOT/logs"
mkdir -p "$OUT_ROOT" "$LOG_DIR"

if [[ ! -f "$SOURCE_MAT" ]]; then
  echo "Missing source MAT: $SOURCE_MAT" >&2
  exit 3
fi

COMMON=(--account=torch_pr_482_general --partition=cpu_short --qos=cpu48)
EXPORTS="ALL,IW_SOURCE_MAT=$SOURCE_MAT,IW_OUT_ROOT=$OUT_ROOT"

SMOKE_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=iw_smoke \
  --cpus-per-task=2 --mem=8G --time=00:30:00 \
  --output="$LOG_DIR/smoke_%j.out" --error="$LOG_DIR/smoke_%j.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_DIR' && matlab -nodisplay -nosplash -batch \"run_inhibition_weight_sweep_smoke\"")

EIG_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=iw_eigs \
  --dependency="afterok:$SMOKE_JOB" --array=1-15%12 \
  --cpus-per-task=4 --mem=24G --time=02:00:00 \
  --output="$LOG_DIR/eigs_%A_%a.out" --error="$LOG_DIR/eigs_%A_%a.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_DIR' && matlab -nodisplay -nosplash -batch \"run_inhibition_weight_eigenspectrum_task\"")

PLOT_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=iw_plot \
  --dependency="afterok:$EIG_JOB" --cpus-per-task=2 --mem=8G --time=00:30:00 \
  --output="$LOG_DIR/plot_%j.out" --error="$LOG_DIR/plot_%j.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_DIR' && matlab -nodisplay -nosplash -batch \"plot_inhibition_weight_eigenspectrum_sweep\"")

MANIFEST="$RUN_ROOT/submission_manifest.tsv"
printf 'stage\tjob_id\tdependency\toutput_root\n' > "$MANIFEST"
printf 'smoke\t%s\t\t%s\n' "$SMOKE_JOB" "$OUT_ROOT" >> "$MANIFEST"
printf 'full_eigenspectrum_array\t%s\tafterok:%s\t%s\n' "$EIG_JOB" "$SMOKE_JOB" "$OUT_ROOT" >> "$MANIFEST"
printf 'plots\t%s\tafterok:%s\t%s\n' "$PLOT_JOB" "$EIG_JOB" "$OUT_ROOT" >> "$MANIFEST"

echo "SMOKE_JOB=$SMOKE_JOB"
echo "EIG_JOB=$EIG_JOB"
echo "PLOT_JOB=$PLOT_JOB"
echo "MANIFEST=$MANIFEST"
