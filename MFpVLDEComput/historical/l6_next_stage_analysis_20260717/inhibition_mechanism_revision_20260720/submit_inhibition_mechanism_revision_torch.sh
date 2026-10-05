#!/bin/bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 RUN_ROOT SOURCE_ROOT" >&2
  exit 2
fi

RUN_ROOT="$1"
SOURCE_ROOT="$2"
CODE_ROOT="$RUN_ROOT/code"
OUTPUT_ROOT="$RUN_ROOT/results"
LOG_ROOT="$RUN_ROOT/logs"
mkdir -p "$OUTPUT_ROOT" "$LOG_ROOT"

if [[ ! -d "$SOURCE_ROOT/geometry_mat" ]]; then
  echo "Missing source geometry directory: $SOURCE_ROOT/geometry_mat" >&2
  exit 3
fi

COMMON=(--account=torch_pr_482_general --partition=cpu_short --qos=cpu48)
EXPORTS="ALL,INHR_SOURCE_ROOT=$SOURCE_ROOT,INHR_OUTPUT_ROOT=$OUTPUT_ROOT"

SWEEP_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=inhrev_sweep \
  --array=1-16%8 --cpus-per-task=4 --mem=8G --time=02:00:00 \
  --output="$LOG_ROOT/sweep_%A_%a.out" --error="$LOG_ROOT/sweep_%A_%a.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_ROOT' && matlab -nodisplay -nosplash -batch \"run_inhibition_revision_sweep_task\"")

MODE_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=inhrev_modes \
  --dependency="afterok:$SWEEP_JOB" --array=1-24%8 --cpus-per-task=4 \
  --mem=12G --time=02:00:00 \
  --output="$LOG_ROOT/modes_%A_%a.out" --error="$LOG_ROOT/modes_%A_%a.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_ROOT' && matlab -nodisplay -nosplash -batch \"run_inhibition_revision_mode_task\"")

PLOT_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=inhrev_plot \
  --dependency="afterok:$SWEEP_JOB:$MODE_JOB" --cpus-per-task=2 --mem=8G \
  --time=00:30:00 \
  --output="$LOG_ROOT/plot_%j.out" --error="$LOG_ROOT/plot_%j.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_ROOT' && matlab -nodisplay -nosplash -batch \"plot_inhibition_mechanism_revision\"")

MANIFEST="$RUN_ROOT/submission_manifest.tsv"
printf 'stage\tjob_id\tdependency\toutput_root\n' > "$MANIFEST"
printf 'expanded_weight_sweep\t%s\t\t%s\n' "$SWEEP_JOB" "$OUTPUT_ROOT" >> "$MANIFEST"
printf 'mode_atlas\t%s\tafterok:%s\t%s\n' "$MODE_JOB" "$SWEEP_JOB" "$OUTPUT_ROOT" >> "$MANIFEST"
printf 'summary_plots\t%s\tafterok:%s:%s\t%s\n' "$PLOT_JOB" "$SWEEP_JOB" "$MODE_JOB" "$OUTPUT_ROOT" >> "$MANIFEST"

echo "SWEEP_JOB=$SWEEP_JOB"
echo "MODE_JOB=$MODE_JOB"
echo "PLOT_JOB=$PLOT_JOB"
echo "MANIFEST=$MANIFEST"
