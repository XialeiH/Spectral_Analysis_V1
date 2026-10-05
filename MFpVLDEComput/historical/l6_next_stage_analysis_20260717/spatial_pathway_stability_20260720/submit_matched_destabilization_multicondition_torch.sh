#!/bin/bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 RUN_ROOT SOURCE_ROOT CODE_ROOT" >&2
  exit 2
fi

RUN_ROOT="$1"
SOURCE_ROOT="$2"
CODE_ROOT="$3"
OUTPUT_ROOT="$RUN_ROOT/results"
LOG_ROOT="$RUN_ROOT/logs"
mkdir -p "$OUTPUT_ROOT/figures" "$OUTPUT_ROOT/data" "$LOG_ROOT"

if [[ ! -d "$SOURCE_ROOT/geometry_mat" ]]; then
  echo "Missing source geometry directory: $SOURCE_ROOT/geometry_mat" >&2
  exit 3
fi

COMMON=(--account=torch_pr_482_general --partition=cpu_short --qos=cpu48)
EXPORTS="ALL,MATCHED_SOURCE_ROOT=$SOURCE_ROOT,MATCHED_OUTPUT_ROOT=$OUTPUT_ROOT"

ARRAY_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=matched_l6i_modes \
  --array=1-16%8 --cpus-per-task=4 --mem=24G --time=03:00:00 \
  --output="$LOG_ROOT/condition_%A_%a.out" --error="$LOG_ROOT/condition_%A_%a.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_ROOT' && matlab -nodisplay -nosplash -batch \"run_matched_destabilization_multicondition_task\"")

SUMMARY_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=matched_l6i_summary \
  --dependency="afterok:$ARRAY_JOB" --cpus-per-task=1 --mem=4G --time=00:20:00 \
  --output="$LOG_ROOT/summary_%j.out" --error="$LOG_ROOT/summary_%j.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_ROOT' && matlab -nodisplay -nosplash -batch \"aggregate_matched_destabilization_multicondition\"")

MANIFEST="$RUN_ROOT/submission_manifest.tsv"
printf 'stage\tjob_id\tarray\toutput_root\n' > "$MANIFEST"
printf 'matched_conditions\t%s\t1-16%%8\t%s\n' "$ARRAY_JOB" "$OUTPUT_ROOT" >> "$MANIFEST"
printf 'matched_summary\t%s\t-\t%s\n' "$SUMMARY_JOB" "$OUTPUT_ROOT" >> "$MANIFEST"

echo "ARRAY_JOB=$ARRAY_JOB"
echo "SUMMARY_JOB=$SUMMARY_JOB"
echo "MANIFEST=$MANIFEST"
