#!/bin/bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 RUN_ROOT SOURCE_ROOT SWEEP_ROOT" >&2
  exit 2
fi

RUN_ROOT="$1"
SOURCE_ROOT="$2"
SWEEP_ROOT="$3"
CODE_ROOT="$RUN_ROOT/code"
OUTPUT_ROOT="$RUN_ROOT/results"
LOG_ROOT="$RUN_ROOT/logs"
mkdir -p "$OUTPUT_ROOT" "$LOG_ROOT"

for required in "$SOURCE_ROOT/geometry_mat" "$SOURCE_ROOT/l6_exact_boundary" "$SWEEP_ROOT/sweep_data"; do
  if [[ ! -e "$required" ]]; then
    echo "Missing required input: $required" >&2
    exit 3
  fi
done

COMMON=(--account=torch_pr_482_general --partition=cpu_short --qos=cpu48)
EXPORTS="ALL,L6ATLAS_SOURCE_ROOT=$SOURCE_ROOT,L6ATLAS_SWEEP_ROOT=$SWEEP_ROOT,L6ATLAS_OUTPUT_ROOT=$OUTPUT_ROOT"

ATLAS_JOB=$(sbatch --parsable "${COMMON[@]}" --job-name=l6_top4_atlas \
  --array=1-4%1 --cpus-per-task=4 --mem=16G --time=02:00:00 \
  --output="$LOG_ROOT/atlas_%A_%a.out" --error="$LOG_ROOT/atlas_%A_%a.err" \
  --export="$EXPORTS" \
  --wrap="module load matlab/2025b && cd '$CODE_ROOT' && matlab -nodisplay -nosplash -batch \"run_l6_multiangle_top4_atlas_task\"")

MANIFEST="$RUN_ROOT/submission_manifest.tsv"
printf 'stage\tjob_id\tarray\toutput_root\n' > "$MANIFEST"
printf 'multiangle_top4_atlas\t%s\t1-4%%1\t%s\n' "$ATLAS_JOB" "$OUTPUT_ROOT" >> "$MANIFEST"

echo "ATLAS_JOB=$ATLAS_JOB"
echo "MANIFEST=$MANIFEST"
