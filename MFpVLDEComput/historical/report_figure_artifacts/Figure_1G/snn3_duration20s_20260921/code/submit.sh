#!/bin/bash
set -euo pipefail
export OUTPUT_ROOT=/scratch/xh2906/librarySCI_runs/figure1b_snn3_duration20s_20260921
export BASE_ROOT=/scratch/xh2906/librarySCI_runs/figure1g2_dense_cold_equilibrium_20260829
export EXACT_ROOT="$BASE_ROOT/snn_paper2_original_exact_20260901_160003"
export DRIVER_SHA=b3958be55c41dba44f674de73e3dfb8bf37e2b552c7f4b145b824872197f4c6b
test ! -e "$OUTPUT_ROOT/submission.txt"
echo "11981ea1ae5e8927e171930770ff7858997b4c78d01755fc0f294b7f21fccd2b  $BASE_ROOT/snn_repo/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
echo "$DRIVER_SHA  $OUTPUT_ROOT/code/Paper2_Fig7Comp_NW_LDE.m" | sha256sum -c -
python3 -c 'import pathlib,os; b=pathlib.Path(os.environ["BASE_ROOT"])/"snn_repo/Paper2_Fig7Comp_NW_LDE.m"; p=pathlib.Path(os.environ["OUTPUT_ROOT"])/"code/Paper2_Fig7Comp_NW_LDE.m"; s=b.read_bytes(); old=b"T = 10000; dt = 0.1; TPar = 10;"; assert s.count(old)==1; assert p.read_bytes().rstrip(b"\n")==s.replace(old,b"T = 20000; dt = 0.1; TPar = 20;").rstrip(b"\n"); print("Verified: only duration/partition-count line and final newline differ from archived original.")'
bash -n "$OUTPUT_ROOT/code/snn_20s.sbatch"
sbatch --parsable --export=ALL --output="$OUTPUT_ROOT/logs/snn_%A_%a.out" --error="$OUTPUT_ROOT/logs/snn_%A_%a.err" "$OUTPUT_ROOT/code/snn_20s.sbatch" | tee "$OUTPUT_ROOT/submission.txt"
