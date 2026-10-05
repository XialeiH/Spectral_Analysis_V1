#!/bin/zsh

set -u

cd /Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717

rsync -a --partial --progress \
  -e 'ssh -S ~/.ssh/torch-codex-cm' \
  torch:/scratch/xh2906/librarySCI_runs/l6ns_branch_atlas_20260803/branch_atlas_20260803/ \
  branch_atlas_20260803/ \
  > branch_atlas_20260803.transfer.log 2>&1

print -r -- "$?" > branch_atlas_20260803.transfer.exit
