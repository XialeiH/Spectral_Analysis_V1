# Upstream numerical computations

This folder preserves numerical MATLAB drivers and Slurm launch scripts;
the user-facing plotting interface is the root-level `.mlx` collection.
`source_manifest.json` identifies archived files by SHA-256. Python collection
and analysis sources have `python_source_manifest.json`. The
`unavailable_*_sources.json` lists record any unresolved original files;
an empty list means that archive pass recovered every listed source. Archived
scripts retain their original paths for provenance and are **not** a
portable, one-command Slurm installation. Configure paths, modules, account,
partition, resource requests and output locations before running them.

The final Main 1B rate-model drivers are under:

```
historical/report_figure_artifacts/Figure_1G/hcnorm_residual_allfields_20260925/code/
```

The archived 610-task mapping uses CG repeat counts
`[100,100,20,20,20,20,5,5]` and DNN counts
`[100,100,20,20,20,20,20,20]` at field sizes
`[3,4,6,8,10,20,30,40]`. The stopping test is
`HCnorm(next-current)/max(HCnorm(current),eps) < 0.005`
for three consecutive iterations, with `p=0.33` and at most 200 iterations.
Use the archived accelerated field-specific DNN drivers, not a generic
unaccelerated substitute. Reference comparisons are after the timer.

The 500-second SNN bounded-memory driver and separate variability-analysis
code are in `SNN500/`. The full-process timer ends before the later custom
variability calculation. Deleted raw SNN spike histories are not restored
or inferred by this archive. Repeating the simulation requires the original
network connectivity and initialization inputs in addition to this driver.

`Data/timing_runtime` and `Utils/timing_runtime` contain the retained small-field
runtime inputs/helpers that were available locally. Other field bundles can
be constructed by the archived bundle builder once its original network
parameter inputs are supplied. HPC runtime reproducibility is therefore
separate from reproducing the published plots from retained numerical data.
