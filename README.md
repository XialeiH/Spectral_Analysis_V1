## Figure numbering

The main figure identifiers follow **main_01.pdf through main_05.pdf** and
**main_06.key**, not the older filenames used during analysis. Figure 6 was
checked against the Keynote A-F preview, as requested by the author.
For example, the historical `Figure_3D.m` singular-value panel is now
`Main_04_B.mlx`, and historical `Figure_4A.m` is now `Main_05_A.mlx`.
`panels.json` gives the complete source-to-panel mapping.

All plotting entry points are `.mlx` Live Scripts in the sectioned style of
`Paper3_Fig6.mlx` and `Paper3_Fig2_Fig7_8.mlx`: workspace setup, retained data,
figure construction, and export. They contain the plotting code, not just
links to another machine's scripts. Local helper functions remain in the
Live Script; shared functions are under `Utils/`.

## Run a panel

Open MATLAB R2025b with a JVM, change to this repository, and run, for example:

```matlab
run('Main_04_B.mlx')
```

Or open that Live Script in the MATLAB editor and choose Run. Outputs are
written beneath `Outputs/<panel identifier>/`. Inputs in `Data/` are not
overwritten. Like the Paper 3 plotting scripts, each entry starts with
`clearvars`; use a clean MATLAB workspace. Scripts that share a common original construction are grouped:
`Main_01_FG.mlx` produces the DNN/CG rows, and `Main_03_AB.mlx` produces
separate L6 and inhibition panels. Supplementary identifiers refer to the
current supplementary manuscript. `Additional_*` files are explicitly
historical alternatives, not additional main-figure labels.

To run all current panels in a dedicated MATLAB session:

```matlab
addpath('Utils');
reproduce_all;
verify_published_numbers;
```

Use `reproduce_all(true)` to also run the historical alternatives. Verify
the downloaded code/data bytes with `python3 Utils/verify_manifest.py`.

The exact typography uses Arial. On macOS the scripts use the system Arial
fonts. Font substitution on another platform can change text placement;
see `Utils/repro_fontdir.m`. The scripts require the MATLAB functions and
toolboxes used by their retained source; the execution report records the
tested environment and any unresolved panels.

## Repository layout

- Root `.mlx`: complete panel plotting entry points.
- `Utils/`: shared functions, retaining study-specific subdirectories to
  prevent different versions of the same function name from colliding.
- `Data/`: MAT, FIG, CSV and TSV inputs with checksums. Retained FIG files
  are graph-object/data inputs, not raster screenshots.
- `MFpVLDEComput/`: numerical drivers, Slurm scripts and provenance for
  upstream HPC computations. Historical launch scripts are not automatically
  submitted and require installation-specific paths/accounts.
- `Reference/`: current labeled figure assemblies used to verify numbering.
- `VALIDATION.md`: actual execution results, limitations and coverage.

The optional large upstream CG response table is split into byte chunks to
respect GitHub's single-file size limit. Before an upstream timing rerun,
restore it with `python3 Utils/restore_large_data.py`. Plotting does not need
this restoration. The script verifies the reconstructed SHA-256.

