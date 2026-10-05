#!/usr/bin/env python3
"""Plot the two contrast-100 tuning-distance maps and selected conditions."""

from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path("/Users/xialeihuang/Desktop/Neuroscience_Project")
DATA_ROOT = PROJECT_ROOT / (
    "Spectral_Analysis/matlab-inserting_into_CG_model/Complete_Code_for_Paper3/"
    "NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/"
    "L6 and Inhibition/L6_and_Inhibition_offset/"
    "Figure6_0_FourContrast_FiveMetric"
)
DATASET = DATA_ROOT / "figure6_0_contrast100_dataset.mat"
REGRESSION = DATA_ROOT / "contrast_100_metric_regression.tsv"
OUTPUT = PROJECT_ROOT / (
    "Spectral_Analysis_Report/6_Pathway Compensation/"
    "6.16.0_Contrast100_Selected_Tuning_Heatmap_Conditions.pdf"
)

SELECTED_BETA6 = np.array(
    [0.056, 0.024, 0.240, 0.368, 0.320, 0.168, 0.168, 0.176]
)
SELECTED_BETAI = np.array(
    [0.016, 0.104, 0.096, 0.112, 0.352, 0.136, 0.112, 0.064]
)
SELECTIONS = (np.arange(0, 5), np.arange(5, 8))
METRIC_NAMES = ("pixel5x10CircularW1Deg", "pixel1x10CircularW1Deg")
PANEL_TITLES = (
    "E pixel (5,10) tuning-map distance",
    "E pixel (1,10) tuning-map distance",
)


def read_regression_slopes() -> dict[str, float]:
    table = np.genfromtxt(REGRESSION, delimiter="\t", names=True, dtype=None, encoding=None)
    return {
        str(row["metric"]): float(row["throughOriginSlopeBetaIOverBeta6"])
        for row in table
    }


def main() -> None:
    with h5py.File(DATASET, "r") as handle:
        beta_grid = np.asarray(handle["betaGrid"]).reshape(-1)
        # HDF5 reverses MATLAB dimensions: pixel, beta6, betaI.
        maps = np.asarray(handle["tuningWassersteinCircular"]).transpose(0, 2, 1)

    slopes = read_regression_slopes()
    cmap = mpl.colormaps["jet"].copy()
    cmap.set_bad((0.72, 0.72, 0.72, 1.0))
    step = beta_grid[1] - beta_grid[0]
    extent = (
        beta_grid[0] - step / 2,
        beta_grid[-1] + step / 2,
        beta_grid[0] - step / 2,
        beta_grid[-1] + step / 2,
    )

    mpl.rcParams.update(
        {
            "font.size": 13,
            "axes.titlesize": 15,
            "axes.titleweight": "bold",
            "axes.labelsize": 15,
            "axes.labelweight": "bold",
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "pdf.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(1, 2, figsize=(20, 8.5), constrained_layout=True)
    figure.suptitle(
        "Contrast 100: selected tuning-curve conditions",
        fontsize=20,
        fontweight="bold",
    )

    for panel, axis in enumerate(axes):
        values = np.ma.masked_invalid(maps[panel])
        values = np.ma.masked_less_equal(values, 0)
        finite_values = values.compressed()
        image = axis.imshow(
            values,
            origin="lower",
            extent=extent,
            interpolation="nearest",
            aspect="equal",
            cmap=cmap,
            norm=mpl.colors.LogNorm(
                vmin=float(finite_values.min()), vmax=float(finite_values.max())
            ),
            rasterized=True,
        )
        slope = slopes[METRIC_NAMES[panel]]
        x_line = beta_grid
        y_line = slope * x_line
        valid_line = (y_line >= beta_grid[0]) & (y_line <= beta_grid[-1])
        axis.plot(x_line[valid_line], y_line[valid_line], color="0.2", linewidth=1.6)

        for selection in SELECTIONS[panel]:
            axis.plot(
                SELECTED_BETA6[selection],
                SELECTED_BETAI[selection],
                marker="s",
                markersize=13,
                markerfacecolor="none",
                markeredgecolor="0.2",
                markeredgewidth=2.5,
            )
            axis.text(
                SELECTED_BETA6[selection],
                SELECTED_BETAI[selection],
                str(selection + 1),
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color="white",
                clip_on=True,
            )

        axis.set_title(f"{PANEL_TITLES[panel]}\nEmpirical slope = {slope:.4f}", pad=12)
        axis.set_xlabel(r"$\beta_6$ (L6 increase)")
        if panel == 0:
            axis.set_ylabel(r"$\beta_I$ (inhibition increase)")
        axis.set_xlim(0, 0.4)
        axis.set_ylim(0, 0.4)
        axis.set_xticks(np.arange(0, 0.41, 0.1))
        axis.set_yticks(np.arange(0, 0.41, 0.1))
        colorbar = figure.colorbar(image, ax=axis, fraction=0.048, pad=0.025)
        colorbar.set_label(r"Circular $W_1$", fontsize=13)
        colorbar.ax.set_title("deg", fontsize=11, pad=9)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", bbox_inches="tight")
    plt.close(figure)
    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    main()
