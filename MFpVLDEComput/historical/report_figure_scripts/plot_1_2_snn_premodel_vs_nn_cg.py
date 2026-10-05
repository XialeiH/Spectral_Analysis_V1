#!/usr/bin/env python3
"""Replot the published Paper-2 SNN/CG pixelwise firing-rate comparison.

The public Paper-2 repository contains the plotting and simulation code, but
not the four NWSimulationPix_* files needed to reevaluate the current neural
surrogate on the original pixelwise inputs.  This script therefore digitizes
the published paired CG-target points from the supplied reference image.  The
vertical-axis label states that these are the targets represented by the NN
surrogate; it must not be interpreted as a fresh NN evaluation.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REFERENCE_IMAGE = Path(
    "/var/folders/gb/k9rrr4cs2f3cnn_z4rf_v7080000gn/T/"
    "codex-clipboard-3db30f39-90ee-4591-b6e8-1a3b65af95ba.png"
)
OUTPUT_DIR = Path(
    "/Users/xialeihuang/Desktop/Neuroscience_Project/"
    "Spectral_Analysis_Report/1_Neural Network Surrogate"
)
OUTPUT_PDF = OUTPUT_DIR / "1.2_SNN_premodel_vs_NN_CG_pixelwise_firing_rates.pdf"
OUTPUT_PNG = OUTPUT_DIR / "1.2_SNN_premodel_vs_NN_CG_pixelwise_firing_rates.png"


def digitized_pixels(image, crop, population):
    x0, y0, x1, y1 = crop
    rgb = image[y0:y1, x0:x1].astype(float)
    red, green, blue = np.moveaxis(rgb, -1, 0)
    if population == "E":
        mask = (
            (red > 180)
            & (red > 1.4 * green)
            & (red > 1.4 * blue)
            & (green < 170)
        )
    else:
        mask = (
            (blue > 150)
            & (blue > 1.3 * red)
            & (blue > 1.15 * green)
            & (red < 180)
        )
    rows, columns = np.nonzero(mask)
    return columns + x0, rows + y0


def pixels_to_data(columns, rows, crop, limits):
    x0, y0, x1, y1 = crop
    low, high = limits
    x = low + (columns - x0) * (high - low) / (x1 - x0)
    y = high - (rows - y0) * (high - low) / (y1 - y0)
    return x, y


def draw_panel(axis, x, y, limits, population):
    color = "#ff2020" if population == "E" else "#355cff"
    axis.scatter(
        x, y, s=0.55, marker=".", color=color, alpha=0.48,
        linewidths=0, rasterized=True
    )

    low, high = limits
    line_x = np.linspace(low, high, 500)
    axis.plot(line_x, line_x, color="black", linewidth=1.4, linestyle=(0, (10, 8)))
    axis.plot(line_x, (4 / 3) * line_x, color="black", linewidth=2.0, linestyle="--")
    axis.plot(line_x, (3 / 4) * line_x, color="black", linewidth=2.0, linestyle="--")

    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_title(f"Pixelwise {population} firing rates", fontsize=20, fontweight="bold", pad=12)
    axis.set_xlabel("Premodel Fr (Hz)", fontsize=18, labelpad=8)
    if population == "E":
        axis.set_ylabel("NN-CG Model Fr (Hz)", fontsize=18, labelpad=8)
    axis.tick_params(axis="both", labelsize=15, width=1.2, length=5, direction="in")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.3)
    axis.spines["bottom"].set_linewidth(1.3)


def main():
    image = np.asarray(Image.open(REFERENCE_IMAGE).convert("RGB"))
    e_crop = (169, 187, 891, 910)
    i_crop = (1110, 187, 1832, 910)

    e_columns, e_rows = digitized_pixels(image, e_crop, "E")
    i_columns, i_rows = digitized_pixels(image, i_crop, "I")
    e_x, e_y = pixels_to_data(e_columns, e_rows, e_crop, (5, 40))
    i_x, i_y = pixels_to_data(i_columns, i_rows, i_crop, (25, 110))

    figure, axes = plt.subplots(1, 2, figsize=(13.2, 6.1), constrained_layout=False)
    draw_panel(axes[0], e_x, e_y, (5, 40), "E")
    draw_panel(axes[1], i_x, i_y, (25, 110), "I")
    axes[0].set_xticks([10, 20, 30, 40])
    axes[0].set_yticks([5, 10, 15, 20, 25, 30, 35, 40])
    axes[1].set_xticks([40, 60, 80, 100])
    axes[1].set_yticks([30, 40, 50, 60, 70, 80, 90, 100, 110])
    figure.subplots_adjust(left=0.075, right=0.985, bottom=0.14, top=0.88, wspace=0.28)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PDF, bbox_inches="tight")
    figure.savefig(OUTPUT_PNG, dpi=400, bbox_inches="tight")
    plt.close(figure)
    print(OUTPUT_PDF)
    print(OUTPUT_PNG)


if __name__ == "__main__":
    main()
