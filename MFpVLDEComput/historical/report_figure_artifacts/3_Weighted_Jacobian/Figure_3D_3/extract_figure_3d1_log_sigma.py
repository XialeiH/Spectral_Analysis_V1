#!/usr/bin/env python3
"""Recover the plotted log10(sigma_j) values from the vector Figure 3D.1 PDF."""

from pathlib import Path

import numpy as np
import pdfplumber


SOURCE = Path(
    "/Users/xialeihuang/Desktop/Neuroscience_Project/"
    "Spectral_Analysis_Report/Figures/Figure 3D.1.pdf"
)
OUTPUT = Path(__file__).with_name("figure_3d3_log_sigma.tsv")

BLUE = np.array((0.262745098, 0.576470588, 0.764705882))
RED = np.array((0.698039215, 0.094117647, 0.16862745))


def is_marker(curve):
    color = curve.get("non_stroking_color")
    if color is None:
        return False
    color = np.asarray(color, dtype=float)
    width = float(curve.get("width", np.nan))
    return (
        abs(width - 5.47722525) < 1e-3
        and np.linalg.norm(color - BLUE) < 1e-4
    ) or (
        abs(width - 7.348469) < 1e-3
        and np.linalg.norm(color - RED) < 1e-4
    )


with pdfplumber.open(SOURCE) as pdf:
    centers = np.array(
        [
            (curve["x0"] + curve["x1"]) / 2
            for curve in pdf.pages[0].curves
            if is_marker(curve)
        ]
    )

if centers.size != 9600:
    raise RuntimeError(f"Expected 9,600 markers, found {centers.size}.")

left_pdf_ticks = np.array([175.634, 273.831, 372.027, 470.223, 568.420])
right_pdf_ticks = np.array([777.720, 876.007, 974.294, 1072.581, 1170.868])
data_ticks = np.array([-3.0, -2.0, -1.0, 0.0, 1.0])

rows = []
for panel, mask, pdf_ticks in (
    ("Baseline Gaussian Connectivity", centers < 650, left_pdf_ticks),
    ("Row-Sum-Matched Flattened Connectivity", centers > 650, right_pdf_ticks),
):
    slope, intercept = np.polyfit(pdf_ticks, data_ticks, 1)
    values = slope * centers[mask] + intercept
    if values.size != 4800:
        raise RuntimeError(f"Expected 4,800 values for {panel}, found {values.size}.")
    rows.extend((panel, value) for value in values)

with OUTPUT.open("w", encoding="ascii") as stream:
    stream.write("panel\tlog10_sigma\n")
    for panel, value in rows:
        stream.write(f"{panel}\t{value:.12g}\n")

print(OUTPUT)
