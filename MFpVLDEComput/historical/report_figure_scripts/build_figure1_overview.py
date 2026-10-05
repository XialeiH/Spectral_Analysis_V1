#!/usr/bin/env python3
"""Build the panel-by-panel Figure 1 overview from existing report assets.

Panels sourced from measured/current report outputs retain their original data.
Panels that describe analyses not yet present in the report are explicitly marked
as protocol schematics so that the figure does not imply uncomputed results.
"""

from __future__ import annotations

import io
from pathlib import Path

import fitz
import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from PIL import Image


ROOT = Path("/Users/xialeihuang/Desktop/Neuroscience_Project")
REPORT = ROOT / "Spectral_Analysis_Report"
SCRIPT_ROOT = ROOT / "Spectral_Analysis/matlab-inserting_into_CG_model"
ARTIFACTS = SCRIPT_ROOT / "report_figure_artifacts/1_Neural_Network_Surrogate"
OUT = REPORT / "Figure_1_Final"
WORK = ROOT / "tmp/pdfs/figure1_final_sources"

NAVY = "#17365D"
BLUE = "#2878B5"
TEAL = "#2A9D8F"
ORANGE = "#E76F51"
GOLD = "#E9C46A"
PURPLE = "#7A5195"
GRAY = "#5B6573"
LIGHT = "#F4F7FA"
GRID = "#D9E1E8"


def render_pdf(pdf: Path, name: str, clip_fraction=None, dpi: int = 180) -> Path:
    """Render one PDF page, optionally clipping by page fractions."""
    WORK.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(pdf)
    page = doc[0]
    rect = page.rect
    if clip_fraction is not None:
        x0, y0, x1, y1 = clip_fraction
        rect = fitz.Rect(
            page.rect.x0 + x0 * page.rect.width,
            page.rect.y0 + y0 * page.rect.height,
            page.rect.x0 + x1 * page.rect.width,
            page.rect.y0 + y1 * page.rect.height,
        )
    pix = page.get_pixmap(matrix=fitz.Matrix(dpi / 72, dpi / 72), clip=rect, alpha=False)
    out = WORK / f"{name}.png"
    pix.save(out)
    doc.close()
    return out


def panel_frame(ax, letter: str, title: str, status: str | None = None):
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#C8D1DA")
        spine.set_linewidth(0.8)
    ax.text(0.018, 0.975, letter, transform=ax.transAxes, va="top", ha="left",
            fontsize=12, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.22", facecolor=NAVY, edgecolor=NAVY))
    ax.text(0.086, 0.972, title, transform=ax.transAxes, va="top", ha="left",
            fontsize=10.2, fontweight="bold", color=NAVY)
    if status:
        ax.text(0.982, 0.975, status, transform=ax.transAxes, va="top", ha="right",
                fontsize=7.2, color=GRAY, style="italic")


def add_box(ax, xy, wh, text, color, fontsize=8.2, subtitle=None):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                           facecolor=color, edgecolor="white", linewidth=1.2,
                           transform=ax.transAxes)
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h * 0.60, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize, fontweight="bold", color="white")
    if subtitle:
        ax.text(x + w / 2, y + h * 0.27, subtitle, transform=ax.transAxes,
                ha="center", va="center", fontsize=fontsize - 1.7, color="white")


def arrow(ax, start, end, label=None, color=GRAY, curve=0.0):
    arr = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12,
                          linewidth=1.4, color=color,
                          connectionstyle=f"arc3,rad={curve}", transform=ax.transAxes)
    ax.add_patch(arr)
    if label:
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2 + (0.05 if curve >= 0 else -0.05)
        ax.text(mx, my, label, transform=ax.transAxes, ha="center", va="center",
                fontsize=6.8, color=color,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1.0))


def draw_a(ax):
    panel_frame(ax, "A", "Model lineage", "schematic")
    xs = [0.04, 0.285, 0.53, 0.775]
    labels = [
        ("Conductance SNN", "spikes, synapses"),
        ("Local CG model", "population rates"),
        ("DNN response map", "differentiable Ψ"),
        ("Spatial surrogate", "fixed points + J"),
    ]
    colors = [PURPLE, ORANGE, TEAL, BLUE]
    for x, lab, col in zip(xs, labels, colors):
        add_box(ax, (x, 0.34), (0.18, 0.28), lab[0], col, subtitle=lab[1])
    for i in range(3):
        arrow(ax, (xs[i] + 0.18, 0.48), (xs[i + 1], 0.48), color=NAVY)
    ax.text(0.5, 0.16, "Biological spatial organization retained; local input–output evaluation is amortized",
            transform=ax.transAxes, ha="center", fontsize=7.5, color=GRAY)
    ax.text(0.5, 0.08, "The DNN is trained on coarse-grained local-response targets, not directly on experimental data.",
            transform=ax.transAxes, ha="center", fontsize=7.2, color=NAVY, fontweight="bold")


def draw_b(ax):
    panel_frame(ax, "B", "State and pathway diagram", "schematic")
    add_box(ax, (0.04, 0.58), (0.15, 0.22), "LGN", GOLD, subtitle="feedforward")
    add_box(ax, (0.32, 0.63), (0.14, 0.18), "S", BLUE, subtitle="simple E")
    add_box(ax, (0.32, 0.34), (0.14, 0.18), "C", TEAL, subtitle="complex E")
    add_box(ax, (0.66, 0.49), (0.14, 0.20), "I", ORANGE, subtitle="inhibitory")
    add_box(ax, (0.66, 0.16), (0.14, 0.18), "L6", PURPLE, subtitle="feedback")
    arrow(ax, (0.19, 0.69), (0.32, 0.72), "drive")
    arrow(ax, (0.19, 0.65), (0.32, 0.43))
    arrow(ax, (0.46, 0.71), (0.66, 0.61), "E kernel")
    arrow(ax, (0.46, 0.43), (0.66, 0.56))
    arrow(ax, (0.66, 0.52), (0.46, 0.67), "I kernel", curve=0.22)
    arrow(ax, (0.66, 0.49), (0.46, 0.40), curve=-0.18)
    arrow(ax, (0.46, 0.38), (0.66, 0.25), "L4E→L6")
    arrow(ax, (0.66, 0.28), (0.46, 0.62), "L6→L4", curve=-0.28)
    ax.text(0.87, 0.54, r"$E_u=w_S S+w_C C$", transform=ax.transAxes,
            ha="center", va="center", fontsize=10, color=NAVY,
            bbox=dict(boxstyle="round,pad=0.3", facecolor=LIGHT, edgecolor="#B7C4D0"))
    ax.text(0.87, 0.41, "mixed-E readout", transform=ax.transAxes,
            ha="center", fontsize=7, color=GRAY)


def load_paired():
    path = ARTIFACTS / "1.2_1.3_paper3_snn_h96_paired_data.mat"
    with h5py.File(path, "r") as f:
        g = f["paired"]
        return {k: np.asarray(g[k]).ravel() for k in g.keys()}


def draw_c(ax, paired):
    panel_frame(ax, "C", "SNN–NN-CG response parity", "measured: E/I")
    ax.set_axis_off()
    positions = [(0.07, 0.17, 0.25, 0.64), (0.38, 0.17, 0.25, 0.64), (0.70, 0.17, 0.25, 0.64)]
    for (x, y, w, h), pop, color in zip(positions[:2], ["E", "I"], [ORANGE, BLUE]):
        a = ax.inset_axes([x, y, w, h])
        xx = paired[f"{pop}_premodel"]
        yy = paired[f"{pop}_h96"]
        lo, hi = np.percentile(np.r_[xx, yy], [0.2, 99.8])
        a.scatter(xx, yy, s=1.3, color=color, alpha=0.58, edgecolors="none")
        a.plot([lo, hi], [lo, hi], color="#30343B", lw=1)
        a.set_xlim(lo, hi); a.set_ylim(lo, hi)
        a.set_title(f"{pop} parity", fontsize=7.5, pad=2)
        a.set_xlabel("SNN rate", fontsize=6.5); a.set_ylabel("NN-CG rate", fontsize=6.5)
        a.tick_params(labelsize=5.8, length=2)
        corr = np.corrcoef(xx, yy)[0, 1]
        rmse = np.sqrt(np.mean((yy - xx) ** 2))
        a.text(0.04, 0.95, f"r={corr:.3f}\nRMSE={rmse:.2f} Hz", transform=a.transAxes,
               va="top", fontsize=6.2, color=NAVY)
    a = ax.inset_axes(positions[2])
    for pop, color in [("E", ORANGE), ("I", BLUE)]:
        xx = paired[f"{pop}_premodel"]
        yy = paired[f"{pop}_h96"]
        rel = np.abs(yy - xx) / np.maximum(xx, 1e-6)
        bins = np.linspace(np.percentile(xx, 1), np.percentile(xx, 99), 8)
        ctr = 0.5 * (bins[:-1] + bins[1:])
        med = [np.median(rel[(xx >= bins[i]) & (xx < bins[i + 1])]) for i in range(len(ctr))]
        a.plot(ctr, med, "o-", ms=2.5, lw=1.1, color=color, label=pop)
    a.set_title("Rate-stratified error", fontsize=7.5, pad=2)
    a.set_xlabel("SNN rate (Hz)", fontsize=6.5); a.set_ylabel("median relative error", fontsize=6.5)
    a.tick_params(labelsize=5.8, length=2); a.legend(frameon=False, fontsize=6)
    ax.text(0.5, 0.06, "4 orientations × 2 stochastic replicates; 4,680 paired pixels per population",
            transform=ax.transAxes, ha="center", fontsize=6.8, color=GRAY)


def image_panel(ax, letter, title, image_path, status="existing report result", fit="contain"):
    panel_frame(ax, letter, title, status)
    ax.set_axis_off()
    img = Image.open(image_path)
    iax = ax.inset_axes([0.02, 0.04, 0.96, 0.87])
    iax.imshow(img)
    iax.axis("off")


def draw_e(ax, contrast_img, paired):
    panel_frame(ax, "E", "Contrast response and map-wide error", "existing + derived")
    ax.set_axis_off()
    left = ax.inset_axes([0.025, 0.08, 0.58, 0.81])
    left.imshow(Image.open(contrast_img)); left.axis("off")
    right = ax.inset_axes([0.66, 0.19, 0.31, 0.62])
    for pop, color in [("E", ORANGE), ("I", BLUE)]:
        x = paired[f"{pop}_premodel"]
        y = paired[f"{pop}_h96"]
        rel = 100 * (y - x) / np.maximum(x, 1e-6)
        right.hist(rel, bins=np.linspace(-35, 35, 40), histtype="step", density=True,
                   lw=1.3, color=color, label=pop)
    right.axvline(0, color="#333333", lw=0.8)
    right.set_xlabel("signed rate error (%)", fontsize=6.5)
    right.set_ylabel("density", fontsize=6.5)
    right.tick_params(labelsize=5.8, length=2)
    right.legend(frameon=False, fontsize=6)
    right.set_title("Full-map error", fontsize=7.5)


def draw_f(ax):
    panel_frame(ax, "F", "Derivative validation", "protocol schematic — run required")
    ax.set_axis_off()
    a1 = ax.inset_axes([0.08, 0.20, 0.38, 0.62])
    t = np.linspace(-1, 1, 80)
    a1.plot(t, t, color="#333333", lw=1)
    a1.scatter(t[::5], t[::5], s=12, facecolors="none", edgecolors=TEAL, lw=0.8)
    a1.set_xlabel("AD: $Jv$", fontsize=6.5); a1.set_ylabel("CG finite difference", fontsize=6.5)
    a1.tick_params(labelsize=5.5); a1.set_title("directional agreement", fontsize=7.2)
    a2 = ax.inset_axes([0.58, 0.20, 0.36, 0.62])
    h = np.logspace(-7, -1, 80)
    err = 3e-8 / h + 0.35 * h
    a2.loglog(h, err, color=PURPLE, lw=1.4)
    a2.axvspan(2e-4, 2e-3, color=GOLD, alpha=0.35)
    a2.set_xlabel("step size $h$", fontsize=6.5); a2.set_ylabel("relative error", fontsize=6.5)
    a2.tick_params(labelsize=5.5); a2.set_title("linear regime", fontsize=7.2)
    ax.text(0.5, 0.07, "Layout only: populate with parent-CG centered differences before publication.",
            transform=ax.transAxes, ha="center", fontsize=6.5, color=ORANGE, fontweight="bold")


def draw_g(ax):
    panel_frame(ax, "G", "Speed–accuracy frontier", "protocol schematic — benchmark required")
    ax.set_axis_off()
    a = ax.inset_axes([0.12, 0.18, 0.80, 0.65])
    pts = {"SNN": (3.0, 0.02, PURPLE), "parent CG": (2.0, 0.015, ORANGE), "DNN surrogate": (0.1, 0.022, TEAL)}
    for name, (x, y, col) in pts.items():
        a.scatter([x], [y], s=38, color=col)
        a.text(x, y * 1.12, name, fontsize=6.5, ha="center", color=col)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlim(0.04, 6); a.set_ylim(0.008, 0.05)
    a.set_xlabel("wall time (relative, log)", fontsize=6.5)
    a.set_ylabel("response error (log)", fontsize=6.5)
    a.tick_params(labelsize=5.5)
    ax.text(0.5, 0.07, "Positions are placeholders; benchmark all models on identical hardware and report CIs.",
            transform=ax.transAxes, ha="center", fontsize=6.5, color=ORANGE, fontweight="bold")


def draw_i(ax, img1, img2):
    panel_frame(ax, "I", "Recovery changes under enhanced L6", "existing trajectories; pair matching pending")
    ax.set_axis_off()
    for y, img, label in [(0.52, img1, "baseline L6"), (0.08, img2, "enhanced L6")]:
        sub = ax.inset_axes([0.035, y, 0.93, 0.37])
        sub.imshow(Image.open(img)); sub.axis("off")
        sub.text(0.01, 0.96, label, transform=sub.transAxes, ha="left", va="top",
                 fontsize=7.2, color=NAVY, fontweight="bold",
                 bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.0))
    ax.text(0.5, 0.025, "These existing trajectories are not yet the response-matched compensation pair in H.",
            transform=ax.transAxes, ha="center", fontsize=6.3, color=ORANGE, fontweight="bold")


def save_individual(draw_fn, letter, *args):
    fig, ax = plt.subplots(figsize=(7.2, 4.2), facecolor="white")
    draw_fn(ax, *args)
    fig.savefig(OUT / f"Figure_1{letter}_panel.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"Figure_1{letter}_panel.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    paired = load_paired()

    tuning_pdf = REPORT / "1_Neural Network Surrogate/1.1_H96_NN_contrast_and_orientation_tuning_curves.pdf"
    h_pdf = REPORT / "6_Pathway Compensation/6.0.1_Contrast100_Beta6_0p040_BetaI_0p024_Baseline_Comparison.pdf"
    i1_pdf = REPORT / "7_Perturbation/7.1_Perturbation_Transient_and_Return_Baseline_L6.pdf"
    i2_pdf = REPORT / "7_Perturbation/7.2_Perturbation_Transient_and_Return_Enhanced_L6_1p15x.pdf"
    j_pdf = REPORT / "6_Pathway Compensation/6.0_FourContrast_FiveMetric_L6_Inhibition_Compensation_Heatmaps.pdf"

    # D: retain the orientation-tuning section; E: retain the contrast-tuning section.
    d_img = render_pdf(tuning_pdf, "panel_D_orientation", (0.16, 0.18, 0.84, 0.97), dpi=170)
    e_img = render_pdf(tuning_pdf, "panel_E_contrast", (0.0, 0.0, 1.0, 0.22), dpi=190)
    h_img = render_pdf(h_pdf, "panel_H", dpi=150)
    i1_img = render_pdf(i1_pdf, "panel_I_baseline", (0.03, 0.025, 0.97, 0.22), dpi=190)
    i2_img = render_pdf(i2_pdf, "panel_I_enhanced", (0.03, 0.025, 0.97, 0.22), dpi=190)
    j_img = render_pdf(j_pdf, "panel_J", (0.08, 0.04, 0.92, 0.96), dpi=145)

    # Individual panels.
    save_individual(draw_a, "A")
    save_individual(draw_b, "B")
    save_individual(draw_c, "C", paired)
    save_individual(lambda ax: image_panel(ax, "D", "Orientation tuning across contrasts", d_img,
                                          "existing NN-CG result"), "D")
    save_individual(draw_e, "E", e_img, paired)
    save_individual(draw_f, "F")
    save_individual(draw_g, "G")
    save_individual(lambda ax: image_panel(ax, "H", "Flagship response-matched pair", h_img,
                                          "existing compensation result"), "H")
    save_individual(draw_i, "I", i1_img, i2_img)
    save_individual(lambda ax: image_panel(ax, "J", "Response-matched parameter landscape", j_img,
                                          "existing compensation grid"), "J")

    # Final five-row composite.
    fig = plt.figure(figsize=(20, 24), facecolor="white")
    gs = fig.add_gridspec(5, 6, height_ratios=[0.88, 1.18, 0.93, 1.18, 1.35],
                          hspace=0.16, wspace=0.10, left=0.025, right=0.985,
                          top=0.955, bottom=0.025)
    fig.suptitle("Figure 1 | A differentiable V1 surrogate reveals response-matched but dynamically distinct circuits",
                 fontsize=21, fontweight="bold", color=NAVY, y=0.985)
    fig.text(0.5, 0.963,
             "Measured report results are shown directly; missing derivative and hardware benchmarks are marked as protocols.",
             ha="center", fontsize=10, color=GRAY)

    draw_a(fig.add_subplot(gs[0, 0:3]))
    draw_b(fig.add_subplot(gs[0, 3:6]))
    draw_c(fig.add_subplot(gs[1, 0:3]), paired)
    image_panel(fig.add_subplot(gs[1, 3:6]), "D", "Orientation tuning across contrasts", d_img,
                "existing NN-CG result")
    draw_e(fig.add_subplot(gs[2, 0:2]), e_img, paired)
    draw_f(fig.add_subplot(gs[2, 2:4]))
    draw_g(fig.add_subplot(gs[2, 4:6]))
    image_panel(fig.add_subplot(gs[3, 0:3]), "H", "Flagship response-matched pair", h_img,
                "existing compensation result")
    draw_i(fig.add_subplot(gs[3, 3:6]), i1_img, i2_img)
    image_panel(fig.add_subplot(gs[4, 0:6]), "J", "Response-matched parameter landscape", j_img,
                "existing compensation grid")

    fig.savefig(OUT / "Figure_1_final_A_to_J.pdf")
    fig.savefig(OUT / "Figure_1_final_A_to_J.png", dpi=220)
    plt.close(fig)

    manifest = OUT / "Figure_1_panel_sources.txt"
    manifest.write_text(
        "A: generated model-lineage schematic\n"
        "B: generated state/pathway schematic\n"
        "C: generated from Paper3 SNN vs corrected-h96 paired data (E/I)\n"
        f"D: {tuning_pdf}\n"
        "E: contrast-tuning crop from D source plus paired-data error distribution\n"
        "F: protocol schematic; derivative-validation measurements not found in report\n"
        "G: protocol schematic; identical-hardware benchmark not found in report\n"
        f"H: {h_pdf}\n"
        f"I: {i1_pdf} and {i2_pdf}\n"
        f"J: {j_pdf}\n",
        encoding="utf-8",
    )
    print(OUT / "Figure_1_final_A_to_J.pdf")


if __name__ == "__main__":
    main()
