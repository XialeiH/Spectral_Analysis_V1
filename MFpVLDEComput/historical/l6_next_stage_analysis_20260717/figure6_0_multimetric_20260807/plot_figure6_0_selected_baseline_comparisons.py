#!/usr/bin/env python3
"""Plot baseline comparisons for three selected contrast-100 grid points."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATASET = Path(
    "/Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis/"
    "matlab-inserting_into_CG_model/Complete_Code_for_Paper3/"
    "NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/"
    "L6 and Inhibition/L6_and_Inhibition_offset/"
    "Figure6_0_FourContrast_FiveMetric/figure6_0_contrast100_dataset.mat"
)
SUMMARY = DATASET.with_name("contrast_100_point_summary.tsv")
PROJECT_OUTPUT = DATASET.parent
REPORT_OUTPUT = Path(
    "/Users/xialeihuang/Desktop/Neuroscience_Project/"
    "Spectral_Analysis_Report/6_Pathway Compensation"
)
SELECTED = [
    ("6.0.1", 0.040, 0.024),
    ("6.0.2", 0.120, 0.064),
    ("6.0.3", 0.256, 0.144),
]
MAP_SIDE = 40
POPULATION_SIZE = MAP_SIDE**2
C_WEIGHT = 0.3077
CORRECTED_TUNING = DATASET.with_name(
    "figure6_0_selected_corrected_tuning_curves.mat"
)


def e_map(vector: np.ndarray) -> np.ndarray:
    """Return the project E map using MATLAB column-major map ordering."""
    vector = np.asarray(vector).reshape(-1)
    excitatory = (1 - C_WEIGHT) * vector[:POPULATION_SIZE] + C_WEIGHT * vector[
        POPULATION_SIZE : 2 * POPULATION_SIZE
    ]
    return excitatory.reshape((MAP_SIDE, MAP_SIDE), order="F")


def shared_limits(first: np.ndarray, second: np.ndarray, signed: bool) -> tuple:
    values = np.concatenate((first.ravel(), second.ravel()))
    values = values[np.isfinite(values)]
    if signed:
        bound = np.max(np.abs(values))
        return (-bound, bound)
    return (np.min(values), np.max(values))


def add_heatmap(
    axis, values: np.ndarray, title: str, limits: tuple, signed: bool = False
):
    image = axis.imshow(
        values,
        origin="lower",
        interpolation="nearest",
        cmap="jet",
        vmin=limits[0],
        vmax=limits[1],
        aspect="equal",
    )
    axis.set_title(title, fontsize=13, fontweight="bold", pad=8)
    axis.set_xlabel("Map column", fontsize=11)
    axis.set_ylabel("Map row", fontsize=11)
    axis.set_xticks([0, 9, 19, 29, 39], [1, 10, 20, 30, 40])
    axis.set_yticks([0, 9, 19, 29, 39], [1, 10, 20, 30, 40])
    axis.tick_params(labelsize=9)
    return image


def load_data():
    with h5py.File(DATASET, "r") as handle:
        beta_grid = np.asarray(handle["betaGrid"]).reshape(-1)
        baseline_state = np.asarray(handle["setup/BaselineState"]).reshape(-1)
        baseline_cluster = np.asarray(
            handle["setup/BaselineClusterEnvelope"]
        ).reshape(-1)
        baseline_right = np.asarray(handle["setup/BaselineRightMode"]).reshape(-1)
        baseline_tuning = np.asarray(handle["setup/BaselineTuningCurves"])
        angles = np.asarray(handle["setup/FullAngles"]).reshape(-1)

        selected = []
        for label, requested_beta6, requested_beta_i in SELECTED:
            beta6_index = int(np.argmin(np.abs(beta_grid - requested_beta6)))
            beta_i_index = int(np.argmin(np.abs(beta_grid - requested_beta_i)))
            point_index = beta_i_index + beta6_index * beta_grid.size
            selected.append(
                {
                    "label": label,
                    "beta6": float(beta_grid[beta6_index]),
                    "beta_i": float(beta_grid[beta_i_index]),
                    "point_id": point_index + 1,
                    "state": np.asarray(handle["fixedPointStates"][point_index, :]),
                    "cluster": np.asarray(
                        handle["leadingEigenclusterEnvelopes"][point_index, :]
                    ),
                    "right": np.asarray(
                        handle["leadingRightSingularModes"][point_index, :]
                    ),
                    "tuning": np.asarray(handle["tuningCurves"][point_index, :, :]),
                }
            )

    result = {
        "angles": angles,
        "baseline_state": baseline_state,
        "baseline_cluster": baseline_cluster,
        "baseline_right": baseline_right,
        "baseline_tuning": baseline_tuning,
        "selected": selected,
    }
    if CORRECTED_TUNING.is_file():
        with h5py.File(CORRECTED_TUNING, "r") as corrected:
            result["baseline_tuning"] = np.asarray(
                corrected["baselineTuning"]
            )
            corrected_selected = np.asarray(corrected["selectedTuning"])
            for index, entry in enumerate(result["selected"]):
                entry["tuning"] = corrected_selected[index, :, :]
    return result


def plot_one(data: dict, selected: dict, metrics: pd.DataFrame):
    baseline_state_map = e_map(data["baseline_state"])
    selected_state_map = e_map(selected["state"])
    baseline_cluster_map = e_map(data["baseline_cluster"])
    selected_cluster_map = e_map(selected["cluster"])

    selected_right = selected["right"].copy()
    if np.dot(selected_right, data["baseline_right"]) < 0:
        selected_right *= -1
    baseline_right_map = np.abs(e_map(data["baseline_right"]))
    selected_right_map = np.abs(e_map(selected_right))

    limits = [
        shared_limits(baseline_state_map, selected_state_map, signed=False),
        shared_limits(baseline_cluster_map, selected_cluster_map, signed=False),
        (0, max(np.max(baseline_right_map), np.max(selected_right_map))),
    ]

    figure, axes = plt.subplots(
        2,
        4,
        figsize=(22, 10.5),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1, 1, 1, 1.45]},
    )
    images = []
    images.append(
        add_heatmap(
            axes[0, 0], baseline_state_map, "Baseline E fixed point", limits[0]
        )
    )
    add_heatmap(
        axes[1, 0], selected_state_map, "Selected E fixed point", limits[0]
    )
    images.append(
        add_heatmap(
            axes[0, 1],
            baseline_cluster_map,
            "Baseline E eigencluster envelope",
            limits[1],
        )
    )
    add_heatmap(
        axes[1, 1],
        selected_cluster_map,
        "Selected E eigencluster envelope",
        limits[1],
    )
    images.append(
        add_heatmap(
            axes[0, 2],
            baseline_right_map,
            "Baseline |E top right singular mode|",
            limits[2],
        )
    )
    add_heatmap(
        axes[1, 2],
        selected_right_map,
        "Selected |E top right singular mode|",
        limits[2],
    )

    colorbar_labels = [
        "E firing rate (Hz)",
        "Eigencluster envelope (a.u.)",
        "|Right singular mode| (a.u.)",
    ]
    for column, (image, label) in enumerate(zip(images, colorbar_labels)):
        colorbar = figure.colorbar(
            image, ax=[axes[0, column], axes[1, column]], shrink=0.88, pad=0.02
        )
        colorbar.set_label(label, fontsize=11)
        colorbar.ax.tick_params(labelsize=9)

    tuning_titles = ["E pixel (5,10)", "E pixel (1,10)"]
    selected_color = (0.85, 0.12, 0.10)
    for row in range(2):
        axis = axes[row, 3]
        axis.plot(
            data["angles"],
            data["baseline_tuning"][row, :],
            color="black",
            linestyle="--",
            linewidth=2.3,
            label=r"Baseline: $\beta_6=0$, $\beta_I=0$",
        )
        axis.plot(
            data["angles"],
            selected["tuning"][row, :],
            color=selected_color,
            linewidth=2.5,
            label=(
                rf"Selected: $\beta_6={selected['beta6']:.3f}$, "
                rf"$\beta_I={selected['beta_i']:.3f}$"
            ),
        )
        axis.set_title(
            f"{tuning_titles[row]} tuning comparison",
            fontsize=13,
            fontweight="bold",
            pad=8,
        )
        axis.set_xlim(0, 180)
        axis.set_xticks(np.arange(0, 181, 22.5))
        axis.set_xlabel("Orientation (deg)", fontsize=11)
        axis.set_ylabel("E firing rate (Hz)", fontsize=11)
        axis.grid(True, color="0.85", linewidth=0.8)
        axis.tick_params(labelsize=9)
        axis.legend(loc="best", fontsize=9, frameon=True)

    metric_row = metrics.loc[metrics["pointId"] == selected["point_id"]].iloc[0]
    figure.suptitle(
        (
            f"{selected['label']}  |  Contrast 100  |  "
            rf"$\beta_6={selected['beta6']:.3f}$, "
            rf"$\beta_I={selected['beta_i']:.3f}$  |  "
            f"fixed-point HC norm = {metric_row.fixedPointHC:.3f} Hz"
        ),
        fontsize=18,
        fontweight="bold",
    )

    beta6_text = f"{selected['beta6']:.3f}".replace(".", "p")
    beta_i_text = f"{selected['beta_i']:.3f}".replace(".", "p")
    base_name = (
        f"{selected['label']}_Contrast100_Beta6_{beta6_text}_"
        f"BetaI_{beta_i_text}_Baseline_Comparison"
    )
    for output_root in (PROJECT_OUTPUT, REPORT_OUTPUT):
        output_root.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_root / f"{base_name}.pdf", bbox_inches="tight")
        figure.savefig(
            output_root / f"{base_name}.png", dpi=300, bbox_inches="tight"
        )
    plt.close(figure)
    return base_name, metric_row


def main():
    data = load_data()
    metrics = pd.read_csv(SUMMARY, sep="\t")
    rows = []
    for selected in data["selected"]:
        base_name, metric_row = plot_one(data, selected, metrics)
        rows.append(
            {
                "figure": selected["label"],
                "fileBase": base_name,
                "pointId": selected["point_id"],
                "beta6": selected["beta6"],
                "betaI": selected["beta_i"],
                "fixedPointHC": metric_row.fixedPointHC,
                "singularModeHC": metric_row.singularModeHC,
                "eigenclusterHC": metric_row.eigenclusterHC,
                "pixel5x10CircularW1Deg": metric_row.pixel5x10CircularW1Deg,
                "pixel1x10CircularW1Deg": metric_row.pixel1x10CircularW1Deg,
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(
        PROJECT_OUTPUT / "6.0.1-6.0.3_Selected_Condition_Manifest.tsv",
        sep="\t",
        index=False,
    )
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
