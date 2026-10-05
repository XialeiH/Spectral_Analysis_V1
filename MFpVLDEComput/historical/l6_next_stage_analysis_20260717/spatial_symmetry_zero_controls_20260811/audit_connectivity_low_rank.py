#!/usr/bin/env python3
"""Audit symmetry, algebraic rank, and effective rank of spatial connectivity."""

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
SETUP_FILE = ROOT / "results_global_bifurcation_20260722/global_bifurcation_setup.mat"
SINGULAR_FILE = ROOT / (
    "results_full/section2_shared_modes/full_singular_spectrum_w_p0p000000.tsv"
)
MATRIX_FILE = ROOT / (
    "near_zero_clustering_20260722/results_canonical/canonical_matrix_summary.tsv"
)
OUTPUT_ROOT = Path(__file__).resolve().parent / "low_rank_audit"


def load_sparse(group, dimension):
    return csc_matrix(
        (
            np.asarray(group["data"]),
            np.asarray(group["ir"], dtype=np.int64),
            np.asarray(group["jc"], dtype=np.int64),
        ),
        shape=(dimension, dimension),
    )


def permutation(kind):
    index_map = np.arange(1600).reshape((40, 40), order="F")
    if kind == "translate_row1":
        transformed = np.roll(index_map, 1, axis=0)
    elif kind == "translate_column1":
        transformed = np.roll(index_map, 1, axis=1)
    elif kind == "rotate90":
        transformed = np.rot90(index_map, 1)
    elif kind == "rotate180":
        transformed = np.rot90(index_map, 2)
    elif kind == "reflect_rows":
        transformed = np.flipud(index_map)
    elif kind == "reflect_columns":
        transformed = np.fliplr(index_map)
    else:
        raise ValueError(kind)
    return transformed.reshape(-1, order="F")


def fft_singular_values(matrix):
    kernel = matrix[:, 0].toarray().ravel().reshape((40, 40), order="F")
    return np.sort(np.abs(np.fft.fft2(kernel).ravel()))[::-1]


def rank_metrics(values):
    squared = values**2
    energy = squared / squared.sum()
    tolerance = 1600 * np.finfo(float).eps * values[0]
    return {
        "algebraicRankAtMatlabTolerance": int(np.sum(values > tolerance)),
        "largestSingularValue": values[0],
        "smallestSingularValue": values[-1],
        "stableRank": squared.sum() / squared[0],
        "entropyEffectiveRank": np.exp(-np.sum(energy * np.log(np.maximum(energy, 1e-300)))),
        "countSingularLE0p01": int(np.sum(values <= 0.01)),
        "countSingularLE0p025": int(np.sum(values <= 0.025)),
        "countSingularLE0p05": int(np.sum(values <= 0.05)),
        "countSingularLE1PercentMax": int(np.sum(values <= 0.01 * values[0])),
    }


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    block_names = ["SS", "SC", "SI", "CS", "CC", "CI", "IS", "IC", "II"]
    symmetry_names = [
        "translate_row1",
        "translate_column1",
        "rotate90",
        "rotate180",
        "reflect_rows",
        "reflect_columns",
    ]
    rank_rows = []
    symmetry_rows = []
    frequency_rows = []
    matrices = {}
    frequency_axis = np.fft.fftfreq(40)
    frequency_radius = np.sqrt(
        frequency_axis[:, None] ** 2 + frequency_axis[None, :] ** 2
    ).ravel()
    with h5py.File(SETUP_FILE, "r") as handle:
        context = handle["setup/Context"]
        for name in block_names:
            matrix = load_sparse(context[f"C_{name}"], 1600)
            matrices[name] = matrix
            values = fft_singular_values(matrix)
            row = {"block": name, **rank_metrics(values)}
            row["rowSumCoefficientOfVariation"] = np.std(np.asarray(matrix.sum(1)).ravel()) / abs(
                np.mean(np.asarray(matrix.sum(1)).ravel())
            )
            rank_rows.append(row)
            kernel = matrix[:, 0].toarray().ravel().reshape((40, 40), order="F")
            fourier_gain = np.abs(np.fft.fft2(kernel)).ravel()
            relative_gain = fourier_gain / np.max(fourier_gain)
            frequency_rows.append(
                {
                    "block": name,
                    "spearmanRadiusVsLog10Gain": spearmanr(
                        frequency_radius, np.log10(np.maximum(relative_gain, 1e-300))
                    ).statistic,
                    "medianRelativeGainLowFrequency": np.median(
                        relative_gain[frequency_radius <= 0.15]
                    ),
                    "medianRelativeGainMidFrequency": np.median(
                        relative_gain[(frequency_radius > 0.15) & (frequency_radius <= 0.35)]
                    ),
                    "medianRelativeGainHighFrequency": np.median(
                        relative_gain[frequency_radius > 0.35]
                    ),
                    "highFrequencyFractionBelow1Percent": np.mean(
                        relative_gain[frequency_radius > 0.35] <= 0.01
                    ),
                }
            )
            denominator = np.linalg.norm(matrix.data)
            for symmetry in symmetry_names:
                order = permutation(symmetry)
                residual = matrix[order, :][:, order] - matrix
                symmetry_rows.append(
                    {
                        "block": name,
                        "symmetry": symmetry,
                        "relativeCommutatorError": np.linalg.norm(residual.data) / denominator,
                    }
                )

        kernel = np.asarray(context["L6Kernel"])
        padded = np.zeros((40, 40))
        padded[: kernel.shape[0], : kernel.shape[1]] = kernel
        rank_rows.append({"block": "L6Kernel", **rank_metrics(np.sort(np.abs(np.fft.fft2(padded).ravel()))[::-1])})

    proportionality_rows = []
    for left_name in block_names:
        left = matrices[left_name]
        denominator = left.multiply(left).sum()
        for right_name in block_names:
            right = matrices[right_name]
            scale = left.multiply(right).sum() / denominator
            residual = right - scale * left
            proportionality_rows.append(
                {
                    "referenceBlock": left_name,
                    "comparedBlock": right_name,
                    "bestScale": scale,
                    "relativeResidual": np.linalg.norm(residual.data) / np.linalg.norm(right.data),
                }
            )

    pd.DataFrame(rank_rows).to_csv(OUTPUT_ROOT / "connectivity_effective_rank.tsv", sep="\t", index=False)
    pd.DataFrame(symmetry_rows).to_csv(OUTPUT_ROOT / "connectivity_symmetry_audit.tsv", sep="\t", index=False)
    pd.DataFrame(frequency_rows).to_csv(
        OUTPUT_ROOT / "connectivity_fourier_attenuation.tsv", sep="\t", index=False
    )
    pd.DataFrame(proportionality_rows).to_csv(
        OUTPUT_ROOT / "connectivity_block_proportionality.tsv", sep="\t", index=False
    )

    singular = pd.read_csv(SINGULAR_FILE, sep="\t")["singularValue"].to_numpy()
    matrix_summary = pd.read_csv(MATRIX_FILE, sep="\t")
    baseline = matrix_summary.loc[matrix_summary["matrix"] == "J_baseline"].iloc[0]
    energy = singular**2 / np.sum(singular**2)
    comparison = pd.DataFrame(
        [
            {
                "matrix": "J_baseline",
                "algebraicRank": 4800,
                "stableRank": np.sum(singular**2) / singular[0] ** 2,
                "entropyEffectiveRank": np.exp(
                    -np.sum(energy * np.log(np.maximum(energy, 1e-300)))
                ),
                "countAbsEigLE0p01": int(baseline["countAbsEigLE1e2"]),
                "countSingularLE0p01": int(np.sum(singular <= 0.01)),
                "countAbsEigLE0p025": int(baseline["countAbsEigLE2p5e2"]),
                "countSingularLE0p025": int(np.sum(singular <= 0.025)),
                "countAbsEigLE0p05": int(baseline["countAbsEigLE5e2"]),
                "countSingularLE0p05": int(np.sum(singular <= 0.05)),
                "conditionNumber2": singular[0] / singular[-1],
            }
        ]
    )
    comparison.to_csv(OUTPUT_ROOT / "jacobian_singular_eigen_comparison.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
