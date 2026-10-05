"""Compare independent 500 s repeats and disjoint windows within each repeat."""
import argparse
import json
from pathlib import Path

import numpy as np


def hc_variability(states):
    # states: independent observations x S/C/I x pixel
    e = 0.6923 * states[:, 0] + 0.3077 * states[:, 1]
    i = states[:, 2]
    energy = np.mean(0.8 * e.mean(axis=0) ** 2 + 0.2 * i.mean(axis=0) ** 2)
    variance = np.mean(0.8 * e.var(axis=0, ddof=1) + 0.2 * i.var(axis=0, ddof=1))
    return {
        "HCnormSpS": float(np.sqrt(energy)),
        "HCstdSpS": float(np.sqrt(variance)),
        "relativeHCstd": float(np.sqrt(variance / energy)),
        "selectedPixelSD_SCI_SpS": states[:, :, 584].std(axis=0, ddof=1).tolist(),
    }


parser = argparse.ArgumentParser()
parser.add_argument("--root", type=Path, required=True)
parser.add_argument("--repeats", type=int, nargs="+", default=[1, 2, 3, 4, 5])
args = parser.parse_args()

counts = []
rates = []
for repeat in args.repeats:
    folder = args.root / f"500s_{repeat:03d}"
    check = json.loads((folder / "summary.json").read_text())
    assert check["duration"] == 500 and check["start"] == 0
    assert max(check["last9_max_abs"].values()) < 1e-10
    with np.load(folder / "counts.npz") as data:
        c = data["counts"]
        neurons = data["neurons"]
        assert c.shape == (3, 900, 500) and neurons.shape == (3, 900)
        assert int(data["duration"]) == 500 and int(data["start"]) == 0
        counts.append(c)
        rates.append(c / neurons[:, :, None])
for first in range(len(counts)):
    for second in range(first + 1, len(counts)):
        assert not np.array_equal(counts[first], counts[second]), "Identical spike counts"

rates = np.stack(rates)
report = {
    "durationSeconds": 500,
    "repeats": args.repeats,
    "targetRelativeHCstd": 0.005,
    "windows": {},
    "note": "Within-run disjoint-block SD is primary; across-trial SD is separate. Neither is the deterministic step residual. The original SNN output still averages only the last 9 s.",
}
for width in (100, 150, 200, 250, 300, 350, 400, 450):
    block_count = (500 - 50) // width
    entry = {
        "windowSeconds": width,
        "acrossTrialLastWindow": hc_variability(rates[:, :, :, -width:].mean(axis=-1)),
        "nonoverlappingBlocksPerTrial": block_count,
    }
    if block_count >= 2:
        block_means = rates[:, :, :, -block_count * width:].reshape(
            len(args.repeats), 3, 900, block_count, width
        ).mean(axis=-1)
        within = [
            hc_variability(block_means[trial].transpose(2, 0, 1))
            for trial in range(len(args.repeats))
        ]
        relative = np.array([item["relativeHCstd"] for item in within])
        entry["withinTrialDisjointBlocks"] = within
        entry["withinTrialMeanRelativeHCstd"] = float(relative.mean())
        entry["withinTrialRmsRelativeHCstd"] = float(np.sqrt(np.mean(relative ** 2)))
    report["windows"][str(width)] = entry

output = args.root / "500s_variability_report.json"
output.write_text(json.dumps(report, indent=2) + "\n")
print(f"Saved {output}")
for width, entry in report["windows"].items():
    within = entry.get("withinTrialMeanRelativeHCstd")
    print(
        f"{width} s: across-trial {100 * entry['acrossTrialLastWindow']['relativeHCstd']:.4f}%"
        + (f", mean within-trial {100 * within:.4f}%" if within is not None else "")
    )
