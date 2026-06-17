from __future__ import annotations

from pathlib import Path

import hist
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import pandas as pd
import uproot

from utils import TREE_NAME, branch_values, hist_bins, safe_name


PARTS = (
    ("sdst", "SDST_"),
    ("raw_sdst", "RAWSDST_"),
    ("raw_fadana", "RAWFADANA_"),
)
GROUP_ALIASES = {"TracRaw": "TrackRaw"}

SampleFile = tuple[str, Path]
OpenSample = tuple[str, object, set[str], dict[str, str]]


def select_branches(samples: list[SampleFile], prefix: str = "") -> list[str]:
    branches = set()
    for _, path in samples:
        with uproot.open(path) as root_file:
            branches.update(
                branch for branch in root_file[TREE_NAME].keys()
                if branch.startswith(prefix)
            )
    return sorted(branches)


def open_samples(samples: list[SampleFile]) -> list[OpenSample]:
    opened = []
    for label, path in samples:
        tree = uproot.open(path)[TREE_NAME]
        opened.append((label, tree, set(tree.keys()), tree.typenames()))
    return opened


def variable_name(branch: str, prefix: str) -> str:
    return branch[len(prefix):] if branch.startswith(prefix) else branch


def group_name(branch: str, prefix: str) -> str:
    group = variable_name(branch, prefix).split("_", 1)[0]
    if group.startswith("n") and len(group) > 1 and group[1].isupper():
        group = group[1:]
    return GROUP_ALIASES.get(group, group)


def branch_summary(samples: list[OpenSample], branches: list[str]) -> pd.DataFrame:
    rows = []
    for branch in branches:
        missing = []
        empty = []
        typenames_seen = set()
        chunks = []

        for label, tree, keys, typenames in samples:
            if branch not in keys:
                missing.append(label)
                continue

            values, _ = branch_values(tree, branch, typenames)
            typenames_seen.add(typenames[branch])
            if values.size:
                chunks.append(values)
            else:
                empty.append(label)

        values = np.concatenate(chunks) if chunks else np.array([])
        rows.append(
            {
                "branch": branch,
                "typename": " | ".join(sorted(typenames_seen)),
                "missing_files": ", ".join(missing),
                "empty_files": ", ".join(empty),
                "entries": int(values.size),
                "min": float(np.min(values)) if values.size else np.nan,
                "max": float(np.max(values)) if values.size else np.nan,
                "mean": float(np.mean(values)) if values.size else np.nan,
            }
        )
    return pd.DataFrame(rows)


def save_branch_plot(branch: str, samples: list[OpenSample], output_dir: Path, prefix: str) -> dict[str, object]:
    variable = variable_name(branch, prefix)
    values_by_sample = {}

    for label, tree, keys, typenames in samples:
        if branch in keys:
            values, _ = branch_values(tree, branch, typenames)
            values_by_sample[label] = values

    non_empty_values = [values for values in values_by_sample.values() if values.size]
    values = np.concatenate(non_empty_values) if non_empty_values else np.array([])

    fig, ax = plt.subplots(figsize=(12, 9))
    if values.size:
        nbins, lo, hi = hist_bins(values)
        for label, sample_values in values_by_sample.items():
            if sample_values.size:
                h = hist.Hist.new.Reg(nbins, lo, hi, name="value", label=variable).Double()
                h.fill(value=sample_values)
                mh.histplot(h, ax=ax, histtype="step", label=label, yerr=True)
        ax.legend()

    ax.set_xlabel(variable)
    ax.set_ylabel("Entries")
    mh.label.exp_label(
        exp="DELPHI",
        text="Private Work",
        rlabel=r"LEP1, $\sqrt{s}\sim 91.25$ GeV",
        loc=0,
        ax=ax,
    )
    fig.tight_layout()

    path = output_dir / group_name(branch, prefix) / f"{safe_name(branch)}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)

    return {
        "branch": branch,
        "variable": variable,
        "entries": int(values.size),
        "empty": not bool(values.size),
        "path": str(path),
    }


def plot_part(part: str, prefix: str, input_samples: list[SampleFile], output_root: Path) -> dict[str, object]:
    output_dir = output_root / part
    output_dir.mkdir(parents=True, exist_ok=True)

    branches = select_branches(input_samples, prefix)
    samples = open_samples(input_samples)
    summary = branch_summary(samples, branches)
    manifest = pd.DataFrame(save_branch_plot(branch, samples, output_dir, prefix) for branch in branches)

    part_summary = pd.DataFrame([
        {
            "part": part,
            "branches": len(branches),
            "plotted": len(manifest),
            "empty": int(manifest["empty"].sum()),
            "output_dir": str(output_dir),
        }
    ])

    summary.to_csv(output_dir / "branch_summary.csv", index=False)
    manifest.to_csv(output_dir / "plot_manifest.csv", index=False)
    part_summary.to_csv(output_dir / "plot_variables_summary.csv", index=False)

    row = part_summary.iloc[0]
    print(f"[{part}] plotted {row.plotted} / {row.branches} branches; {row.empty} empty plots")
    return {"branches": branches, "branch_table": summary, "manifest": manifest, "summary": part_summary}


def plot_variables(input_samples: list[SampleFile], output_root: Path) -> dict[str, object]:
    mh.style.use(mh.styles.CMS)
    output_root.mkdir(parents=True, exist_ok=True)

    all_branches = select_branches(input_samples)
    all_branches_table = pd.DataFrame({"branch": all_branches})
    all_branches_table.to_csv(output_root / "all_branches.csv", index=False)

    print(f"Found {len(all_branches)} branches across {len(input_samples)} ROOT files")
    print(f"Writing plots and tables under {output_root}")

    results = {
        part: plot_part(part, prefix, input_samples, output_root)
        for part, prefix in PARTS
    }
    plot_summaries = pd.concat([result["summary"] for result in results.values()], ignore_index=True)
    plot_summaries.to_csv(output_root / "plot_variables_summary.csv", index=False)

    return {
        "all_branches": all_branches,
        "all_branches_table": all_branches_table,
        "parts": results,
        "summary": plot_summaries,
    }
