from __future__ import annotations

import argparse
import math
import re
from collections import Counter
from pathlib import Path

import awkward as ak
import numpy as np


NUMERIC_WORDS = ("bool", "float", "double", "int", "short", "long")
INTEGER_WORDS = ("bool", "char", "int", "short", "long")
SAMPLE_STYLE_CYCLE = (
    ("C0", "/"),
    ("C1", "\\"),
    ("C2", "x"),
    ("C3", "."),
)
CHARGED_ABS_PDG = {11, 13, 211, 321, 2212}
MATCH_SOURCES = (
    ("nanoaod_raw_sdst", "nanoaod_raw_sdst.root"),
)
GEN_P4_BRANCHES = {
    "nanoaod_raw_sdst": (
        "GenPart_fourMomentum.fCoordinates.fX",
        "GenPart_fourMomentum.fCoordinates.fY",
        "GenPart_fourMomentum.fCoordinates.fZ",
    ),
}
MATCH_ALPHA_CUTS = np.linspace(0.0, 0.20, 201)
MATCH_PT_BINS = np.linspace(0.0, 60.0, 31)
MATCH_P_BINS = np.linspace(0.0, 60.0, 31)
MATCH_THETA_BINS = np.linspace(0.0, np.pi, 33)
MATCH_PHI_BINS = np.linspace(0.0, 2.0 * np.pi, 33)
MATCH_AXES = (
    ("pt", MATCH_PT_BINS, r"$p_T$ [GeV]"),
    ("p", MATCH_P_BINS, r"$p$ [GeV]"),
    ("theta", MATCH_THETA_BINS, r"$\theta$ [rad]"),
    ("phi", MATCH_PHI_BINS, r"$\phi$ [rad]"),
)
SOURCE_STYLES = {
    "nanoaod": ("C0", "."),
    "nanoaod_raw_sdst": ("C1", "/"),
    "nanoaod_raw_fadana": ("C2", "\\"),
}
GEN_COMPARE_PAIR_STYLES = {
    ("nanoaod", "nanoaod_raw_sdst"): SOURCE_STYLES["nanoaod"],
    ("nanoaod", "nanoaod_raw_fadana"): SOURCE_STYLES["nanoaod_raw_sdst"],
    ("nanoaod_raw_sdst", "nanoaod_raw_fadana"): SOURCE_STYLES["nanoaod_raw_fadana"],
}
FILE_GROUPS = (
    ("nanoaod", ("nanoaod.root",)),
    ("nanoaod_raw_sdst", ("nanoaod_raw_sdst.root",)),
    ("nanoaod_raw_fadana", ("nanoaod_raw_fadana.root",)),
    ("pythiastandalone", ("pythiastandalone.root", "nanoaod_pythiastandalone.root")),
)
TREE_NAMES_BY_SOURCE = {
    "nanoaod": ("Events", "t", "tgen", "tgenBefore"),
}
COMPARE_FILE_GROUPS = (
    ("nanoaod_raw_fadana", ("nanoaod_raw_fadana.root",)),
    ("nanoaod_raw_sdst", ("nanoaod_raw_sdst.root",)),
    ("nanoaod", ("nanoaod.root",)),
)


def add_samples_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--samples",
        nargs="+",
        help="sample directories to process (default: all directories under the input root)",
    )


def default_plot_root(project_root: Path, input_root: Path) -> Path:
    return project_root / "plots" / input_root.expanduser().resolve().name


def default_check_root(project_root: Path, input_root: Path) -> Path:
    return project_root / "data" / "check" / input_root.expanduser().resolve().name


def clean_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "branch"


def discover_samples(input_root: Path) -> tuple[str, ...]:
    return tuple(path.name for path in sorted(input_root.iterdir()) if path.is_dir())


def resolve_samples(input_root: Path, samples: tuple[str, ...] | list[str] | None = None) -> tuple[str, ...]:
    if samples:
        return tuple(samples)
    discovered = discover_samples(input_root)
    if not discovered:
        raise ValueError(f"No sample directories found under {input_root}")
    return discovered


def sample_styles(samples: tuple[str, ...] | list[str]) -> dict[str, tuple[str, str]]:
    return {
        sample: SAMPLE_STYLE_CYCLE[index % len(SAMPLE_STYLE_CYCLE)]
        for index, sample in enumerate(samples)
    }


def infer_expected_abs_pdg(status_events, pdg_events, parent_events) -> int | None:
    if hasattr(status_events, "layout"):
        status_events = ak.to_list(status_events)
    if hasattr(pdg_events, "layout"):
        pdg_events = ak.to_list(pdg_events)
    if hasattr(parent_events, "layout"):
        parent_events = ak.to_list(parent_events)

    counts = Counter()
    for statuses, pdgs, parents in zip(status_events, pdg_events, parent_events, strict=True):
        for idx, (status, pdg) in enumerate(zip(statuses, pdgs, strict=True)):
            abs_pdg = abs(int(pdg))
            if int(status) != 1 or abs_pdg not in CHARGED_ABS_PDG:
                continue
            if has_ancestor_pdg(pdgs, parents, idx, 22):
                continue
            counts[abs_pdg] += 1

    if not counts:
        return None
    return counts.most_common(1)[0][0]


def sample_grid(samples: tuple[str, ...] | list[str], subplot_width: float = 8.0, subplot_height: float = 7.0):
    import matplotlib.pyplot as plt

    n_samples = max(len(samples), 1)
    ncols = 2 if n_samples > 1 else 1
    nrows = math.ceil(n_samples / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(subplot_width * ncols, subplot_height * nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    flat_axes = list(axes.ravel())
    for ax in flat_axes[len(samples) :]:
        ax.set_visible(False)
    return fig, flat_axes[: len(samples)]


def build_branch_summary(sample_dir: Path, sample: str) -> dict:
    import uproot

    summary = {"sample": sample, "trees": {}}

    for source, file_names in FILE_GROUPS:
        try:
            path = resolve_file(sample_dir, file_names)
        except FileNotFoundError:
            continue

        with uproot.open(path) as root_file:
            for tree_name in TREE_NAMES_BY_SOURCE.get(source, ("Events",)):
                if tree_name not in root_file:
                    continue

                tree = root_file[tree_name]
                tree_summary = summary["trees"].setdefault(tree_name, {"branches": {}})
                for branch in sorted(tree.typenames()):
                    branch_summary = tree_summary["branches"].setdefault(branch, {"source": []})
                    branch_summary["source"].append(source)

    return summary


def resolve_file(sample_dir: Path, file_names: tuple[str, ...]) -> Path:
    for file_name in file_names:
        path = sample_dir / file_name
        if path.exists():
            return path
    raise FileNotFoundError(f"None of {file_names} exists under {sample_dir}")


def read_numeric_values(tree, branch: str, typename: str) -> np.ndarray:
    array = tree[branch].array(library="ak")
    if ak.fields(array):
        return np.array([], dtype=float)

    typename_lower = typename.lower()
    is_number = "char*" not in typename_lower and "string" not in typename_lower and any(
        word in typename_lower for word in NUMERIC_WORDS
    )
    if is_number:
        values = np.asarray(ak.to_numpy(ak.drop_none(ak.flatten(array, axis=None))), dtype=float).reshape(-1)
        return values[np.isfinite(values)]

    if typename.startswith("std::vector"):
        values = np.asarray(ak.to_numpy(ak.num(array, axis=1)), dtype=float)
        return values[np.isfinite(values)]

    return np.array([], dtype=float)


def is_integer_dtype(dtype: str) -> bool:
    dtype_lower = dtype.lower()
    if "char*" in dtype_lower or "string" in dtype_lower:
        return False
    if "float" in dtype_lower or "double" in dtype_lower:
        return False
    return any(word in dtype_lower for word in INTEGER_WORDS)


def is_probability_like(stats: list[dict]) -> bool:
    return all(
        not is_integer_dtype(item["dtype"])
        and float(item["min"]) >= 0.0
        and float(item["max"]) <= 1.0
        for item in stats
    )


def clamp(value: int, lo: int, hi: int) -> int:
    return max(lo, min(value, hi))


def charge_from_pdg(pdg_id: int) -> int:
    abs_pdg_id = abs(pdg_id)
    if abs_pdg_id in (11, 13):
        return -1 if pdg_id > 0 else 1
    if abs_pdg_id in (211, 321, 2212):
        return 1 if pdg_id > 0 else -1
    return 0


def valid_index(index: int, size: int) -> bool:
    return 0 <= index < size


def has_ancestor_pdg(pdgs, parents, start_idx: int, target_pdg: int, max_depth: int = 20) -> bool:
    seen = set()
    current = int(parents[start_idx])
    for _ in range(max_depth):
        if not valid_index(current, len(pdgs)) or current in seen:
            return False
        if int(pdgs[current]) == target_pdg:
            return True
        seen.add(current)
        current = int(parents[current])
    return False


def wrap_phi(phi: float) -> float:
    return (phi + np.pi) % (2.0 * np.pi) - np.pi


def phi_0_2pi(phi: float) -> float:
    return phi % (2.0 * np.pi)


def opening_angle(theta1: float, phi1: float, theta2: float, phi2: float) -> float:
    cos_angle = (
        np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
        + np.cos(theta1) * np.cos(theta2)
    )
    return float(np.arccos(np.clip(cos_angle, -1.0, 1.0)))


def is_track_fiducial(theta: float, pt: float, fiducial_cos_max: float, pt_min: float) -> bool:
    return np.isfinite(theta) and np.isfinite(pt) and pt > pt_min and abs(np.cos(theta)) < fiducial_cos_max


def divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan, dtype=float), where=denominator > 0)


def match_cut_dir(alpha_cut: float, fiducial_cos_max: float, pt_min: float) -> str:
    text = f"alpha_{alpha_cut:.3g}_cos_{fiducial_cos_max:.3g}_pt_{pt_min:.3g}"
    return text.replace(".", "p").replace("-", "m")


def add_delphi_label(ax) -> None:
    import mplhep as mh

    mh.label.exp_label(exp="DELPHI", llabel="Simulation", rlabel="LEP 1 (91.2 GeV)", loc=3, ax=ax)


def use_log_scale(counts_list: list[np.ndarray], threshold: float = 0.95) -> bool:
    for counts in counts_list:
        total = float(np.sum(counts))
        if total <= 0:
            continue
        biggest_bin = float(np.max(counts))
        if biggest_bin / total >= threshold and biggest_bin < total:
            return True
    return False


def set_hist_yaxis(ax, counts_list: list[np.ndarray], y_values: list[np.ndarray] | None = None) -> None:
    if not use_log_scale(counts_list):
        _, ymax = ax.get_ylim()
        ax.set_ylim(top=ymax * 1.3)
        return

    positive = []
    if y_values is None:
        positive = [float(value) for counts in counts_list for value in counts if value > 0]
    else:
        positive = [float(value) for values in y_values for value in values if value > 0]

    if not positive:
        _, ymax = ax.get_ylim()
        ax.set_ylim(top=ymax * 1.3)
        return

    ax.set_yscale("log")
    ax.set_ylim(bottom=max(min(positive) * 0.5, max(positive) * 1e-6), top=max(positive) * 1e2)


def histogram_bins(stats_by_label: dict) -> tuple[int, float, float]:
    stats = list(stats_by_label.values())

    if all(is_integer_dtype(item["dtype"]) for item in stats):
        lo = min(float(item["min"]) for item in stats)
        hi = max(float(item["max"]) for item in stats)
        lo = math.floor(lo)
        hi = math.ceil(hi)
        nbins = min(int(hi - lo + 1), 1000)
        return nbins, lo - 0.5, hi + 0.5

    if is_probability_like(stats):
        lo, hi = 0.0, 1.0
    else:
        lo = min(float(item["q005"]) for item in stats)
        hi = max(float(item["q995"]) for item in stats)

    if lo == hi:
        value = lo
        pad = max(abs(value) * 0.1, 0.5)
        return 100, value - pad, value + pad

    pad = max((hi - lo) * 1e-9, 1e-12)
    plot_lo = lo - pad
    plot_hi = hi + pad
    n_object = sum(int(item["n_object"]) for item in stats)
    iqr = sum((float(item["q75"]) - float(item["q25"])) * int(item["n_object"]) for item in stats) / n_object
    if iqr <= 0:
        return 100, plot_lo, plot_hi

    bin_width = 2.0 * iqr / (n_object ** (1.0 / 3.0))
    nbins = clamp(math.ceil((hi - lo) / bin_width), 20, 150)
    return nbins, plot_lo, plot_hi


def tree_branch_plot_path(output_root: Path, tree_name: str, branch: str) -> Path:
    path_parts = []
    for part in branch.split("."):
        path_parts.extend(part.split("_") if tree_name == "Events" else [part])
    path_parts = [clean_name(part) for part in path_parts if part]
    return output_root / tree_name / Path(*path_parts[:-1]) / f"{path_parts[-1]}.png"
