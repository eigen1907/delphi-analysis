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
    ("C0", "*"),
    ("C1", "."),
    ("C2", "/"),
    ("C3", "\\"),
    ("C4", "x"),
)
CHARGED_ABS_PDG = {11, 13, 211, 321, 2212}
MATCH_PT_BINS = np.linspace(0.0, 60.0, 31)
MATCH_P_BINS = np.linspace(0.0, 60.0, 31)
MATCH_THETA_BINS = np.linspace(0.0, np.pi, 33)
MATCH_PHI_BINS = np.linspace(0.0, 2.0 * np.pi, 33)
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


def add_color_scale_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--color-scale",
        choices=("shared", "per-sample"),
        default="per-sample",
        help="2D plot color scale mode (default: per-sample)",
    )


def default_plot_root(project_root: Path, input_root: Path, split: str = "train") -> Path:
    """Keep validation outputs under the campaign and numbered script stage."""
    path = input_root.expanduser().resolve()
    if path.suffix == ".root":
        path = path.parent
    if (path / f"{split}.root").is_file() or path.name in ("chunks", "merged"):
        path = path.parent
    return project_root / "output" / "plots" / path.name / "03_sample_validation" / split


def default_check_root(project_root: Path, input_root: Path) -> Path:
    name = input_root.expanduser().resolve().name
    if name in ("chunks", "merged"):
        name = input_root.expanduser().resolve().parent.name
    return project_root / "output" / "dataset" / "checks" / name


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


def normalize_counts(counts: np.ndarray) -> np.ndarray:
    total = np.sum(counts)
    return counts / total if total else counts.astype(float, copy=False)


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


def add_heatmap_label(ax, label: str, *, fontsize: float = 14, family: str | None = None) -> None:
    import matplotlib.patheffects as path_effects

    text_options = {
        "transform": ax.transAxes,
        "ha": "right",
        "va": "top",
        "multialignment": "left",
        "fontsize": fontsize,
        "color": "white",
        "path_effects": [path_effects.withStroke(linewidth=3.0, foreground="black")],
        "bbox": {
            "boxstyle": "round,pad=0.25",
            "facecolor": "black",
            "edgecolor": "none",
            "alpha": 0.45,
        },
    }
    if family is not None:
        text_options["family"] = family

    ax.text(0.94, 0.92, label, **text_options)


def build_branch_summary(sample_dir: Path, sample: str, split: str = "train") -> dict:
    import uproot
    from .dataset import SOURCE_TREES

    summary = {"sample": sample, "trees": {}}
    with uproot.open(sample_dir / f"{split}.root") as root_file:
        branches = summary["trees"].setdefault("Events", {"branches": {}})["branches"]
        for source, tree_path in SOURCE_TREES.items():
            tree = root_file[tree_path]
            for branch in sorted(tree.typenames()):
                branch_summary = branches.setdefault(branch, {"source": []})
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


def phi_0_2pi(phi: float) -> float:
    return phi % (2.0 * np.pi)


def divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan, dtype=float), where=denominator > 0)


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


# The three object studies share these small plotting functions.
FIGURE_SIZE = (11, 9)
SAMPLE_LABELS = dict(Zee=r'$\boldsymbol{Z\to e^+e^-}$',
                     Zmumu=r'$\boldsymbol{Z\to\mu^+\mu^-}$',
                     Ztautau=r'$\boldsymbol{Z\to\tau^+\tau^-}$',
                     ZKK=r'$\boldsymbol{Z\to K^+K^-}$',
                     Zpipi=r'$\boldsymbol{Z\to\pi^+\pi^-}$')


def set_style():
    import matplotlib.pyplot as plt
    import mplhep as mh

    mh.style.use(mh.styles.CMS)
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['patch.linewidth'] = 1.5
    plt.rcParams['legend.fontsize'] = 20
    plt.rcParams['legend.title_fontsize'] = 22


def finish(fig, ax, path, sample, xlabel, ylabel):
    import matplotlib.pyplot as plt
    import mplhep as mh

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ax.get_legend() is None and ax.get_legend_handles_labels()[0]:
        ax.legend(loc='upper right', frameon=False)
    label = next((label for name, label in SAMPLE_LABELS.items() if name in sample), sample)
    ax.text(0.02, 0.97, label, transform=ax.transAxes, va='top', fontweight='bold')
    ax.grid(alpha=0.2)
    mh.label.exp_label(exp='DELPHI', llabel='Simulation', rlabel='LEP 1 (91.2 GeV)', loc=0, ax=ax)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def clopper_pearson(numerator, denominator):
    """68.27% exact binomial intervals."""
    from scipy.stats import beta

    numerator, denominator = np.asarray(numerator), np.asarray(denominator)
    alpha = 1 - 0.6826894921370859
    lower = np.zeros_like(denominator, dtype=float)
    upper = np.ones_like(denominator, dtype=float)
    nonzero, not_all = numerator > 0, numerator < denominator
    lower[nonzero] = beta.ppf(alpha / 2, numerator[nonzero], denominator[nonzero] - numerator[nonzero] + 1)
    upper[not_all] = beta.ppf(1 - alpha / 2, numerator[not_all] + 1, denominator[not_all] - numerator[not_all])
    return lower, upper


def plot_efficiency(ax, x, matched, edges, label, color):
    denominator = np.histogram(x, bins=edges)[0]
    numerator = np.histogram(np.asarray(x)[matched], bins=edges)[0]
    shown = denominator > 0
    ratio = divide(numerator, denominator)
    lower, upper = clopper_pearson(numerator, denominator)
    centers = (edges[:-1] + edges[1:]) / 2
    ax.stairs(ratio, edges, baseline=None, color=color, label=label)
    ax.errorbar(centers[shown], ratio[shown],
                yerr=np.maximum([ratio[shown] - lower[shown], upper[shown] - ratio[shown]], 0),
                fmt='none', color=color, capsize=2)


def sample_files(input_root, filename):
    """Read merged samples or job files; keep the sample's dataset name in outputs."""
    found = False
    samples = [input_root] if (input_root / 'final_root').is_dir() else sorted(input_root.iterdir())
    for sample in samples:
        if not sample.is_dir():
            continue
        merged = sample / 'final_root' / filename
        if not merged.exists():
            merged = sample / filename
        files = [merged] if merged.exists() else sorted(sample.glob(f'final_root/job_*/{filename}'))
        if files:
            found = True
            dataset = sample.name if sample.name.startswith('20') else f'{input_root.name}_{sample.name}'
            yield dataset, files
    if not found:
        raise FileNotFoundError(f'No {filename} samples in {input_root}')


def plot_path(root, study, stage, dataset, source, collection, section=None, cut=None):
    """Stage first, then the dataset, source file, Reco definition and optional cut."""
    path = root / study / stage
    if section:
        path /= section
    path = path / dataset / source / collection
    return path / cut if cut else path
