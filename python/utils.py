from __future__ import annotations

import math
import re
from pathlib import Path

import awkward as ak
import numpy as np


TREE_NAME = "Events"
META_TREE_NAME = "Meta"
NUMERIC_WORDS = ("bool", "float", "double", "int", "short", "long")


def project_root(start: Path | str | None = None) -> Path:
    path = Path.cwd() if start is None else Path(start)
    path = path.resolve()
    for candidate in (path, *path.parents):
        if (candidate / "python").is_dir() and (candidate / "scripts").is_dir() and (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"Could not find project root from {path}")


def sample_label(path: Path | str) -> str:
    name = Path(path).name
    return name.split("_")[-1]


def job_number(path: Path) -> int:
    return int(path.name.removeprefix("job_"))


def job_dirs(sample_dir: Path) -> list[Path]:
    final_root = sample_dir / "final_root"
    return sorted(final_root.glob("job_*"), key=job_number)


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")


def is_numeric(typename: str) -> bool:
    typename = typename.lower()
    return "char*" not in typename and "string" not in typename and any(word in typename for word in NUMERIC_WORDS)


def branch_values(tree, branch: str, typenames: dict[str, str]) -> tuple[np.ndarray, str]:
    array = tree[branch].array(library="ak")
    typename = typenames[branch]

    if ak.fields(array):
        return np.array([], dtype=float), "record"

    if is_numeric(typename):
        values = np.asarray(ak.to_numpy(ak.drop_none(ak.flatten(array, axis=None))), dtype=float).reshape(-1)
        return values[np.isfinite(values)], "values"

    if typename.startswith("std::vector"):
        values = np.asarray(ak.to_numpy(ak.num(array, axis=1)), dtype=float)
        return values[np.isfinite(values)], "multiplicity"

    return np.array([], dtype=float), "non_numeric"


def hist_bins(values: np.ndarray) -> tuple[int, float, float]:
    lo = float(np.min(values))
    hi = float(np.max(values))
    unique = np.unique(values)

    if np.allclose(values, np.round(values)) and unique.size <= 80:
        return max(1, unique.size), math.floor(lo) - 0.5, math.ceil(hi) + 0.5
    if lo == hi:
        pad = max(abs(lo) * 0.1, 0.5)
        return 40, lo - pad, hi + pad

    mean = float(np.mean(values))
    std = float(np.std(values))
    if std > 0:
        plot_lo = max(lo, mean - 3.0 * std)
        plot_hi = min(hi, mean + 3.0 * std)
        if plot_lo < plot_hi:
            return 80, plot_lo, plot_hi
    return 80, lo, hi
