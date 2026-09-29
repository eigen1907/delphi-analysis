from __future__ import annotations

import csv
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from delphi_analysis.plot_utils import (
    add_delphi_label,
    add_heatmap_label,
    resolve_samples,
    sample_grid,
    sample_styles,
    set_hist_yaxis,
)


SOURCE_FILE = "nanoaod_raw_sdst.root"
REFERENCE_FILE = "nanoaod.root"
TREE_NAME = "Events"
RICH_BRANCHES = (
    "HaidRaw_richQuality",
    "HaidRaw_thetaGas",
    "HaidRaw_sigmaGas",
    "HaidRaw_nphGas",
    "HaidRaw_nepGas",
    "HaidRaw_flagGas",
    "HaidRaw_thetaLiq",
    "HaidRaw_sigmaLiq",
    "HaidRaw_nphLiq",
    "HaidRaw_nepLiq",
    "HaidRaw_flagLiq",
)
FLOAT_BRANCHES = {
    "HaidRaw_thetaGas",
    "HaidRaw_sigmaGas",
    "HaidRaw_nepGas",
    "HaidRaw_thetaLiq",
    "HaidRaw_sigmaLiq",
    "HaidRaw_nepLiq",
}
EXPECTED_CONTINUOUS_BRANCHES = {
    "HaidRaw_thetaGas",
    "HaidRaw_sigmaGas",
    "HaidRaw_nepGas",
    "HaidRaw_thetaLiq",
    "HaidRaw_sigmaLiq",
    "HaidRaw_nepLiq",
}
REFERENCE_BRANCHES = {
    "HaidRaw_richQuality": "Haid_richQuality",
    "HaidRaw_thetaGas": "Rich_theg",
    "HaidRaw_sigmaGas": "Rich_sigg",
    "HaidRaw_nphGas": "Rich_nphg",
    "HaidRaw_nepGas": "Rich_nepg",
    "HaidRaw_flagGas": "Rich_flagg",
    "HaidRaw_thetaLiq": "Rich_thel",
    "HaidRaw_sigmaLiq": "Rich_sigl",
    "HaidRaw_nphLiq": "Rich_nphl",
    "HaidRaw_nepLiq": "Rich_nepl",
    "HaidRaw_flagLiq": "Rich_flagl",
}
RADIATORS = {
    "gas": {
        "theta": "HaidRaw_thetaGas",
        "sigma": "HaidRaw_sigmaGas",
        "nph": "HaidRaw_nphGas",
        "nep": "HaidRaw_nepGas",
        "flag": "HaidRaw_flagGas",
    },
    "liquid": {
        "theta": "HaidRaw_thetaLiq",
        "sigma": "HaidRaw_sigmaLiq",
        "nph": "HaidRaw_nphLiq",
        "nep": "HaidRaw_nepLiq",
        "flag": "HaidRaw_flagLiq",
    },
}


def flatten(array: ak.Array) -> np.ndarray:
    return np.asarray(ak.to_numpy(ak.flatten(array, axis=None)))


def ratio(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def branch_status(branch: str, dtype: str, stats: dict[str, float | int | str]) -> str:
    if branch == "HaidRaw_richQuality" and "int8" in dtype:
        return "invalid_narrow_dtype"
    if float(stats["subnormal_fraction"]) > 0.001:
        return "invalid_subnormal"
    if branch in EXPECTED_CONTINUOUS_BRANCHES and float(stats["integer_like_fraction"]) > 0.999:
        return "review_integer_like"
    if branch.endswith(("flagGas", "flagLiq")) and int(stats["n_nonzero"]) == 0:
        return "review_all_zero"
    return "ok"


def summarize_branch(sample: str, branch: str, dtype: str, values: np.ndarray) -> dict[str, float | int | str]:
    finite = values[np.isfinite(values)]
    nonzero = finite[finite != 0]
    is_float = branch in FLOAT_BRANCHES
    if is_float:
        subnormal = (np.abs(nonzero) < np.finfo(np.float32).tiny)
        integer_like = np.isclose(nonzero, np.rint(nonzero), rtol=0.0, atol=1e-6)
    else:
        subnormal = np.zeros(len(nonzero), dtype=bool)
        integer_like = np.ones(len(nonzero), dtype=bool)

    quantiles = np.quantile(finite, [0.01, 0.5, 0.99]) if len(finite) else [np.nan] * 3
    stats: dict[str, float | int | str] = {
        "sample": sample,
        "branch": branch,
        "dtype": dtype,
        "n_value": int(len(values)),
        "n_finite": int(len(finite)),
        "n_nonzero": int(len(nonzero)),
        "zero_fraction": ratio(int(np.count_nonzero(finite == 0)), len(finite)),
        "subnormal_fraction": ratio(int(np.count_nonzero(subnormal)), len(nonzero)),
        "integer_like_fraction": ratio(int(np.count_nonzero(integer_like)), len(nonzero)),
        "n_unique": int(len(np.unique(finite))),
        "min": float(np.min(finite)) if len(finite) else float("nan"),
        "q01": float(quantiles[0]),
        "median": float(quantiles[1]),
        "q99": float(quantiles[2]),
        "max": float(np.max(finite)) if len(finite) else float("nan"),
    }
    stats["status"] = branch_status(branch, dtype, stats)
    return stats


def histogram_edges(values_by_sample: dict[str, np.ndarray], integer: bool) -> np.ndarray:
    finite = np.concatenate(
        [values[np.isfinite(values)] for values in values_by_sample.values()]
    )
    if integer:
        lo = int(np.min(finite))
        hi = int(np.max(finite))
        if hi - lo <= 100:
            return np.arange(lo - 0.5, hi + 1.5, 1.0)

    lo, hi = np.quantile(finite, [0.001, 0.999])
    if lo == hi:
        lo = float(np.min(finite))
        hi = float(np.max(finite))
    if lo == hi:
        pad = abs(float(lo)) * 0.1 if lo else 0.5
    else:
        pad = float(hi - lo) * 0.05
    return np.linspace(float(lo) - pad, float(hi) + pad, 101)


def plot_distribution(
    output_dir: Path,
    branch: str,
    values_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
    integer: bool,
) -> None:
    bins = histogram_edges(values_by_sample, integer)
    fig, ax = plt.subplots(figsize=(12, 10))
    counts_list = []
    fractions_list = []

    for sample in samples:
        values = values_by_sample[sample]
        values = values[np.isfinite(values)]
        counts, _ = np.histogram(values, bins=bins)
        fractions = counts / len(values)
        counts_list.append(counts)
        fractions_list.append(fractions)
        color, hatch = styles[sample]
        weights = np.full(len(values), 1.0 / len(values))
        ax.hist(
            values,
            bins=bins,
            weights=weights,
            histtype="stepfilled",
            facecolor="none",
            edgecolor=color,
            hatch=hatch,
            linewidth=0.0,
            label=f"{sample:<8} N={len(values)}",
        )
        ax.hist(
            values,
            bins=bins,
            weights=weights,
            histtype="step",
            color=color,
            linewidth=2,
            label="_nolegend_",
        )

    set_hist_yaxis(ax, counts_list, fractions_list)
    if not integer and max(abs(bins[0]), abs(bins[-1])) < 1e-3:
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    ax.set_xlabel(branch)
    ax.set_ylabel("Fraction / bin")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 15})
    ax.grid(alpha=0.25)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{branch}.png", dpi=150)
    plt.close(fig)


def plot_raw_bits(
    output_dir: Path,
    branch: str,
    values_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    bits_by_sample = {
        sample: np.asarray(values[np.isfinite(values)], dtype=np.float32).view(np.uint32)
        for sample, values in values_by_sample.items()
    }
    max_bit = max(int(np.max(values)) for values in bits_by_sample.values())
    bins = np.arange(-0.5, max_bit + 1.5, 1.0) if max_bit <= 200 else np.linspace(0, max_bit, 101)
    fig, ax = plt.subplots(figsize=(12, 10))
    for sample in samples:
        values = bits_by_sample[sample]
        weights = np.full(len(values), 1.0 / len(values))
        color, _ = styles[sample]
        ax.hist(values, bins=bins, weights=weights, histtype="step", linewidth=2, color=color, label=sample)
    ax.set_xlabel(f"{branch} raw float32 bits interpreted as uint32")
    ax.set_ylabel("Fraction / bin")
    ax.legend()
    ax.grid(alpha=0.25)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{branch}_raw_bits.png", dpi=150)
    plt.close(fig)


def relation_range(values_by_sample: dict[str, np.ndarray]) -> tuple[float, float]:
    values = np.concatenate([value[np.isfinite(value)] for value in values_by_sample.values()])
    lo, hi = np.quantile(values, [0.0, 0.995])
    if lo == hi:
        hi = float(np.max(values))
    if lo == hi:
        pad = abs(float(lo)) * 0.1 if lo else 0.5
        return float(lo) - pad, float(hi) + pad
    return float(lo), float(hi)


def plot_relation(
    output_dir: Path,
    name: str,
    x_branch: str,
    y_branch: str,
    sample_values: dict[str, dict[str, np.ndarray]],
    samples: tuple[str, ...],
) -> None:
    x_by_sample = {sample: sample_values[sample][x_branch] for sample in samples}
    y_by_sample = {sample: sample_values[sample][y_branch] for sample in samples}
    x_range = relation_range(x_by_sample)
    y_range = relation_range(y_by_sample)
    fig, axes = sample_grid(samples, subplot_width=7.0, subplot_height=6.0)

    for ax, sample in zip(axes, samples, strict=True):
        x = x_by_sample[sample]
        y = y_by_sample[sample]
        valid = np.isfinite(x) & np.isfinite(y) & ((x != 0) | (y != 0))
        hist, x_edges, y_edges = np.histogram2d(
            x[valid], y[valid], bins=80, range=(x_range, y_range)
        )
        total = np.sum(hist)
        values = hist / total if total else hist
        mesh = ax.pcolormesh(x_edges, y_edges, np.ma.masked_where(values.T == 0, values.T), cmap="viridis")
        fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.02, label="Fraction")
        add_heatmap_label(ax, f"Sample: {sample}\nPairs: {int(np.count_nonzero(valid))}", family="monospace")
        ax.set_xlabel(x_branch)
        ax.set_ylabel(y_branch)
        ax.grid(alpha=0.2)
        if max(abs(x_range[0]), abs(x_range[1])) < 1e-3:
            ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        if max(abs(y_range[0]), abs(y_range[1])) < 1e-3:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def summarize_consistency(
    sample: str,
    radiator: str,
    values: dict[str, np.ndarray],
) -> dict[str, float | int | str]:
    branches = RADIATORS[radiator]
    theta = values[branches["theta"]]
    sigma = values[branches["sigma"]]
    nph = values[branches["nph"]]
    nep = values[branches["nep"]]
    flag = values[branches["flag"]]
    normal_min = np.finfo(np.float32).tiny
    has_angle = np.isfinite(theta) & (theta > 0)
    valid_sigma = np.isfinite(sigma) & (sigma >= normal_min)
    has_photons = np.isfinite(nph) & (nph > 0)
    valid_expected = np.isfinite(nep) & (nep >= normal_min)
    return {
        "sample": sample,
        "radiator": radiator,
        "n_object": int(len(theta)),
        "n_angle_positive": int(np.count_nonzero(has_angle)),
        "valid_sigma_given_angle": ratio(int(np.count_nonzero(has_angle & valid_sigma)), int(np.count_nonzero(has_angle))),
        "n_photon_positive": int(np.count_nonzero(has_photons)),
        "valid_expected_given_photons": ratio(int(np.count_nonzero(has_photons & valid_expected)), int(np.count_nonzero(has_photons))),
        "flag_nonzero_fraction": ratio(int(np.count_nonzero(flag)), len(flag)),
    }


def plot_consistency(
    output_dir: Path,
    rows: list[dict[str, float | int | str]],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    labels = ("gas sigma", "gas expected", "liquid sigma", "liquid expected")
    fig, ax = plt.subplots(figsize=(12, 8))
    x = np.arange(len(labels))
    width = 0.8 / len(samples)
    for index, sample in enumerate(samples):
        sample_rows = {row["radiator"]: row for row in rows if row["sample"] == sample}
        values = (
            sample_rows["gas"]["valid_sigma_given_angle"],
            sample_rows["gas"]["valid_expected_given_photons"],
            sample_rows["liquid"]["valid_sigma_given_angle"],
            sample_rows["liquid"]["valid_expected_given_photons"],
        )
        color, hatch = styles[sample]
        ax.bar(x + (index - (len(samples) - 1) / 2) * width, values, width, label=sample, color=color, hatch=hatch)
    ax.set_xticks(x, labels)
    ax.set_ylim(0.0, 1.1)
    ax.set_ylabel("Conditional valid fraction")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / "consistency_fractions.png", dpi=150)
    plt.close(fig)


def collection_summary(
    sample: str,
    arrays: dict[str, ak.Array],
) -> dict[str, int | str]:
    pa_idx = ak.to_list(arrays["HaidRaw_paIdx"])
    track_pa_idx = ak.to_list(arrays["TracRaw_paIdx"])
    counts = np.asarray(ak.to_numpy(ak.num(arrays["HaidRaw_paIdx"], axis=1)))
    reported = np.asarray(ak.to_numpy(arrays["nHaidRaw"]))
    length_mismatch = np.zeros(len(counts), dtype=bool)
    for branch in RICH_BRANCHES:
        length_mismatch |= np.asarray(ak.to_numpy(ak.num(arrays[branch], axis=1))) != counts
    return {
        "sample": sample,
        "n_event": int(len(counts)),
        "n_object": int(np.sum(counts)),
        "reported_count_mismatch_events": int(np.count_nonzero(counts != reported)),
        "vector_length_mismatch_events": int(np.count_nonzero(length_mismatch)),
        "duplicate_paidx_events": int(sum(len(ids) != len(set(ids)) for ids in pa_idx)),
        "direct_track_index_mismatch_events": int(
            sum(not set(track_ids).issubset(set(ids)) for ids, track_ids in zip(pa_idx, track_pa_idx, strict=True))
        ),
        "one_based_track_index_mismatch_events": int(
            sum(
                not {int(track_id) + 1 for track_id in track_ids}.issubset(set(ids))
                for ids, track_ids in zip(pa_idx, track_pa_idx, strict=True)
            )
        ),
    }


def compare_reference(
    sample: str,
    source_tree,
    reference_tree,
) -> tuple[list[dict[str, float | int | str]], tuple[np.ndarray, np.ndarray] | None]:
    rows = []
    quality_pair = None
    source_types = source_tree.typenames()
    reference_types = reference_tree.typenames()
    for source_branch, reference_branch in REFERENCE_BRANCHES.items():
        if source_branch not in source_tree or reference_branch not in reference_tree:
            continue
        source = flatten(source_tree[source_branch].array(library="ak"))
        reference = flatten(reference_tree[reference_branch].array(library="ak"))
        same_size = len(source) == len(reference)
        if same_size:
            equal = np.isclose(source, reference, rtol=0.0, atol=0.0, equal_nan=True)
            difference = np.abs(source.astype(np.float64) - reference.astype(np.float64))
        else:
            equal = np.zeros(0, dtype=bool)
            difference = np.zeros(0, dtype=float)
        rows.append({
            "sample": sample,
            "source_branch": source_branch,
            "reference_branch": reference_branch,
            "source_dtype": source_types[source_branch],
            "reference_dtype": reference_types[reference_branch],
            "source_count": int(len(source)),
            "reference_count": int(len(reference)),
            "equal_fraction": ratio(int(np.count_nonzero(equal)), len(equal)),
            "max_abs_difference": float(np.max(difference)) if len(difference) else float("nan"),
        })
        if source_branch == "HaidRaw_richQuality" and same_size:
            quality_pair = source, reference
    return rows, quality_pair


def plot_quality_mapping(
    output_dir: Path,
    quality_by_sample: dict[str, tuple[np.ndarray, np.ndarray]],
    samples: tuple[str, ...],
) -> None:
    fig, axes = sample_grid(samples, subplot_width=7.0, subplot_height=6.0)
    for ax, sample in zip(axes, samples, strict=True):
        source, reference = quality_by_sample[sample]
        step = max(len(source) // 30000, 1)
        ax.scatter(reference[::step], source[::step], s=3, alpha=0.15)
        ax.set_xlabel("Haid_richQuality [int32]")
        ax.set_ylabel("HaidRaw_richQuality [int8]")
        add_heatmap_label(ax, f"Sample: {sample}", family="monospace")
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "rich_quality_source_vs_reference.png", dpi=150)
    plt.close(fig)


def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_rich(
    input_root: Path,
    output_root: Path,
    data_root: Path,
    samples: list[str] | tuple[str, ...] | None = None,
) -> None:
    mh.style.use(mh.styles.CMS)
    samples = resolve_samples(input_root, samples)
    styles = sample_styles(samples)
    output_dir = output_root / "rich"
    distribution_dir = output_dir / "distributions"
    bit_dir = output_dir / "raw_bits"
    relation_dir = output_dir / "relations"
    distribution_dir.mkdir(parents=True, exist_ok=True)
    bit_dir.mkdir(parents=True, exist_ok=True)
    relation_dir.mkdir(parents=True, exist_ok=True)

    required = ("nHaidRaw", "HaidRaw_paIdx", "TracRaw_paIdx", *RICH_BRANCHES)
    sample_values: dict[str, dict[str, np.ndarray]] = {}
    branch_rows = []
    consistency_rows = []
    collection_rows = []
    comparison_rows = []
    quality_by_sample = {}

    for sample in samples:
        source_path = input_root / sample / SOURCE_FILE
        reference_path = input_root / sample / REFERENCE_FILE
        with uproot.open(source_path) as source_file:
            source_tree = source_file[TREE_NAME]
            missing = [branch for branch in required if branch not in source_tree]
            if missing:
                raise KeyError(f"Missing RICH study branches in {source_path}: {missing}")
            arrays = {branch: source_tree[branch].array(library="ak") for branch in required}
            values = {branch: flatten(arrays[branch]) for branch in RICH_BRANCHES}
            sample_values[sample] = values
            types = source_tree.typenames()
            branch_rows.extend(
                summarize_branch(sample, branch, types[branch], values[branch])
                for branch in RICH_BRANCHES
            )
            consistency_rows.extend(
                summarize_consistency(sample, radiator, values)
                for radiator in RADIATORS
            )
            collection_rows.append(collection_summary(sample, arrays))

            if reference_path.exists():
                with uproot.open(reference_path) as reference_file:
                    rows, quality_pair = compare_reference(sample, source_tree, reference_file[TREE_NAME])
                    comparison_rows.extend(rows)
                    if quality_pair is not None:
                        quality_by_sample[sample] = quality_pair

    for branch in RICH_BRANCHES:
        values_by_sample = {sample: sample_values[sample][branch] for sample in samples}
        plot_distribution(
            distribution_dir,
            branch,
            values_by_sample,
            samples,
            styles,
            branch not in FLOAT_BRANCHES,
        )
        branch_stats = [row for row in branch_rows if row["branch"] == branch]
        if branch in FLOAT_BRANCHES and any(float(row["subnormal_fraction"]) > 0.001 for row in branch_stats):
            plot_raw_bits(bit_dir, branch, values_by_sample, samples, styles)

    for radiator, branches in RADIATORS.items():
        plot_relation(
            relation_dir,
            f"{radiator}_theta_vs_sigma",
            branches["theta"],
            branches["sigma"],
            sample_values,
            samples,
        )
        plot_relation(
            relation_dir,
            f"{radiator}_observed_vs_expected_photons",
            branches["nph"],
            branches["nep"],
            sample_values,
            samples,
        )

    plot_consistency(output_dir, consistency_rows, samples, styles)
    if set(quality_by_sample) == set(samples):
        plot_quality_mapping(output_dir, quality_by_sample, samples)

    write_csv(data_root / "rich" / "branch_summary.csv", branch_rows)
    write_csv(data_root / "rich" / "consistency_summary.csv", consistency_rows)
    write_csv(data_root / "rich" / "collection_summary.csv", collection_rows)
    if comparison_rows:
        write_csv(data_root / "rich" / "source_comparison.csv", comparison_rows)

    suspect = sorted({str(row["branch"]) for row in branch_rows if row["status"] != "ok"})
    print(f"RICH plots: {output_dir}")
    print(f"RICH diagnostics: {data_root / 'rich'}")
    print(f"branches requiring review: {', '.join(suspect) if suspect else 'none'}")
