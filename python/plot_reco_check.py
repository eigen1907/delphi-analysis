from __future__ import annotations

import shutil
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from plot_utils import (
    MATCH_P_BINS,
    MATCH_PHI_BINS,
    MATCH_PT_BINS,
    MATCH_THETA_BINS,
    add_delphi_label,
    phi_0_2pi,
    resolve_samples,
    sample_grid,
    sample_styles,
    set_hist_yaxis,
)


FILE_NAME = "nanoaod_raw_sdst.root"
TREE_NAME = "Events"
TRACK_NAME = "TrackRaw"
TRACK_BRANCHES = (
    "Event_bFieldGevCm",
    "TracRaw_theta",
    "TracRaw_phi",
    "TracRaw_invR",
    "TracRaw_charge",
)
HIST_PLOTS = (
    ("inv_r", r"TrackRaw $1/R$", (-0.05, 0.05)),
    ("pt", r"TrackRaw $p_T$ [GeV]", (0.0, 60.0)),
    ("p", r"TrackRaw $p$ [GeV]", (0.0, 60.0)),
    ("theta", r"TrackRaw $\theta$ [rad]", (0.0, np.pi)),
    ("phi", r"TrackRaw $\phi$ [rad]", (0.0, 2.0 * np.pi)),
    ("px", r"TrackRaw $p_x$ [GeV]", (-60.0, 60.0)),
    ("py", r"TrackRaw $p_y$ [GeV]", (-60.0, 60.0)),
    ("pz", r"TrackRaw $p_z$ [GeV]", (-120.0, 120.0)),
)
PLOTS_2D = (
    ("pt_theta", "pt_theta_pt", "pt_theta_theta", r"$p_T$ [GeV]", r"$\theta$ [rad]", MATCH_PT_BINS, MATCH_THETA_BINS),
    ("p_theta", "p_theta_p", "p_theta_theta", r"$p$ [GeV]", r"$\theta$ [rad]", MATCH_P_BINS, MATCH_THETA_BINS),
    ("theta_phi", "theta_phi_theta", "theta_phi_phi", r"$\theta$ [rad]", r"$\phi$ [rad]", MATCH_THETA_BINS, MATCH_PHI_BINS),
)
LEGEND_HEADER = f"{'Sample':<8} {'Event':<8} {'Object':<8} {'UF':<6} {'OF':<6}"


def flatten(array) -> np.ndarray:
    values = np.asarray(ak.to_numpy(ak.flatten(array, axis=None)), dtype=float)
    return values[np.isfinite(values)]


def finite(values: np.ndarray) -> np.ndarray:
    return values[np.isfinite(values)]


def plot_hist(
    output_dir: Path,
    name: str,
    xlabel: str,
    value_range: tuple[float, float],
    values_by_sample: dict[str, np.ndarray],
    n_events: dict[str, int],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    bins = np.linspace(value_range[0], value_range[1], 101)

    ax.plot([], [], color="none", label=LEGEND_HEADER)
    counts_list = []
    density_list = []
    for sample in samples:
        values = values_by_sample[sample]
        color, hatch = styles[sample]
        underflow = int(np.count_nonzero(values < bins[0]))
        overflow = int(np.count_nonzero(values > bins[-1]))
        counts, _ = np.histogram(values, bins=bins)
        counts_list.append(counts)
        density_list.append(counts / (np.sum(counts) * np.diff(bins)) if np.sum(counts) else counts)
        label = f"{sample:<8} {n_events[sample]:<8} {len(values):<8} {underflow:<6} {overflow:<6}"
        if len(values):
            ax.hist(values, bins=bins, histtype="stepfilled", density=True, label=label, facecolor="none", edgecolor=color, hatch=hatch, linewidth=0.0)
            ax.hist(values, bins=bins, histtype="step", density=True, label="_nolegend_", color=color, alpha=0.9, linewidth=2)
        else:
            ax.plot([], [], color=color, label=label)

    set_hist_yaxis(ax, counts_list, density_list)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalized")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 15})
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_discrete(
    output_dir: Path,
    name: str,
    xlabel: str,
    values_by_sample: dict[str, np.ndarray],
    n_events: dict[str, int],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    unique_values = sorted({int(value) for values in values_by_sample.values() for value in values})
    if not unique_values:
        unique_values = [0]
    x = np.arange(len(unique_values))

    ax.plot([], [], color="none", label=LEGEND_HEADER)
    counts_list = []
    fraction_list = []
    for sample in samples:
        values = values_by_sample[sample].astype(int)
        raw_counts = np.asarray([np.count_nonzero(values == value) for value in unique_values], dtype=float)
        n_object = int(np.sum(raw_counts))
        counts = raw_counts / n_object if n_object else raw_counts
        counts_list.append(raw_counts)
        fraction_list.append(counts)
        color, hatch = styles[sample]
        label = f"{sample:<8} {n_events[sample]:<8} {n_object:<8} {0:<6} {0:<6}"
        ax.bar(x, counts, width=0.8, label=label, facecolor="none", edgecolor=color, hatch=hatch, alpha=0.9)

    set_hist_yaxis(ax, counts_list, fraction_list)
    ax.set_xticks(x)
    ax.set_xticklabels([str(value) for value in unique_values])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalized")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 15})
    ax.grid(axis="y", alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_int_hist(
    output_dir: Path,
    name: str,
    xlabel: str,
    values_by_sample: dict[str, np.ndarray],
    n_events: dict[str, int],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    all_values = np.concatenate([values for values in values_by_sample.values() if len(values)])
    lo = int(np.min(all_values)) if len(all_values) else 0
    hi = int(np.max(all_values)) if len(all_values) else 0
    bins = np.arange(lo - 0.5, hi + 1.5, 1.0)

    ax.plot([], [], color="none", label=LEGEND_HEADER)
    counts_list = []
    density_list = []
    for sample in samples:
        values = values_by_sample[sample]
        color, hatch = styles[sample]
        underflow = int(np.count_nonzero(values < bins[0]))
        overflow = int(np.count_nonzero(values > bins[-1]))
        counts, _ = np.histogram(values, bins=bins)
        counts_list.append(counts)
        density_list.append(counts / (np.sum(counts) * np.diff(bins)) if np.sum(counts) else counts)
        label = f"{sample:<8} {n_events[sample]:<8} {len(values):<8} {underflow:<6} {overflow:<6}"
        ax.hist(values, bins=bins, histtype="stepfilled", density=True, label=label, facecolor="none", edgecolor=color, hatch=hatch, linewidth=0.0)
        ax.hist(values, bins=bins, histtype="step", density=True, label="_nolegend_", color=color, alpha=0.9, linewidth=2)

    set_hist_yaxis(ax, counts_list, density_list)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalized")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 15})
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_2d(
    output_dir: Path,
    name: str,
    x_label: str,
    y_label: str,
    x_bins: np.ndarray,
    y_bins: np.ndarray,
    x_values_by_sample: dict[str, np.ndarray],
    y_values_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
) -> None:
    counts_by_sample = {}
    for sample in samples:
        counts, _, _ = np.histogram2d(x_values_by_sample[sample], y_values_by_sample[sample], bins=(x_bins, y_bins))
        counts_by_sample[sample] = counts

    vmax = max(1.0, max(float(np.max(counts)) for counts in counts_by_sample.values()))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="white", alpha=0.0)
    fig, axes = sample_grid(samples)
    mesh = None

    for ax, sample in zip(axes, samples, strict=True):
        counts = counts_by_sample[sample]
        values = np.ma.masked_where(counts <= 0, counts)
        mesh = ax.pcolormesh(x_bins, y_bins, values.T, vmin=0.0, vmax=vmax, cmap=cmap)
        ax.text(
            0.94,
            0.92,
            f"{sample}\nObject {int(np.sum(counts))}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=14,
            family="monospace",
        )
        ax.grid(alpha=0.2)

    for ax in axes:
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    fig.subplots_adjust(right=0.86, hspace=0.08, wspace=0.08)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.025, 0.70])
    cbar = fig.colorbar(mesh, cax=cbar_ax)
    cbar.set_label("Entries")
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_reco_check(input_root: Path, output_root: Path, samples: list[str] | tuple[str, ...] | None = None) -> None:
    mh.style.use(mh.styles.CMS)
    samples = resolve_samples(input_root, samples)
    styles = sample_styles(samples)
    output_dir = output_root / "reco_check" / TRACK_NAME
    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_events = {}
    values = {
        "n_track_raw": {},
        "charge": {},
        "inv_r": {},
        "pt": {},
        "p": {},
        "theta": {},
        "phi": {},
        "px": {},
        "py": {},
        "pz": {},
        "pt_theta_pt": {},
        "pt_theta_theta": {},
        "p_theta_p": {},
        "p_theta_theta": {},
        "theta_phi_theta": {},
        "theta_phi_phi": {},
    }

    for sample in samples:
        tree = uproot.open(input_root / sample / FILE_NAME)[TREE_NAME]
        arrays = {branch: tree[branch].array(library="ak") for branch in TRACK_BRANCHES}
        n_events[sample] = int(tree.num_entries)

        theta = flatten(arrays["TracRaw_theta"])
        phi = phi_0_2pi(flatten(arrays["TracRaw_phi"]))
        inv_r = flatten(arrays["TracRaw_invR"])
        charge = flatten(arrays["TracRaw_charge"])
        bfield_gev_cm = flatten(ak.broadcast_arrays(arrays["Event_bFieldGevCm"], arrays["TracRaw_invR"])[0])

        has_pt = inv_r != 0
        pt = bfield_gev_cm[has_pt] / np.abs(inv_r[has_pt])
        theta_with_pt = theta[has_pt]
        phi_with_pt = phi[has_pt]
        p = pt / np.sin(theta_with_pt)
        px = pt * np.cos(phi_with_pt)
        py = pt * np.sin(phi_with_pt)
        pz = pt / np.tan(theta_with_pt)

        values["n_track_raw"][sample] = np.asarray(ak.to_numpy(ak.num(arrays["TracRaw_theta"], axis=1)), dtype=float)
        values["charge"][sample] = charge
        values["inv_r"][sample] = inv_r
        values["pt"][sample] = finite(pt)
        values["p"][sample] = finite(p)
        values["theta"][sample] = theta
        values["phi"][sample] = phi
        values["px"][sample] = finite(px)
        values["py"][sample] = finite(py)
        values["pz"][sample] = finite(pz)
        values["pt_theta_pt"][sample] = finite(pt)
        values["pt_theta_theta"][sample] = theta_with_pt[np.isfinite(pt)]
        values["p_theta_p"][sample] = finite(p)
        values["p_theta_theta"][sample] = theta_with_pt[np.isfinite(p)]
        values["theta_phi_theta"][sample] = theta
        values["theta_phi_phi"][sample] = phi

    plot_int_hist(output_dir, "n_track_raw", r"$N_{\mathrm{TrackRaw}}$", values["n_track_raw"], n_events, samples, styles)
    plot_discrete(output_dir, "charge", "TrackRaw charge", values["charge"], n_events, samples, styles)
    for plot_name, xlabel, value_range in HIST_PLOTS:
        plot_hist(output_dir, plot_name, xlabel, value_range, values[plot_name], n_events, samples, styles)
    for plot_name, x_values, y_values, x_label, y_label, x_bins, y_bins in PLOTS_2D:
        plot_2d(output_dir, plot_name, x_label, y_label, x_bins, y_bins, values[x_values], values[y_values], samples)

    print(f"wrote {output_dir}")
