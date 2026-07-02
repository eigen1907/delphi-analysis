from __future__ import annotations

from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from plot_utils import (
    CHARGED_ABS_PDG,
    GEN_P4_BRANCHES,
    MATCH_ALPHA_CUTS,
    MATCH_AXES,
    MATCH_SOURCES,
    charge_from_pdg,
    divide,
    has_ancestor_pdg,
    infer_expected_abs_pdg,
    is_track_fiducial,
    match_cut_dir,
    opening_angle,
    phi_0_2pi,
    resolve_samples,
    sample_grid,
    sample_styles,
    set_hist_yaxis,
)


GEN_AXIS_LABELS = {
    "pt": r"Stable charged GenPart $p_T$ [GeV]",
    "p": r"Stable charged GenPart $p$ [GeV]",
    "theta": r"Stable charged GenPart $\theta$ [rad]",
    "phi": r"Stable charged GenPart $\phi$ [rad]",
}
RECO_AXIS_LABELS = {
    "pt": r"Reco track $p_T$ [GeV]",
    "p": r"Reco track $p$ [GeV]",
    "theta": r"Reco track $\theta$ [rad]",
    "phi": r"Reco track $\phi$ [rad]",
}
AXIS_BINS = {axis: bins for axis, bins, _ in MATCH_AXES}
PLAIN_AXIS_LABELS = {axis: label for axis, _, label in MATCH_AXES}
GEN_2D_AXES = (
    ("pt_theta", "pt", "theta"),
    ("p_theta", "p", "theta"),
    ("theta_phi", "theta", "phi"),
)
RECO_2D_AXES = GEN_2D_AXES


def fill_1d(hist: np.ndarray, value: float, bins: np.ndarray) -> None:
    index = np.searchsorted(bins, value, side="right") - 1
    if value == bins[-1]:
        index = len(hist) - 1
    if 0 <= index < len(hist):
        hist[index] += 1


def fill_2d(hist: np.ndarray, x_value: float, y_value: float, x_bins: np.ndarray, y_bins: np.ndarray) -> None:
    x_index = np.searchsorted(x_bins, x_value, side="right") - 1
    y_index = np.searchsorted(y_bins, y_value, side="right") - 1
    if x_value == x_bins[-1]:
        x_index = hist.shape[0] - 1
    if y_value == y_bins[-1]:
        y_index = hist.shape[1] - 1
    if 0 <= x_index < hist.shape[0] and 0 <= y_index < hist.shape[1]:
        hist[x_index, y_index] += 1


def plot_scan(
    output_dir: Path,
    name: str,
    ylabel: str,
    values_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.axhline(1.0, color="black", linestyle=":", linewidth=2, alpha=0.8)

    for sample in samples:
        color, _ = styles[sample]
        ax.plot(MATCH_ALPHA_CUTS, values_by_sample[sample], label=sample, color=color, linewidth=3, alpha=0.7)

    ax.set_xlim(MATCH_ALPHA_CUTS[0], MATCH_ALPHA_CUTS[-1])
    ax.set_ylim(0.0, 1.25)
    ax.set_xlabel(r"$\alpha$ cut [rad]")
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper right", fontsize=16)
    ax.grid(alpha=0.8)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}_vs_alpha.png", dpi=150)
    plt.close(fig)


def plot_metric(
    output_dir: Path,
    name: str,
    xlabel: str,
    ylabel: str,
    bins: np.ndarray,
    numerator_by_sample: dict[str, np.ndarray],
    denominator_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    centers = 0.5 * (bins[:-1] + bins[1:])
    ax.axhline(1.0, color="black", linestyle=":", linewidth=2, alpha=0.8)

    for sample in samples:
        color, _ = styles[sample]
        values = divide(numerator_by_sample[sample], denominator_by_sample[sample])
        ax.step(centers, values, where="mid", label=sample, color=color, linewidth=3, alpha=0.7)
        #ax.scatter(centers, values, color=color, s=25, alpha=0.7)

    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0.0, 1.25)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper right", fontsize=16)
    ax.grid(alpha=0.8)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_metric_2d(
    output_dir: Path,
    name: str,
    x_label: str,
    y_label: str,
    x_bins: np.ndarray,
    y_bins: np.ndarray,
    values_label: str,
    numerator_by_sample: dict[str, np.ndarray],
    denominator_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
) -> None:
    fig, axes = sample_grid(samples)
    mesh = None

    for ax, sample in zip(axes, samples, strict=True):
        values = divide(numerator_by_sample[sample], denominator_by_sample[sample])
        mesh = ax.pcolormesh(x_bins, y_bins, values.T, vmin=0.0, vmax=1.0, cmap="viridis")
        ax.text(0.94, 0.92, sample, transform=ax.transAxes, ha="right", fontsize=22)
        ax.grid(alpha=0.2)

    for ax in axes:
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    fig.subplots_adjust(right=0.86, hspace=0.08, wspace=0.08)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.025, 0.70])
    cbar = fig.colorbar(mesh, cax=cbar_ax)
    cbar.set_label(values_label)
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_count(
    output_dir: Path,
    name: str,
    xlabel: str,
    bins: np.ndarray,
    counts_by_sample: dict[str, np.ndarray],
    matched_by_sample: dict[str, np.ndarray],
    missed_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    centers = 0.5 * (bins[:-1] + bins[1:])
    counts_list = []

    ax.plot([], [], color="none", label=f"{'Sample':<8} {'Object':<8} {'Matched':<8} {'Missed':<8}")
    for sample in samples:
        counts = counts_by_sample[sample]
        matched = int(np.sum(matched_by_sample[sample]))
        missed = int(np.sum(missed_by_sample[sample]))
        n_object = matched + missed
        counts_list.append(counts)
        color, _ = styles[sample]
        label = f"{sample:<8} {n_object:<8} {matched:<8} {missed:<8}"
        ax.step(centers, counts, where="mid", label=label, color=color, linewidth=3, alpha=0.7)

    set_hist_yaxis(ax, counts_list)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entries")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 15})
    ax.grid(alpha=0.8)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_count_2d(
    output_dir: Path,
    name: str,
    x_label: str,
    y_label: str,
    x_bins: np.ndarray,
    y_bins: np.ndarray,
    counts_by_sample: dict[str, np.ndarray],
    matched_by_sample: dict[str, np.ndarray],
    missed_by_sample: dict[str, np.ndarray],
    samples: tuple[str, ...],
) -> None:
    vmax = max(1.0, max((float(np.max(counts)) for counts in counts_by_sample.values() if counts.size), default=1.0))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="white", alpha=0.0)
    fig, axes = sample_grid(samples)
    mesh = None

    for ax, sample in zip(axes, samples, strict=True):
        counts = counts_by_sample[sample]
        matched = int(np.sum(matched_by_sample[sample]))
        missed = int(np.sum(missed_by_sample[sample]))
        n_object = matched + missed
        values = np.ma.masked_where(counts <= 0, counts)
        mesh = ax.pcolormesh(x_bins, y_bins, values.T, vmin=0.0, vmax=vmax, cmap=cmap)
        label = f"{sample}\n{'Object':>7} {'Matched':>7} {'Missed':>7}\n{n_object:>7} {matched:>7} {missed:>7}"
        ax.text(
            0.94,
            0.92,
            label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
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


def plot_gen_reco_track_match_cut(
    input_root: Path,
    output_root: Path,
    alpha_cut: float = 0.05,
    fiducial_cos_max: float = 1.00,
    pt_min: float = 0.00,
    samples: list[str] | tuple[str, ...] | None = None,
) -> None:
    mh.style.use(mh.styles.CMS)
    samples = resolve_samples(input_root, samples)
    styles = sample_styles(samples)
    output_root = output_root / "gen_reco_track_match_cut"
    scan_output = output_root / "alpha_scan"
    cut_output = output_root / match_cut_dir(alpha_cut, fiducial_cos_max, pt_min)
    scan_output.mkdir(parents=True, exist_ok=True)
    cut_output.mkdir(parents=True, exist_ok=True)

    scan_values = {
        "efficiency": {},
        "fake_rate": {},
        "duplicate_rate": {},
        "ancestor_not_gamma_status_1_expected_efficiency": {},
        "ancestor_not_gamma_status_1_expected_pair_efficiency": {},
    }
    metric_num = {
        "efficiency": {axis: {} for axis, _, _ in MATCH_AXES},
        "fake_rate": {axis: {} for axis, _, _ in MATCH_AXES},
        "duplicate_rate": {axis: {} for axis, _, _ in MATCH_AXES},
        "ancestor_not_gamma_status_1_expected_efficiency": {axis: {} for axis, _, _ in MATCH_AXES},
    }
    metric_den = {
        "efficiency": {axis: {} for axis, _, _ in MATCH_AXES},
        "fake_rate": {axis: {} for axis, _, _ in MATCH_AXES},
        "duplicate_rate": {axis: {} for axis, _, _ in MATCH_AXES},
        "ancestor_not_gamma_status_1_expected_efficiency": {axis: {} for axis, _, _ in MATCH_AXES},
    }
    metric_num_2d = {
        "efficiency": {name: {} for name, _, _ in GEN_2D_AXES},
        "fake_rate": {name: {} for name, _, _ in RECO_2D_AXES},
        "duplicate_rate": {name: {} for name, _, _ in GEN_2D_AXES},
        "ancestor_not_gamma_status_1_expected_efficiency": {name: {} for name, _, _ in GEN_2D_AXES},
    }
    metric_den_2d = {
        "efficiency": {name: {} for name, _, _ in GEN_2D_AXES},
        "fake_rate": {name: {} for name, _, _ in RECO_2D_AXES},
        "duplicate_rate": {name: {} for name, _, _ in GEN_2D_AXES},
        "ancestor_not_gamma_status_1_expected_efficiency": {name: {} for name, _, _ in GEN_2D_AXES},
    }
    expected_matched_count = {axis: {} for axis, _, _ in MATCH_AXES}
    expected_missed_count = {axis: {} for axis, _, _ in MATCH_AXES}
    expected_matched_count_2d = {name: {} for name, _, _ in GEN_2D_AXES}
    expected_missed_count_2d = {name: {} for name, _, _ in GEN_2D_AXES}

    source, file_name = MATCH_SOURCES[0]
    gen_px, gen_py, gen_pz = GEN_P4_BRANCHES[source]

    for sample in samples:
        tree = uproot.open(input_root / sample / file_name)["Events"]
        branches = [
            "GenPart_status",
            "GenPart_pdgId",
            "GenPart_parentIdx",
            gen_px,
            gen_py,
            gen_pz,
            "Event_bFieldGevCm",
            "TracRaw_theta",
            "TracRaw_phi",
            "TracRaw_invR",
            "TracRaw_charge",
        ]
        arrays = {branch: ak.to_list(tree[branch].array(library="ak")) for branch in branches}
        expected_pdg = infer_expected_abs_pdg(
            arrays["GenPart_status"],
            arrays["GenPart_pdgId"],
            arrays["GenPart_parentIdx"],
        )

        n_gen = 0
        n_reco = 0
        n_expected = 0
        n_expected_pair_event = 0
        n_matched_by_cut = np.zeros(len(MATCH_ALPHA_CUTS), dtype=int)
        n_duplicate_by_cut = np.zeros(len(MATCH_ALPHA_CUTS), dtype=int)
        n_expected_matched_by_cut = np.zeros(len(MATCH_ALPHA_CUTS), dtype=int)
        n_expected_pair_matched_by_cut = np.zeros(len(MATCH_ALPHA_CUTS), dtype=int)

        gen_den = {axis: np.zeros(len(bins) - 1, dtype=int) for axis, bins, _ in MATCH_AXES}
        gen_matched = {axis: np.zeros_like(values) for axis, values in gen_den.items()}
        gen_duplicate = {axis: np.zeros_like(values) for axis, values in gen_den.items()}
        expected_den = {axis: np.zeros_like(values) for axis, values in gen_den.items()}
        expected_matched = {axis: np.zeros_like(values) for axis, values in gen_den.items()}
        reco_den = {axis: np.zeros(len(bins) - 1, dtype=int) for axis, bins, _ in MATCH_AXES}
        reco_fake = {axis: np.zeros_like(values) for axis, values in reco_den.items()}
        gen_den_2d = {}
        gen_matched_2d = {}
        gen_duplicate_2d = {}
        expected_den_2d = {}
        expected_matched_2d = {}
        reco_den_2d = {}
        reco_fake_2d = {}
        for plane, x_axis, y_axis in GEN_2D_AXES:
            shape = (len(AXIS_BINS[x_axis]) - 1, len(AXIS_BINS[y_axis]) - 1)
            gen_den_2d[plane] = np.zeros(shape, dtype=int)
            gen_matched_2d[plane] = np.zeros(shape, dtype=int)
            gen_duplicate_2d[plane] = np.zeros(shape, dtype=int)
            expected_den_2d[plane] = np.zeros(shape, dtype=int)
            expected_matched_2d[plane] = np.zeros(shape, dtype=int)
            reco_den_2d[plane] = np.zeros(shape, dtype=int)
            reco_fake_2d[plane] = np.zeros(shape, dtype=int)

        for event_idx in range(tree.num_entries):
            gen_particles = []
            pdgs = [int(value) for value in arrays["GenPart_pdgId"][event_idx]]
            parents = [int(value) for value in arrays["GenPart_parentIdx"][event_idx]]

            for idx, (status, pdg_id, px, py, pz) in enumerate(
                zip(
                    arrays["GenPart_status"][event_idx],
                    arrays["GenPart_pdgId"][event_idx],
                    arrays[gen_px][event_idx],
                    arrays[gen_py][event_idx],
                    arrays[gen_pz][event_idx],
                    strict=True,
                )
            ):
                charge = charge_from_pdg(int(pdg_id))
                if status != 1 or charge == 0 or abs(pdg_id) not in CHARGED_ABS_PDG:
                    continue
                pt = float(np.hypot(px, py))
                p = float(np.sqrt(px * px + py * py + pz * pz))
                theta = float(np.arctan2(pt, pz))
                if not is_track_fiducial(theta, pt, fiducial_cos_max, pt_min):
                    continue
                gen_particles.append({
                    "pt": pt,
                    "p": p,
                    "theta": theta,
                    "phi": phi_0_2pi(float(np.arctan2(py, px))),
                    "charge": charge,
                    "is_expected": expected_pdg is not None and abs(int(pdg_id)) == expected_pdg and not has_ancestor_pdg(pdgs, parents, idx, 22),
                })

            reco_tracks = []
            for theta, phi, inv_r, charge in zip(
                arrays["TracRaw_theta"][event_idx],
                arrays["TracRaw_phi"][event_idx],
                arrays["TracRaw_invR"][event_idx],
                arrays["TracRaw_charge"][event_idx],
                strict=True,
            ):
                if charge == 0 or inv_r == 0:
                    continue
                theta = float(theta)
                pt = float(arrays["Event_bFieldGevCm"][event_idx]) / abs(float(inv_r))
                if is_track_fiducial(theta, pt, fiducial_cos_max, pt_min):
                    reco_tracks.append({
                        "pt": pt,
                        "p": pt / np.sin(theta),
                        "theta": theta,
                        "phi": phi_0_2pi(float(phi)),
                        "charge": int(charge),
                    })

            n_gen += len(gen_particles)
            n_reco += len(reco_tracks)
            n_expected += sum(gen["is_expected"] for gen in gen_particles)
            expected_indices = [idx for idx, gen in enumerate(gen_particles) if gen["is_expected"]]
            if len(expected_indices) >= 2:
                n_expected_pair_event += 1

            candidates = []
            distances_by_gen = [[] for _ in gen_particles]
            for gen_idx, gen in enumerate(gen_particles):
                for reco_idx, reco in enumerate(reco_tracks):
                    if gen["charge"] != reco["charge"]:
                        continue
                    alpha = opening_angle(gen["theta"], gen["phi"], reco["theta"], reco["phi"])
                    if alpha < MATCH_ALPHA_CUTS[-1]:
                        candidates.append((alpha, gen_idx, reco_idx))
                        distances_by_gen[gen_idx].append(alpha)

            used_gen_scan = set()
            used_reco_scan = set()
            expected_match_alpha = []
            for alpha, gen_idx, reco_idx in sorted(candidates):
                if gen_idx in used_gen_scan or reco_idx in used_reco_scan:
                    continue
                used_gen_scan.add(gen_idx)
                used_reco_scan.add(reco_idx)
                cut_index = np.searchsorted(MATCH_ALPHA_CUTS, alpha, side="right")
                n_matched_by_cut[cut_index:] += 1
                if gen_particles[gen_idx]["is_expected"]:
                    n_expected_matched_by_cut[cut_index:] += 1
                    expected_match_alpha.append(alpha)

            if len(expected_match_alpha) >= 2:
                pair_alpha = sorted(expected_match_alpha)[1]
                n_expected_pair_matched_by_cut[np.searchsorted(MATCH_ALPHA_CUTS, pair_alpha, side="right") :] += 1

            for distances in distances_by_gen:
                if len(distances) >= 2:
                    second_match = sorted(distances)[1]
                    n_duplicate_by_cut[np.searchsorted(MATCH_ALPHA_CUTS, second_match, side="right") :] += 1

            used_gen = set()
            used_reco = set()
            for alpha, gen_idx, reco_idx in sorted(candidates):
                if alpha >= alpha_cut or gen_idx in used_gen or reco_idx in used_reco:
                    continue
                used_gen.add(gen_idx)
                used_reco.add(reco_idx)

            for gen_idx, gen in enumerate(gen_particles):
                for plane, x_axis, y_axis in GEN_2D_AXES:
                    fill_2d(gen_den_2d[plane], gen[x_axis], gen[y_axis], AXIS_BINS[x_axis], AXIS_BINS[y_axis])
                for axis, bins, _ in MATCH_AXES:
                    fill_1d(gen_den[axis], gen[axis], bins)
                    if gen_idx in used_gen:
                        fill_1d(gen_matched[axis], gen[axis], bins)
                    if sum(alpha < alpha_cut for alpha in distances_by_gen[gen_idx]) >= 2:
                        fill_1d(gen_duplicate[axis], gen[axis], bins)
                    if gen["is_expected"]:
                        fill_1d(expected_den[axis], gen[axis], bins)
                        if gen_idx in used_gen:
                            fill_1d(expected_matched[axis], gen[axis], bins)
                if gen_idx in used_gen:
                    for plane, x_axis, y_axis in GEN_2D_AXES:
                        fill_2d(gen_matched_2d[plane], gen[x_axis], gen[y_axis], AXIS_BINS[x_axis], AXIS_BINS[y_axis])
                if sum(alpha < alpha_cut for alpha in distances_by_gen[gen_idx]) >= 2:
                    for plane, x_axis, y_axis in GEN_2D_AXES:
                        fill_2d(gen_duplicate_2d[plane], gen[x_axis], gen[y_axis], AXIS_BINS[x_axis], AXIS_BINS[y_axis])
                if gen["is_expected"]:
                    for plane, x_axis, y_axis in GEN_2D_AXES:
                        fill_2d(expected_den_2d[plane], gen[x_axis], gen[y_axis], AXIS_BINS[x_axis], AXIS_BINS[y_axis])
                    if gen_idx in used_gen:
                        for plane, x_axis, y_axis in GEN_2D_AXES:
                            fill_2d(expected_matched_2d[plane], gen[x_axis], gen[y_axis], AXIS_BINS[x_axis], AXIS_BINS[y_axis])

            for reco_idx, reco in enumerate(reco_tracks):
                for plane, x_axis, y_axis in RECO_2D_AXES:
                    fill_2d(reco_den_2d[plane], reco[x_axis], reco[y_axis], AXIS_BINS[x_axis], AXIS_BINS[y_axis])
                for axis, bins, _ in MATCH_AXES:
                    fill_1d(reco_den[axis], reco[axis], bins)
                    if reco_idx not in used_reco:
                        fill_1d(reco_fake[axis], reco[axis], bins)
                if reco_idx not in used_reco:
                    for plane, x_axis, y_axis in RECO_2D_AXES:
                        fill_2d(reco_fake_2d[plane], reco[x_axis], reco[y_axis], AXIS_BINS[x_axis], AXIS_BINS[y_axis])

        scan_values["efficiency"][sample] = n_matched_by_cut / n_gen if n_gen else np.zeros(len(MATCH_ALPHA_CUTS))
        scan_values["fake_rate"][sample] = (n_reco - n_matched_by_cut) / n_reco if n_reco else np.zeros(len(MATCH_ALPHA_CUTS))
        scan_values["duplicate_rate"][sample] = n_duplicate_by_cut / n_gen if n_gen else np.zeros(len(MATCH_ALPHA_CUTS))
        scan_values["ancestor_not_gamma_status_1_expected_efficiency"][sample] = n_expected_matched_by_cut / n_expected if n_expected else np.zeros(len(MATCH_ALPHA_CUTS))
        scan_values["ancestor_not_gamma_status_1_expected_pair_efficiency"][sample] = n_expected_pair_matched_by_cut / n_expected_pair_event if n_expected_pair_event else np.zeros(len(MATCH_ALPHA_CUTS))

        metric_pairs = {
            "efficiency": (gen_matched, gen_den),
            "fake_rate": (reco_fake, reco_den),
            "duplicate_rate": (gen_duplicate, gen_den),
            "ancestor_not_gamma_status_1_expected_efficiency": (expected_matched, expected_den),
        }
        for metric, (num_by_axis, den_by_axis) in metric_pairs.items():
            for axis, _, _ in MATCH_AXES:
                metric_num[metric][axis][sample] = num_by_axis[axis]
                metric_den[metric][axis][sample] = den_by_axis[axis]

        for plane, _, _ in GEN_2D_AXES:
            metric_num_2d["efficiency"][plane][sample] = gen_matched_2d[plane]
            metric_den_2d["efficiency"][plane][sample] = gen_den_2d[plane]
            metric_num_2d["duplicate_rate"][plane][sample] = gen_duplicate_2d[plane]
            metric_den_2d["duplicate_rate"][plane][sample] = gen_den_2d[plane]
            metric_num_2d["ancestor_not_gamma_status_1_expected_efficiency"][plane][sample] = expected_matched_2d[plane]
            metric_den_2d["ancestor_not_gamma_status_1_expected_efficiency"][plane][sample] = expected_den_2d[plane]
            expected_matched_count_2d[plane][sample] = expected_matched_2d[plane]
            expected_missed_count_2d[plane][sample] = expected_den_2d[plane] - expected_matched_2d[plane]
        for plane, _, _ in RECO_2D_AXES:
            metric_num_2d["fake_rate"][plane][sample] = reco_fake_2d[plane]
            metric_den_2d["fake_rate"][plane][sample] = reco_den_2d[plane]
        for axis, _, _ in MATCH_AXES:
            expected_matched_count[axis][sample] = expected_matched[axis]
            expected_missed_count[axis][sample] = expected_den[axis] - expected_matched[axis]

    labels = {
        "efficiency": "Efficiency",
        "fake_rate": "Fake rate",
        "duplicate_rate": "Duplicate rate",
        "ancestor_not_gamma_status_1_expected_efficiency": "Efficiency",
        "ancestor_not_gamma_status_1_expected_pair_efficiency": "Efficiency",
    }
    axis_labels = {
        "efficiency": GEN_AXIS_LABELS,
        "fake_rate": RECO_AXIS_LABELS,
        "duplicate_rate": GEN_AXIS_LABELS,
        "ancestor_not_gamma_status_1_expected_efficiency": GEN_AXIS_LABELS,
    }
    for metric, ylabel in labels.items():
        plot_scan(scan_output, metric, ylabel, scan_values[metric], samples, styles)
        if metric == "ancestor_not_gamma_status_1_expected_pair_efficiency":
            continue
        for axis, bins, _ in MATCH_AXES:
            plot_metric(
                cut_output,
                f"{metric}_vs_{axis}",
                axis_labels[metric][axis],
                ylabel,
                bins,
                metric_num[metric][axis],
                metric_den[metric][axis],
                samples,
                styles,
            )

        for plane, x_axis, y_axis in (RECO_2D_AXES if metric == "fake_rate" else GEN_2D_AXES):
            plot_metric_2d(
                cut_output,
                f"{metric}_vs_{plane}",
                PLAIN_AXIS_LABELS[x_axis],
                PLAIN_AXIS_LABELS[y_axis],
                AXIS_BINS[x_axis],
                AXIS_BINS[y_axis],
                ylabel,
                metric_num_2d[metric][plane],
                metric_den_2d[metric][plane],
                samples,
            )

    for state, counts_by_axis, counts_2d in (
        ("matched", expected_matched_count, expected_matched_count_2d),
        ("missed", expected_missed_count, expected_missed_count_2d),
    ):
        for axis, bins, _ in MATCH_AXES:
            plot_count(
                cut_output,
                f"ancestor_not_gamma_status_1_expected_{state}_count_vs_{axis}",
                GEN_AXIS_LABELS[axis],
                bins,
                counts_by_axis[axis],
                expected_matched_count[axis],
                expected_missed_count[axis],
                samples,
                styles,
            )
        for plane, x_axis, y_axis in GEN_2D_AXES:
            plot_count_2d(
                cut_output,
                f"ancestor_not_gamma_status_1_expected_{state}_count_vs_{plane}",
                PLAIN_AXIS_LABELS[x_axis],
                PLAIN_AXIS_LABELS[y_axis],
                AXIS_BINS[x_axis],
                AXIS_BINS[y_axis],
                counts_2d[plane],
                expected_matched_count_2d[plane],
                expected_missed_count_2d[plane],
                samples,
            )

    print(f"wrote {scan_output}")
    print(f"wrote {cut_output}")
