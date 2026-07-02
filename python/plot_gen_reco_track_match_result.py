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
    MATCH_SOURCES,
    add_delphi_label,
    charge_from_pdg,
    has_ancestor_pdg,
    infer_expected_abs_pdg,
    is_track_fiducial,
    match_cut_dir,
    opening_angle,
    phi_0_2pi,
    resolve_samples,
    sample_styles,
    set_hist_yaxis,
    wrap_phi,
)


RESIDUAL_PLOTS = (
    ("pt_response", r"$(p_{T}^{\mathrm{reco}} - p_{T}^{\mathrm{gen}}) / p_{T}^{\mathrm{gen}}$", (-1.0, 1.0)),
    ("p_response", r"$(p^{\mathrm{reco}} - p^{\mathrm{gen}}) / p^{\mathrm{gen}}$", (-1.0, 1.0)),
    ("pt_residual", r"$p_{T}^{\mathrm{reco}} - p_{T}^{\mathrm{gen}}$ [GeV]", (-30.0, 30.0)),
    ("p_residual", r"$p^{\mathrm{reco}} - p^{\mathrm{gen}}$ [GeV]", (-30.0, 30.0)),
    ("theta_residual", r"$\theta_{\mathrm{reco}} - \theta_{\mathrm{gen}}$ [rad]", (-0.05, 0.05)),
    ("phi_residual", r"$\phi_{\mathrm{reco}} - \phi_{\mathrm{gen}}$ [rad]", (-0.05, 0.05)),
    ("alpha_residual", r"$\alpha$ [rad]", None),
)
EXPECTED_RESIDUAL_PLOTS = (
    ("pt_response", r"$(p_{T}^{\mathrm{reco}} - p_{T}^{\mathrm{gen}}) / p_{T}^{\mathrm{gen}}$", (-1.0, 1.0)),
    ("p_response", r"$(p^{\mathrm{reco}} - p^{\mathrm{gen}}) / p^{\mathrm{gen}}$", (-1.0, 1.0)),
    ("pt_residual", r"$p_{T}^{\mathrm{reco}} - p_{T}^{\mathrm{gen}}$ [GeV]", (-30.0, 30.0)),
    ("p_residual", r"$p^{\mathrm{reco}} - p^{\mathrm{gen}}$ [GeV]", (-30.0, 30.0)),
    ("theta_residual", r"$\theta_{\mathrm{reco}} - \theta_{\mathrm{gen}}$ [rad]", (-0.05, 0.05)),
    ("phi_residual", r"$\phi_{\mathrm{reco}} - \phi_{\mathrm{gen}}$ [rad]", (-0.05, 0.05)),
)


def plot_residual(
    output_dir: Path,
    name: str,
    xlabel: str,
    value_range: tuple[float, float],
    values_by_sample: dict[str, list[float]],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    bins = np.linspace(value_range[0], value_range[1], 101)
    counts_list = []

    for sample in samples:
        values = np.asarray(values_by_sample[sample], dtype=float)
        color, hatch = styles[sample]
        underflow = int(np.count_nonzero(values < bins[0]))
        overflow = int(np.count_nonzero(values > bins[-1]))
        if len(values):
            q16, q84 = np.quantile(values, [0.16, 0.84])
            sigma68 = 0.5 * (q84 - q16)
        else:
            sigma68 = float("nan")
        label = rf"{sample:<8} N={len(values):<7} UF={underflow:<5} OF={overflow:<5} $\sigma_{{68}}$={sigma68:.3g}"
        counts, _ = np.histogram(values, bins=bins)
        counts_list.append(counts)
        ax.hist(values, bins=bins, histtype="stepfilled", label=label, facecolor="none", edgecolor=color, hatch=hatch, linewidth=0.0)
        ax.hist(values, bins=bins, histtype="step", label="_nolegend_", color=color, alpha=0.9, linewidth=2)

    set_hist_yaxis(ax, counts_list)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entries")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 14})
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_event_count(
    output_dir: Path,
    name: str,
    xlabel: str,
    values_by_sample: dict[str, list[int]],
    samples: tuple[str, ...],
    styles: dict[str, tuple[str, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    hi = max([value for values in values_by_sample.values() for value in values], default=0)
    bins = np.arange(-0.5, hi + 1.5, 1.0)
    counts_list = []

    for sample in samples:
        values = np.asarray(values_by_sample[sample], dtype=int)
        color, hatch = styles[sample]
        counts, _ = np.histogram(values, bins=bins)
        counts_list.append(counts)
        ax.hist(values, bins=bins, histtype="stepfilled", label=f"{sample:<8} N_event={len(values)}", facecolor="none", edgecolor=color, hatch=hatch, linewidth=0.0)
        ax.hist(values, bins=bins, histtype="step", label="_nolegend_", color=color, alpha=0.9, linewidth=2)

    set_hist_yaxis(ax, counts_list)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Events")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 15})
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_gen_reco_track_match_result(
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
    output_dir = output_root / "gen_reco_track_match_result" / match_cut_dir(alpha_cut, fiducial_cos_max, pt_min)
    output_dir.mkdir(parents=True, exist_ok=True)

    residual_values = {name: {sample: [] for sample in samples} for name, _, _ in RESIDUAL_PLOTS}
    expected_residual_values = {name: {sample: [] for sample in samples} for name, _, _ in EXPECTED_RESIDUAL_PLOTS}
    n_matched_event = {sample: [] for sample in samples}
    n_matched_expected_event = {sample: [] for sample in samples}

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

            candidates = []
            for gen_idx, gen in enumerate(gen_particles):
                for reco_idx, reco in enumerate(reco_tracks):
                    if gen["charge"] != reco["charge"]:
                        continue
                    alpha = opening_angle(gen["theta"], gen["phi"], reco["theta"], reco["phi"])
                    if alpha < alpha_cut:
                        candidates.append((alpha, gen_idx, reco_idx))

            used_gen = set()
            used_reco = set()
            for alpha, gen_idx, reco_idx in sorted(candidates):
                if gen_idx in used_gen or reco_idx in used_reco:
                    continue
                used_gen.add(gen_idx)
                used_reco.add(reco_idx)

                gen = gen_particles[gen_idx]
                reco = reco_tracks[reco_idx]
                residual_values["pt_response"][sample].append((reco["pt"] - gen["pt"]) / gen["pt"])
                residual_values["p_response"][sample].append((reco["p"] - gen["p"]) / gen["p"])
                residual_values["pt_residual"][sample].append(reco["pt"] - gen["pt"])
                residual_values["p_residual"][sample].append(reco["p"] - gen["p"])
                residual_values["theta_residual"][sample].append(reco["theta"] - gen["theta"])
                residual_values["phi_residual"][sample].append(wrap_phi(reco["phi"] - gen["phi"]))
                residual_values["alpha_residual"][sample].append(alpha)
                if gen["is_expected"]:
                    expected_residual_values["pt_response"][sample].append((reco["pt"] - gen["pt"]) / gen["pt"])
                    expected_residual_values["p_response"][sample].append((reco["p"] - gen["p"]) / gen["p"])
                    expected_residual_values["pt_residual"][sample].append(reco["pt"] - gen["pt"])
                    expected_residual_values["p_residual"][sample].append(reco["p"] - gen["p"])
                    expected_residual_values["theta_residual"][sample].append(reco["theta"] - gen["theta"])
                    expected_residual_values["phi_residual"][sample].append(wrap_phi(reco["phi"] - gen["phi"]))

            n_matched_event[sample].append(len(used_reco))
            n_matched_expected_event[sample].append(sum(gen_particles[gen_idx]["is_expected"] for gen_idx in used_gen))

    for plot_name, xlabel, value_range in RESIDUAL_PLOTS:
        if value_range is None:
            value_range = (0.0, alpha_cut)
        plot_residual(output_dir, plot_name, xlabel, value_range, residual_values[plot_name], samples, styles)
    for plot_name, xlabel, value_range in EXPECTED_RESIDUAL_PLOTS:
        plot_residual(
            output_dir,
            f"ancestor_not_gamma_status_1_expected_track_{plot_name}",
            xlabel,
            value_range,
            expected_residual_values[plot_name],
            samples,
            styles,
        )

    plot_event_count(output_dir, "n_matched_track", "Number of matched reco tracks per event", n_matched_event, samples, styles)
    plot_event_count(
        output_dir,
        "n_matched_ancestor_not_gamma_status_1_expected_track",
        "Number of matched ancestor_not_gamma_status_1_expected tracks per event",
        n_matched_expected_event,
        samples,
        styles,
    )

    print(f"wrote {output_dir}")
