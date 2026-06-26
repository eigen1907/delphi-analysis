from __future__ import annotations

import shutil
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from plot_utils import (
    CHARGED_ABS_PDG,
    EXPECTED_ABS_PDG,
    MATCH_P_BINS,
    MATCH_PHI_BINS,
    MATCH_PT_BINS,
    MATCH_THETA_BINS,
    SAMPLE_LABELS,
    SAMPLE_STYLES,
    add_delphi_label,
    charge_from_pdg,
    has_ancestor_pdg,
    phi_0_2pi,
    set_hist_yaxis,
    valid_index,
)


SAMPLES = SAMPLE_LABELS
Z_MASS = 91.2

BRANCH_PLOTS = (
    ("px", r"$p_x$ [GeV]", "GenPart_vector.fCoordinates.fX", (-60.0, 60.0)),
    ("py", r"$p_y$ [GeV]", "GenPart_vector.fCoordinates.fY", (-60.0, 60.0)),
    ("pz", r"$p_z$ [GeV]", "GenPart_vector.fCoordinates.fZ", (-60.0, 60.0)),
    ("energy", r"$E$ [GeV]", "GenPart_vector.fCoordinates.fT", (0.0, 80.0)),
    ("mass", r"$m$ [GeV]", "GenPart_mass", (0.0, 5.0)),
    ("vertex_x", r"$x_{\mathrm{vtx}}$", "GenPart_vertex.fCoordinates.fX", (-0.2, 0.2)),
    ("vertex_y", r"$y_{\mathrm{vtx}}$", "GenPart_vertex.fCoordinates.fY", (-0.002, 0.002)),
    ("vertex_z", r"$z_{\mathrm{vtx}}$", "GenPart_vertex.fCoordinates.fZ", (-2.5, 2.5)),
    ("vertex_t", r"$t_{\mathrm{vtx}}$", "GenPart_vertex.fCoordinates.fT", (-1e-6, 1e-6)),
)
DERIVED_PLOTS = (
    ("pt", r"$p_T$ [GeV]", (0.0, 60.0)),
    ("p", r"$p$ [GeV]", (0.0, 60.0)),
    ("theta", r"$\theta$ [rad]", (0.0, np.pi)),
    ("phi", r"$\phi$ [rad]", (0.0, 2.0 * np.pi)),
    ("leading_pt", r"Leading $p_T$ [GeV]", (0.0, 60.0)),
    ("leading_p", r"Leading $p$ [GeV]", (0.0, 60.0)),
    ("subleading_pt", r"Subleading $p_T$ [GeV]", (0.0, 60.0)),
    ("subleading_p", r"Subleading $p$ [GeV]", (0.0, 60.0)),
    ("pair_mass", r"Pair mass [GeV]", (0.0, 120.0)),
    ("pair_pt", r"Pair $p_T$ [GeV]", (0.0, 20.0)),
    ("pair_p", r"Pair $p$ [GeV]", (0.0, 60.0)),
    ("pair_theta", r"Pair $\theta$ [rad]", (0.0, np.pi)),
    ("pair_phi", r"Pair $\phi$ [rad]", (0.0, 2.0 * np.pi)),
    ("pair_px", r"Pair $p_x$ [GeV]", (-20.0, 20.0)),
    ("pair_py", r"Pair $p_y$ [GeV]", (-20.0, 20.0)),
    ("pair_pz", r"Pair $p_z$ [GeV]", (-120.0, 120.0)),
)
PLOTS_2D = (
    ("pt_theta", "pt", "theta", r"$p_T$ [GeV]", r"$\theta$ [rad]", MATCH_PT_BINS, MATCH_THETA_BINS),
    ("p_theta", "p", "theta", r"$p$ [GeV]", r"$\theta$ [rad]", MATCH_P_BINS, MATCH_THETA_BINS),
    ("theta_phi", "theta", "phi", r"$\theta$ [rad]", r"$\phi$ [rad]", MATCH_THETA_BINS, MATCH_PHI_BINS),
)
LEGEND_HEADER = f"{'Sample':<8} {'Event':<8} {'Object':<8} {'UF':<6} {'OF':<6}"
SELECTIONS = (
    ("all", "all"),
    ("no_parent", "no parent"),
    ("status_1_charged", "status=1, charged"),
    ("status_1_expected", "status=1, expected species"),
    ("ancestor_gamma_status_1_expected", "photon ancestor, status=1, expected species"),
    ("ancestor_not_gamma_status_1_expected", "no photon ancestor, status=1, expected species"),
    ("ancestor_not_gamma_status_1_expected_origin", "origin of no photon ancestor, status=1 expected species"),
)


def flatten(array) -> np.ndarray:
    values = np.asarray(ak.to_numpy(ak.flatten(array, axis=None)), dtype=float)
    return values[np.isfinite(values)]


def plot_hist(output_root: Path, name: str, xlabel: str, value_range: tuple[float, float], values_by_sample: dict[str, np.ndarray], n_events: dict[str, int]) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    bins = np.linspace(value_range[0], value_range[1], 101)

    ax.plot([], [], color="none", label=LEGEND_HEADER)
    counts_list = []
    density_list = []
    for sample in SAMPLES:
        values = values_by_sample[sample]
        color, hatch = SAMPLE_STYLES[sample]
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
    fig.savefig(output_root / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_discrete(output_root: Path, name: str, xlabel: str, values_by_sample: dict[str, np.ndarray], n_events: dict[str, int]) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    unique_values = sorted({int(value) for values in values_by_sample.values() for value in values})
    if not unique_values:
        unique_values = [0]
    x = np.arange(len(unique_values))

    ax.plot([], [], color="none", label=LEGEND_HEADER)
    counts_list = []
    fraction_list = []
    for sample in SAMPLES:
        values = values_by_sample[sample].astype(int)
        raw_counts = np.asarray([np.count_nonzero(values == value) for value in unique_values], dtype=float)
        n_object = int(np.sum(raw_counts))
        counts = raw_counts / n_object if n_object else raw_counts
        counts_list.append(raw_counts)
        fraction_list.append(counts)
        color, hatch = SAMPLE_STYLES[sample]
        label = f"{sample:<8} {n_events[sample]:<8} {n_object:<8} {0:<6} {0:<6}"
        ax.bar(x, counts, width=0.8, label=label, facecolor="none", edgecolor=color, hatch=hatch, alpha=0.9)

    set_hist_yaxis(ax, counts_list, fraction_list)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [str(value) for value in unique_values],
        rotation=90 if len(unique_values) > 12 else 0,
        fontsize=13 if len(unique_values) > 12 else None,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalized")
    ax.legend(loc="upper right", prop={"family": "monospace", "size": 15})
    ax.grid(axis="y", alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_root / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_int_hist(output_root: Path, name: str, xlabel: str, values_by_sample: dict[str, np.ndarray], n_events: dict[str, int]) -> None:
    fig, ax = plt.subplots(figsize=(12, 10))
    lo = int(min(np.min(values) for values in values_by_sample.values()))
    hi = int(max(np.max(values) for values in values_by_sample.values()))
    bins = np.arange(lo - 0.5, hi + 1.5, 1.0)

    ax.plot([], [], color="none", label=LEGEND_HEADER)
    counts_list = []
    density_list = []
    for sample in SAMPLES:
        values = values_by_sample[sample]
        color, hatch = SAMPLE_STYLES[sample]
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
    fig.savefig(output_root / f"{name}.png", dpi=150)
    plt.close(fig)


def plot_2d(
    output_root: Path,
    name: str,
    x_label: str,
    y_label: str,
    x_bins: np.ndarray,
    y_bins: np.ndarray,
    x_values_by_sample: dict[str, np.ndarray],
    y_values_by_sample: dict[str, np.ndarray],
) -> None:
    counts_by_sample = {}
    for sample in SAMPLES:
        counts, _, _ = np.histogram2d(x_values_by_sample[sample], y_values_by_sample[sample], bins=(x_bins, y_bins))
        counts_by_sample[sample] = counts

    vmax = max(1.0, max(float(np.max(counts)) for counts in counts_by_sample.values()))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="white", alpha=0.0)
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), sharex=True, sharey=True)
    mesh = None

    for ax, sample in zip(axes.ravel(), SAMPLES, strict=True):
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

    for ax in axes[-1, :]:
        ax.set_xlabel(x_label)
    for ax in axes[:, 0]:
        ax.set_ylabel(y_label)

    fig.subplots_adjust(right=0.86, hspace=0.08, wspace=0.08)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.025, 0.70])
    cbar = fig.colorbar(mesh, cax=cbar_ax)
    cbar.set_label("Entries")
    fig.savefig(output_root / f"{name}.png", dpi=150)
    plt.close(fig)


def expected_origin_indices(statuses, pdgs, parents, sample: str) -> set[int]:
    expected_abs_pdg = EXPECTED_ABS_PDG[sample]
    origins = set()

    for idx, status in enumerate(statuses):
        pdg = int(pdgs[idx])
        if int(status) != 1 or abs(pdg) != expected_abs_pdg:
            continue
        if has_ancestor_pdg(pdgs, parents, idx, 22):
            continue

        current = idx
        seen = set()
        while True:
            parent_idx = int(parents[current])
            if not valid_index(parent_idx, len(pdgs)) or parent_idx in seen:
                break
            if abs(int(pdgs[parent_idx])) != expected_abs_pdg:
                break
            seen.add(current)
            current = parent_idx

        origins.add(current)

    return origins


def make_selections(sample_arrays: dict) -> tuple[tuple[str, str], ...]:
    statuses = sorted(
        {
            int(status)
            for sample in SAMPLES
            for event in ak.to_list(sample_arrays[sample]["GenPart_status"])
            for status in event
        }
    )
    return SELECTIONS + tuple((f"status_{status}", f"status={status}") for status in statuses)


def make_selection_mask(arrays: dict, sample: str, selection: str) -> ak.Array:
    status_events = ak.to_list(arrays["GenPart_status"])
    pdg_events = ak.to_list(arrays["GenPart_pdgId"])
    parent_events = ak.to_list(arrays["GenPart_parentIdx"])
    masks = []

    for statuses, pdgs, parents in zip(
        status_events,
        pdg_events,
        parent_events,
        strict=True,
    ):
        event_mask = []
        origin_indices = (
            expected_origin_indices(statuses, pdgs, parents, sample)
            if selection == "ancestor_not_gamma_status_1_expected_origin"
            else set()
        )

        for idx, status in enumerate(statuses):
            status = int(status)
            pdg = int(pdgs[idx])
            parent_idx = int(parents[idx])
            has_parent = valid_index(parent_idx, len(pdgs))
            has_photon_ancestor = has_ancestor_pdg(pdgs, parents, idx, 22)

            if selection == "all":
                keep = True
            elif selection == "no_parent":
                keep = not has_parent
            elif selection == "status_1_charged":
                keep = status == 1 and charge_from_pdg(pdg) != 0
            elif selection == "status_1_expected":
                keep = status == 1 and abs(pdg) == EXPECTED_ABS_PDG[sample]
            elif selection == "ancestor_gamma_status_1_expected":
                keep = has_photon_ancestor and status == 1 and abs(pdg) == EXPECTED_ABS_PDG[sample]
            elif selection == "ancestor_not_gamma_status_1_expected":
                keep = not has_photon_ancestor and status == 1 and abs(pdg) == EXPECTED_ABS_PDG[sample]
            elif selection == "ancestor_not_gamma_status_1_expected_origin":
                keep = idx in origin_indices
            elif selection.startswith("status_"):
                keep = status == int(selection.removeprefix("status_"))
            else:
                keep = False

            event_mask.append(keep)

        masks.append(event_mask)

    return ak.Array(masks)


def plot_gen_check(input_root: Path, output_root: Path) -> None:
    mh.style.use(mh.styles.CMS)
    output_root = output_root / "gen_check"
    output_root.mkdir(parents=True, exist_ok=True)
    for path in output_root.glob("*.png"):
        path.unlink()

    n_events = {}
    sample_arrays = {}

    for sample in SAMPLES:
        tree = uproot.open(input_root / sample / "nanoaod.root")["Events"]
        n_events[sample] = int(tree.num_entries)
        arrays = {branch: tree[branch].array(library="ak") for _, _, branch, _ in BRANCH_PLOTS}
        arrays["nGenPart"] = tree["nGenPart"].array(library="ak")
        arrays["GenPart_status"] = tree["GenPart_status"].array(library="ak")
        arrays["GenPart_pdgId"] = tree["GenPart_pdgId"].array(library="ak")
        arrays["GenPart_parentIdx"] = tree["GenPart_parentIdx"].array(library="ak")
        sample_arrays[sample] = arrays

    selections = make_selections(sample_arrays)
    selection_dirs = {selection_dir for selection_dir, _ in selections}
    for path in output_root.iterdir():
        if path.is_dir() and path.name not in selection_dirs:
            shutil.rmtree(path)

    for selection_dir, selection_label in selections:
        selection_output = output_root / selection_dir
        shutil.rmtree(selection_output, ignore_errors=True)
        selection_output.mkdir(parents=True, exist_ok=True)

        values = {
            "n_gen_part": {},
            "pdg_id": {},
            "status": {},
            "pt": {},
            "p": {},
            "theta": {},
            "phi": {},
            "leading_pt": {},
            "leading_p": {},
            "subleading_pt": {},
            "subleading_p": {},
            "pair_mass": {},
            "pair_pt": {},
            "pair_p": {},
            "pair_theta": {},
            "pair_phi": {},
            "pair_px": {},
            "pair_py": {},
            "pair_pz": {},
        }
        for plot_name, _, _, _ in BRANCH_PLOTS:
            values[plot_name] = {}

        for sample in SAMPLES:
            arrays = sample_arrays[sample]
            status = arrays["GenPart_status"]
            pdg_id = arrays["GenPart_pdgId"]
            mask = make_selection_mask(arrays, sample, selection_dir)

            values["n_gen_part"][sample] = np.asarray(ak.to_numpy(ak.sum(mask, axis=1)), dtype=float)
            values["pdg_id"][sample] = flatten(pdg_id[mask])
            values["status"][sample] = flatten(status[mask])
            for plot_name, _, branch, _ in BRANCH_PLOTS:
                values[plot_name][sample] = flatten(arrays[branch][mask])

            px_events = ak.to_list(arrays["GenPart_vector.fCoordinates.fX"][mask])
            py_events = ak.to_list(arrays["GenPart_vector.fCoordinates.fY"][mask])
            pz_events = ak.to_list(arrays["GenPart_vector.fCoordinates.fZ"][mask])
            e_events = ak.to_list(arrays["GenPart_vector.fCoordinates.fT"][mask])
            pdg_events = ak.to_list(pdg_id[mask])

            pt_values = []
            p_values = []
            theta_values = []
            phi_values = []
            leading_pt_values = []
            leading_p_values = []
            subleading_pt_values = []
            subleading_p_values = []
            pair_mass_values = []
            pair_pt_values = []
            pair_p_values = []
            pair_theta_values = []
            pair_phi_values = []
            pair_px_values = []
            pair_py_values = []
            pair_pz_values = []
            for pxs, pys, pzs, es, pdg_ids in zip(px_events, py_events, pz_events, e_events, pdg_events, strict=True):
                particles = []
                for px, py, pz, energy, pdg in zip(pxs, pys, pzs, es, pdg_ids, strict=True):
                    pt = float(np.hypot(px, py))
                    p = float(np.sqrt(px * px + py * py + pz * pz))
                    pt_values.append(pt)
                    p_values.append(p)
                    theta_values.append(float(np.arctan2(pt, pz)))
                    phi_values.append(phi_0_2pi(float(np.arctan2(py, px))))
                    charge = charge_from_pdg(int(pdg))
                    if charge != 0 and abs(pdg) in CHARGED_ABS_PDG:
                        particles.append({"px": px, "py": py, "pz": pz, "energy": energy, "pt": pt, "p": p, "charge": charge})

                particles = sorted(particles, key=lambda item: item["pt"], reverse=True)
                if len(particles) >= 1:
                    leading_pt_values.append(particles[0]["pt"])
                    leading_p_values.append(particles[0]["p"])
                if len(particles) >= 2:
                    subleading_pt_values.append(particles[1]["pt"])
                    subleading_p_values.append(particles[1]["p"])

                    best_mass = None
                    best_distance = None
                    best_pair = None
                    for i, particle_i in enumerate(particles):
                        for particle_j in particles[i + 1 :]:
                            if particle_i["charge"] * particle_j["charge"] >= 0:
                                continue
                            energy = particle_i["energy"] + particle_j["energy"]
                            px = particle_i["px"] + particle_j["px"]
                            py = particle_i["py"] + particle_j["py"]
                            pz = particle_i["pz"] + particle_j["pz"]
                            mass = float(np.sqrt(max(energy * energy - px * px - py * py - pz * pz, 0.0)))
                            distance = abs(mass - Z_MASS)
                            if best_distance is None or distance < best_distance:
                                best_mass = mass
                                best_distance = distance
                                best_pair = (px, py, pz)

                    if best_mass is not None:
                        pair_mass_values.append(best_mass)
                        pair_px, pair_py, pair_pz = best_pair
                        pair_pt = float(np.hypot(pair_px, pair_py))
                        pair_px_values.append(pair_px)
                        pair_py_values.append(pair_py)
                        pair_pz_values.append(pair_pz)
                        pair_pt_values.append(pair_pt)
                        pair_p_values.append(float(np.sqrt(pair_px * pair_px + pair_py * pair_py + pair_pz * pair_pz)))
                        pair_theta_values.append(float(np.arctan2(pair_pt, pair_pz)))
                        pair_phi_values.append(phi_0_2pi(float(np.arctan2(pair_py, pair_px))))

            values["pt"][sample] = np.asarray(pt_values, dtype=float)
            values["p"][sample] = np.asarray(p_values, dtype=float)
            values["theta"][sample] = np.asarray(theta_values, dtype=float)
            values["phi"][sample] = np.asarray(phi_values, dtype=float)
            values["leading_pt"][sample] = np.asarray(leading_pt_values, dtype=float)
            values["leading_p"][sample] = np.asarray(leading_p_values, dtype=float)
            values["subleading_pt"][sample] = np.asarray(subleading_pt_values, dtype=float)
            values["subleading_p"][sample] = np.asarray(subleading_p_values, dtype=float)
            values["pair_mass"][sample] = np.asarray(pair_mass_values, dtype=float)
            values["pair_pt"][sample] = np.asarray(pair_pt_values, dtype=float)
            values["pair_p"][sample] = np.asarray(pair_p_values, dtype=float)
            values["pair_theta"][sample] = np.asarray(pair_theta_values, dtype=float)
            values["pair_phi"][sample] = np.asarray(pair_phi_values, dtype=float)
            values["pair_px"][sample] = np.asarray(pair_px_values, dtype=float)
            values["pair_py"][sample] = np.asarray(pair_py_values, dtype=float)
            values["pair_pz"][sample] = np.asarray(pair_pz_values, dtype=float)

        plot_discrete(selection_output, "pdg_id", "PDG ID", values["pdg_id"], n_events)
        plot_discrete(selection_output, "status", "Status", values["status"], n_events)
        plot_int_hist(selection_output, "n_gen_part", r"$N_{\mathrm{GenPart}}$", values["n_gen_part"], n_events)

        for plot_name, xlabel, _, value_range in BRANCH_PLOTS:
            plot_hist(selection_output, plot_name, xlabel, value_range, values[plot_name], n_events)
        for plot_name, xlabel, value_range in DERIVED_PLOTS:
            plot_hist(selection_output, plot_name, xlabel, value_range, values[plot_name], n_events)
        for plot_name, x_name, y_name, x_label, y_label, x_bins, y_bins in PLOTS_2D:
            plot_2d(selection_output, plot_name, x_label, y_label, x_bins, y_bins, values[x_name], values[y_name])

    print(f"wrote {output_root}")
