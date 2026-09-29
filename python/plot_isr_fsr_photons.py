from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from plot_utils import add_delphi_label, opening_angle, phi_0_2pi, sample_styles, set_hist_yaxis


SAMPLES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")
TARGET_PDG = {"Zee": 11, "Zmumu": 13, "Ztautau": 15, "ZKK": 321, "Zpipi": 211}
MIN_ENERGY = 0.1  # GeV; gen photon selection for the ISR/FSR study
MAX_ANGLE = 0.05  # rad; same angular cut as the track matching study
ENERGY_BINS = np.geomspace(MIN_ENERGY, 50.0, 41)
RECO_ENERGY_BINS = np.geomspace(0.01, 100.0, 41)
COS_BINS = np.linspace(-1.0, 1.0, 41)
PHI_BINS = np.linspace(0.0, 2.0 * np.pi, 41)
POPULATIONS = (
    "stable_gen",
    "reco",
    "isr",
    "fsr",
    "isr_matched_reco",
    "fsr_matched_reco",
)
POPULATION_LABELS = {
    "stable_gen": "Stable gen photons",
    "reco": "Reco photons",
    "isr": "Gen ISR photons",
    "fsr": "Gen FSR photons",
    "isr_matched_reco": "Reco photons matched to gen ISR",
    "fsr_matched_reco": "Reco photons matched to gen FSR",
}

GEN_BRANCHES = (
    "GenPart_status",
    "GenPart_pdgId",
    "GenPart_parentIdx",
    "GenPart_vector.fCoordinates.fX",
    "GenPart_vector.fCoordinates.fY",
    "GenPart_vector.fCoordinates.fZ",
    "GenPart_vector.fCoordinates.fT",
)
RECO_BRANCHES = (
    "Photon_fourMomentum.fCoordinates.fX",
    "Photon_fourMomentum.fCoordinates.fY",
    "Photon_fourMomentum.fCoordinates.fZ",
    "Photon_fourMomentum.fCoordinates.fT",
)
TRUTH_LINK_BRANCHES = ("Photon_partIdx", "Part_simIdx", "SimPart_genIdx")


def photon_origin(index: int, pdgs: list[int], parents: list[int], hard_parent: int, target_pdg: int) -> str | None:
    ancestors = []
    current = parents[index]
    while 0 <= current < len(pdgs) and current not in ancestors:
        ancestors.append(current)
        current = parents[current]

    emitter = next((ancestor for ancestor in ancestors if pdgs[ancestor] != 22), None)
    if emitter is None:
        return None
    if emitter == hard_parent or (parents[emitter] == hard_parent and abs(pdgs[emitter]) == target_pdg):
        return "fsr"
    if hard_parent not in ancestors and abs(pdgs[emitter]) == 11:
        return "isr"
    return None


def direction(px: float, py: float, pz: float) -> tuple[float, float, float] | None:
    momentum = float(np.sqrt(px * px + py * py + pz * pz))
    if not np.isfinite(momentum) or momentum == 0:
        return None
    cos_theta = float(pz / momentum)
    theta = float(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
    phi = phi_0_2pi(float(np.arctan2(py, px)))
    return cos_theta, theta, phi


def record_photon(values: dict, photon: tuple[float, float, float, float]) -> None:
    energy, cos_theta, _, phi = photon
    values["energy"].append(energy)
    values["cos_theta"].append(cos_theta)
    values["phi"].append(phi)


def linked_gen_index(reco_index: int, photon_parts, part_sims, sim_gens) -> int | None:
    part_index = int(photon_parts[reco_index])
    if not 0 <= part_index < len(part_sims):
        return None
    sim_index = int(part_sims[part_index])
    if not 0 <= sim_index < len(sim_gens):
        return None
    gen_index = int(sim_gens[sim_index])
    return gen_index if gen_index >= 0 else None


def read_sample(input_root: Path, sample: str) -> dict:
    sample_dir = input_root / f"20260828_100kTest_{sample}_photosFSR" / "final_root"
    values = {kind: {name: [] for name in ("multiplicity", "energy", "cos_theta", "phi")} for kind in POPULATIONS}
    for kind in ("isr", "fsr"):
        values[kind]["matched_energy"] = []
        values[kind]["matched_cos_theta"] = []
    values["isr"]["total_energy"] = []
    values["isr"]["nearest_angle"] = []
    values["isr"]["energy_response"] = []
    n_events = 0
    n_linked_isr = 0
    target_pdg = TARGET_PDG[sample]
    target_status = 21 if sample == "Ztautau" else 1

    for path in sorted(sample_dir.glob("job_*/nanoaod.root")):
        with uproot.open(path) as root_file:
            tree = root_file["Events"]
            branches = {name: tree[name].array(library="ak") for name in (*GEN_BRANCHES, *RECO_BRANCHES, *TRUTH_LINK_BRANCHES)}

            for event in range(tree.num_entries):
                n_events += 1
                pdgs = [int(pdg) for pdg in branches["GenPart_pdgId"][event]]
                statuses = [int(status) for status in branches["GenPart_status"][event]]
                parents = [int(parent) for parent in branches["GenPart_parentIdx"][event]]
                daughters = [index for index, (pdg, status) in enumerate(zip(pdgs, statuses, strict=True)) if abs(pdg) == target_pdg and status == target_status]
                if len(daughters) != 2 or parents[daughters[0]] != parents[daughters[1]]:
                    raise ValueError(f"Expected a shared hard parent in {path}, event {event}")
                hard_parent = parents[daughters[0]]

                stable_gen_photons = []
                gen_photons = []
                for index, (pdg, status, px, py, pz, energy) in enumerate(zip(
                    pdgs,
                    statuses,
                    branches[GEN_BRANCHES[3]][event],
                    branches[GEN_BRANCHES[4]][event],
                    branches[GEN_BRANCHES[5]][event],
                    branches[GEN_BRANCHES[6]][event],
                    strict=True,
                )):
                    if pdg != 22 or status != 1:
                        continue
                    angles = direction(float(px), float(py), float(pz))
                    if angles is None:
                        continue
                    photon = (float(energy), *angles)
                    stable_gen_photons.append(photon)
                    record_photon(values["stable_gen"], photon)
                    if energy < MIN_ENERGY:
                        continue
                    origin = photon_origin(index, pdgs, parents, hard_parent, target_pdg)
                    if origin is None:
                        continue
                    gen_photons.append((origin, *photon, index))

                reco_photons = []
                reco_gen_links = []
                photon_parts = branches["Photon_partIdx"][event]
                part_sims = branches["Part_simIdx"][event]
                sim_gens = branches["SimPart_genIdx"][event]
                for original_index, (px, py, pz, energy) in enumerate(zip(*(branches[name][event] for name in RECO_BRANCHES), strict=True)):
                    angles = direction(float(px), float(py), float(pz))
                    if angles is not None:
                        photon = (float(energy), *angles)
                        reco_photons.append(photon)
                        reco_gen_links.append(linked_gen_index(original_index, photon_parts, part_sims, sim_gens))
                        record_photon(values["reco"], photon)

                candidates = []
                nearest_angles = [np.inf] * len(gen_photons)
                for gen_index, (_, _, _, gen_theta, gen_phi, original_index) in enumerate(gen_photons):
                    for reco_index, (_, _, reco_theta, reco_phi) in enumerate(reco_photons):
                        angle = opening_angle(gen_theta, gen_phi, reco_theta, reco_phi)
                        # Keep the nearest angle before either matching requirement.
                        nearest_angles[gen_index] = min(nearest_angles[gen_index], angle)
                        linked_index = reco_gen_links[reco_index]
                        # Only a known disagreement vetoes an angular match.
                        if angle < MAX_ANGLE and (linked_index is None or linked_index == original_index):
                            candidates.append((angle, gen_index, reco_index))
                matched_reco = {}
                used_reco = set()
                for _, gen_index, reco_index in sorted(candidates):
                    if gen_index not in matched_reco and reco_index not in used_reco:
                        matched_reco[gen_index] = reco_index
                        used_reco.add(reco_index)
                assert len(matched_reco) == len(used_reco)

                values["stable_gen"]["multiplicity"].append(len(stable_gen_photons))
                values["reco"]["multiplicity"].append(len(reco_photons))
                values["isr"]["total_energy"].append(sum(photon[1] for photon in gen_photons if photon[0] == "isr"))
                for kind in ("isr", "fsr"):
                    values[kind]["multiplicity"].append(sum(photon[0] == kind for photon in gen_photons))
                    values[f"{kind}_matched_reco"]["multiplicity"].append(sum(
                        photon[0] == kind and index in matched_reco
                        for index, photon in enumerate(gen_photons)
                    ))
                for index, photon in enumerate(gen_photons):
                    kind, energy, cos_theta, _, _, _ = photon
                    record_photon(values[kind], photon[1:5])
                    if kind == "isr" and np.isfinite(nearest_angles[index]):
                        values["isr"]["nearest_angle"].append(nearest_angles[index])
                    if index in matched_reco:
                        values[kind]["matched_energy"].append(energy)
                        values[kind]["matched_cos_theta"].append(cos_theta)
                        reco_index = matched_reco[index]
                        reco_photon = reco_photons[reco_index]
                        record_photon(values[f"{kind}_matched_reco"], reco_photon)
                        if kind == "isr":
                            values["isr"]["energy_response"].append(reco_photon[0] / energy)
                            if reco_gen_links[reco_index] is not None:
                                n_linked_isr += 1

    n_isr = len(values["isr"]["energy"])
    n_matched = len(values["isr"]["matched_energy"])
    assert n_matched <= n_isr
    assert n_matched == len(values["isr_matched_reco"]["energy"])
    print(
        f"{sample}: events={n_events}, stable gen photons={len(values['stable_gen']['energy'])}, "
        f"gen ISR={n_isr}, matched ISR={n_matched} ({n_matched / n_isr:.1%}), "
        f"truth-linked={n_linked_isr}, unlinked={n_matched - n_linked_isr}, "
        f"no reco for ISR={n_isr - len(values['isr']['nearest_angle'])}"
    )
    return values


def plot_distribution(output_dir: Path, kind: str, name: str, values_by_sample: dict, styles: dict) -> None:
    level = "Reco" if kind == "reco" or kind.endswith("matched_reco") else "Gen"
    if name == "energy":
        if kind == "stable_gen":
            minimum = min(energy for sample in SAMPLES for energy in values_by_sample[sample][kind]["energy"] if energy > 0)
            maximum = max(energy for sample in SAMPLES for energy in values_by_sample[sample][kind]["energy"])
            bins = np.geomspace(10 ** np.floor(np.log10(minimum)), 10 ** np.ceil(np.log10(maximum)), 61)
        else:
            bins = RECO_ENERGY_BINS if level == "Reco" else ENERGY_BINS
        xlabel = rf"{level} $E_\gamma$ [GeV]"
    elif name == "cos_theta":
        bins, xlabel = COS_BINS, rf"{level} $\cos\theta_\gamma$"
    else:
        bins, xlabel = PHI_BINS, rf"{level} $\phi_\gamma$ [rad]"
    fig, ax = plt.subplots(figsize=(12, 10))
    counts_list = []
    fractions_list = []
    for sample in SAMPLES:
        values = np.asarray(values_by_sample[sample][kind][name])
        counts = np.histogram(values, bins=bins)[0]
        counts_list.append(counts)
        fractions_list.append(counts / len(values))
        color, _ = styles[sample]
        if len(values):
            if name == "energy":
                ax.hist(values, bins=bins, weights=np.full(len(values), 1.0 / len(values)), histtype="step", color=color, linewidth=2, label=f"{sample} (N={len(values)})")
            else:
                ax.hist(values, bins=bins, density=True, histtype="step", color=color, linewidth=2, label=f"{sample} (N={len(values)})")
        else:
            ax.plot([], [], color=color, label=f"{sample} (N=0)")
    if name == "energy":
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fraction of photons / bin" if name == "energy" else "Density")
    ax.legend(title=POPULATION_LABELS[kind])
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    set_hist_yaxis(ax, counts_list, fractions_list if name == "energy" else None)
    fig.tight_layout()
    fig.savefig(output_dir / f"{kind}_{name}.png", dpi=150)
    plt.close(fig)


def plot_multiplicity(output_dir: Path, kind: str, values_by_sample: dict, styles: dict) -> None:
    maximum = max(max(values_by_sample[sample][kind]["multiplicity"], default=0) for sample in SAMPLES)
    bins = np.arange(-0.5, maximum + 1.5)
    fig, ax = plt.subplots(figsize=(12, 10))
    counts_list = []
    fractions_list = []
    for sample in SAMPLES:
        values = np.asarray(values_by_sample[sample][kind]["multiplicity"])
        counts, _ = np.histogram(values, bins=bins)
        counts_list.append(counts)
        fractions_list.append(counts / len(values))
        color, _ = styles[sample]
        ax.hist(values, bins=bins, weights=np.full(len(values), 1.0 / len(values)), histtype="step", color=color, linewidth=2, label=sample)
    if maximum <= 12:
        ax.set_xticks(range(maximum + 1))
    xlabel = f"{POPULATION_LABELS[kind]} per event"
    if kind in ("isr", "fsr"):
        xlabel += fr" ($E_\gamma \geq {MIN_ENERGY}$ GeV)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fraction of events")
    ax.legend()
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    set_hist_yaxis(ax, counts_list, fractions_list)
    fig.tight_layout()
    fig.savefig(output_dir / f"{kind}_multiplicity.png", dpi=150)
    plt.close(fig)


def plot_total_isr_energy(output_dir: Path, values_by_sample: dict, styles: dict) -> None:
    maximum = max(max(values_by_sample[sample]["isr"]["total_energy"]) for sample in SAMPLES)
    bins = np.linspace(0.0, 5.0 * np.ceil(maximum / 5.0), 61)
    fig, ax = plt.subplots(figsize=(12, 10))
    for sample in SAMPLES:
        values = np.asarray(values_by_sample[sample]["isr"]["total_energy"])
        color, _ = styles[sample]
        ax.hist(values, bins=bins, weights=np.full(len(values), 1.0 / len(values)), histtype="step", color=color, linewidth=2, label=sample)
    ax.set_yscale("log")
    ax.set_xlabel(r"Total gen ISR $E_\gamma$ per event [GeV]")
    ax.set_ylabel("Fraction of events / bin")
    ax.legend()
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / "isr_total_energy_per_event.png", dpi=150)
    plt.close(fig)


def plot_matching_validation(output_dir: Path, name: str, values_by_sample: dict, styles: dict) -> None:
    all_values = [value for sample in SAMPLES for value in values_by_sample[sample]["isr"][name]]
    lower = 10 ** np.floor(np.log10(min(all_values)))
    if name == "nearest_angle":
        bins = np.geomspace(lower, np.pi, 61)
        xlabel = r"Nearest gen ISR–reco photon opening angle $\alpha$ [rad]"
        reference, label = MAX_ANGLE, rf"Match cut $\alpha < {MAX_ANGLE}$ rad"
    else:
        upper = 10 ** np.ceil(np.log10(max(all_values)))
        bins = np.geomspace(lower, upper, 61)
        xlabel = r"Matched ISR $E_\gamma^{\mathrm{reco}} / E_\gamma^{\mathrm{gen}}$"
        reference, label = 1.0, "Unit response"
    fig, ax = plt.subplots(figsize=(12, 10))
    counts_list = []
    fractions_list = []
    for sample in SAMPLES:
        values = np.asarray(values_by_sample[sample]["isr"][name])
        counts = np.histogram(values, bins=bins)[0]
        counts_list.append(counts)
        fractions_list.append(counts / len(values))
        color, _ = styles[sample]
        legend = f"{sample} (N={len(values)})"
        if name == "nearest_angle":
            legend = f"{sample} ({len(values)}/{len(values_by_sample[sample]['isr']['energy'])} with reco)"
        ax.hist(values, bins=bins, weights=np.full(len(values), 1.0 / len(values)), histtype="step", color=color, linewidth=2, label=legend)
    ax.axvline(reference, color="black", linestyle="--", linewidth=2, label=label)
    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fraction of photons / bin")
    ax.legend()
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    set_hist_yaxis(ax, counts_list, fractions_list)
    fig.tight_layout()
    fig.savefig(output_dir / f"isr_{name}.png", dpi=150)
    plt.close(fig)


def plot_efficiency(output_dir: Path, kind: str, name: str, values_by_sample: dict, styles: dict) -> None:
    bins = ENERGY_BINS if name == "energy" else COS_BINS
    xlabel = r"Gen $E_\gamma$ [GeV]" if name == "energy" else r"Gen $\cos\theta_\gamma$"
    fig, ax = plt.subplots(figsize=(12, 10))
    centers = 0.5 * (bins[:-1] + bins[1:])
    for sample in SAMPLES:
        values = values_by_sample[sample][kind]
        denominator = np.histogram(values[name], bins=bins)[0]
        numerator = np.histogram(values[f"matched_{name}"], bins=bins)[0]
        assert np.all(numerator <= denominator)
        efficiency = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=denominator > 0)
        color, _ = styles[sample]
        ax.step(centers, efficiency, where="mid", color=color, linewidth=2, label=f"{sample} ({sum(numerator)}/{sum(denominator)})")
    if name == "energy":
        ax.set_xscale("log")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Reconstruction efficiency")
    ax.legend()
    ax.grid(alpha=0.3)
    add_delphi_label(ax)
    fig.tight_layout()
    fig.savefig(output_dir / f"{kind}_efficiency_vs_{name}.png", dpi=150)
    plt.close(fig)


def plot_isr_fsr_photons(input_root: Path, output_root: Path) -> None:
    mh.style.use(mh.styles.CMS)
    study_dir = output_root / "isr_fsr_photons"
    directories = {
        "gen": study_dir / "01_gen",
        "reco": study_dir / "02_reco",
        "matching": study_dir / "03_matching",
        "efficiency": study_dir / "04_efficiency",
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    styles = sample_styles(SAMPLES)
    values_by_sample = {sample: read_sample(input_root, sample) for sample in SAMPLES}
    for kind in POPULATIONS:
        group = "matching" if kind.endswith("matched_reco") else "reco" if kind == "reco" else "gen"
        plot_multiplicity(directories[group], kind, values_by_sample, styles)
        for name in ("energy", "cos_theta", "phi"):
            plot_distribution(directories[group], kind, name, values_by_sample, styles)
        if kind in ("isr", "fsr"):
            for name in ("energy", "cos_theta"):
                plot_efficiency(directories["efficiency"], kind, name, values_by_sample, styles)
    plot_total_isr_energy(directories["gen"], values_by_sample, styles)
    for name in ("nearest_angle", "energy_response"):
        plot_matching_validation(directories["matching"], name, values_by_sample, styles)
    print(f"plots: {study_dir}")
