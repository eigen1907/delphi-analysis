"""Photon definitions and one shared gen–reco match for the Florian samples."""
from pathlib import Path

import numpy as np
import uproot

SAMPLES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")
TARGET_PDG = {"Zee": 11, "Zmumu": 13, "Ztautau": 15, "ZKK": 321, "Zpipi": 211}
MAX_ANGLE = 0.05  # Opening angle in radians; no energy requirement in matching.
POPULATIONS = (
    "stable_gen", "reco", "isr", "collinear_isr", "noncollinear_isr", "fsr", "others",
    "matched_gen", "matched_reco", "isr_matched_gen", "isr_matched_reco",
    "truth_matched_gen", "isr_truth_matched_gen",
)
GEN_P4 = tuple(f"GenPart_vector.fCoordinates.f{axis}" for axis in "XYZT")
RECO_P4 = tuple(f"Photon_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
BRANCHES = (
    "GenPart_pdgId", "GenPart_status", "GenPart_parentIdx", *GEN_P4, *RECO_P4,
    "Photon_partIdx", "Part_simIdx", "SimPart_genIdx", "Event_cmEnergy",
)


def kinematics(p4):
    """Return parallel (energy, cos(theta), phi) rows and unit momentum vectors.

    Undefined azimuth is retained as NaN. Invalid four-vectors stop the study;
    they are never silently removed from a population or efficiency denominator.
    """
    momentum = np.linalg.norm(p4[:, :3], axis=1)
    if not np.all(np.isfinite(p4)) or np.any(momentum == 0) or np.any(p4[:, 3] <= 0):
        raise ValueError("Photon has nonfinite, zero-momentum, or nonpositive-energy four-vector")
    unit = p4[:, :3] / momentum[:, None]
    phi = np.mod(np.arctan2(p4[:, 1], p4[:, 0]), 2 * np.pi)
    phi[(p4[:, 0] == 0) & (p4[:, 1] == 0)] = np.nan
    return np.column_stack((p4[:, 3], unit[:, 2], phi)), unit


def photon_origin(index, pdgs, parents, hard_parent, target_pdg):
    ancestors = []
    current = parents[index]
    while 0 <= current < len(pdgs) and current not in ancestors:
        ancestors.append(current)
        current = parents[current]
    emitter = next((index for index in ancestors if pdgs[index] != 22), None)
    if emitter is None:
        return "others"
    if emitter == hard_parent or (parents[emitter] == hard_parent and abs(pdgs[emitter]) == target_pdg):
        return "fsr"
    if hard_parent not in ancestors and abs(pdgs[emitter]) == 11:
        return "isr"
    return "others"


def linked_gen_index(reco_index, photon_parts, part_sims, sim_gens):
    part = photon_parts[reco_index]
    if not 0 <= part < len(part_sims):
        return -1
    sim = part_sims[part]
    if not 0 <= sim < len(sim_gens):
        return -1
    return sim_gens[sim]


def match_photons(gen_unit, reco_unit, gen_indices, reco_links):
    """Closest-angle-first, one to one; known disagreeing truth links veto pairs."""
    angles = np.arccos(np.clip(gen_unit @ reco_unit.T, -1.0, 1.0))
    gen_rows, reco_rows = np.where(angles < MAX_ANGLE)
    candidates = sorted(
        (angles[g, r], g, r) for g, r in zip(gen_rows, reco_rows, strict=True)
        if reco_links[r] < 0 or reco_links[r] == gen_indices[g]
    )
    matched = {}
    used_reco = set()
    for _, gen, reco in candidates:
        if gen not in matched and reco not in used_reco:
            matched[gen] = reco
            used_reco.add(reco)
    return matched, angles


def record(values, photons):
    values["multiplicity"].append(len(photons))
    for column, name in enumerate(("energy", "cos_theta", "phi")):
        values[name].extend(photons[:, column])


def read_sample(input_root: Path, sample: str):
    values = {kind: {name: [] for name in ("multiplicity", "energy", "cos_theta", "phi")} for kind in POPULATIONS}
    values["events"] = {name: [] for name in ("isr_energy", "collinear_energy", "matched_reco_energy", "truth_reco_energy", "cm_energy")}
    values["matching"] = {name: [] for name in ("nearest_angle", "energy_response")}
    n_linked = 0
    n_no_reco = 0
    directory = input_root / f"20260828_100kTest_{sample}_photosFSR" / "final_root"
    paths = sorted(directory.glob("job_*/nanoaod.root"))
    if not paths:
        raise FileNotFoundError(directory)
    for path in paths:
        with uproot.open(path) as root_file:
            tree = root_file["Events"]
            branches = {name: tree[name].array(library="ak").to_list() for name in BRANCHES}
            for event in range(tree.num_entries):
                pdgs = np.asarray(branches["GenPart_pdgId"][event])
                statuses = np.asarray(branches["GenPart_status"][event])
                parents = np.asarray(branches["GenPart_parentIdx"][event])
                target_status = 21 if sample == "Ztautau" else 1
                daughters = np.flatnonzero((np.abs(pdgs) == TARGET_PDG[sample]) & (statuses == target_status))
                if len(daughters) != 2 or parents[daughters[0]] != parents[daughters[1]]:
                    raise ValueError(f"Expected the sample's two daughters to share a parent: {path}, entry {event}")
                hard_parent = parents[daughters[0]]
                gen_indices = np.flatnonzero((pdgs == 22) & (statuses == 1))
                gen_p4 = np.column_stack([branches[name][event] for name in GEN_P4])[gen_indices]
                reco_p4 = np.column_stack([branches[name][event] for name in RECO_P4])
                gen, gen_unit = kinematics(gen_p4)
                reco, reco_unit = kinematics(reco_p4)
                origin = np.array([
                    photon_origin(index, pdgs, parents, hard_parent, TARGET_PDG[sample])
                    for index in gen_indices
                ])
                isr, fsr = origin == "isr", origin == "fsr"
                collinear = np.array([
                    isr[row] and gen_p4[row, 0] == 0 and gen_p4[row, 1] == 0
                    and 0 <= parents[index] < len(pdgs)
                    and abs(pdgs[parents[index]]) == 11 and parents[parents[index]] < 0
                    for row, index in enumerate(gen_indices)
                ], dtype=bool)
                others = origin == "others"
                reco_links = np.array([
                    linked_gen_index(index, branches["Photon_partIdx"][event],
                                     branches["Part_simIdx"][event], branches["SimPart_genIdx"][event])
                    for index in range(len(reco))
                ], dtype=int)
                matched, angles = match_photons(gen_unit, reco_unit, gen_indices, reco_links)
                gen_rows = np.array(list(matched), dtype=int)
                reco_rows = np.array(list(matched.values()), dtype=int)
                isr_matched = isr[gen_rows]
                truth_matched = reco_links[reco_rows] >= 0
                for kind, rows in (
                    ("stable_gen", gen), ("reco", reco), ("isr", gen[isr]),
                    ("collinear_isr", gen[collinear]), ("noncollinear_isr", gen[isr & ~collinear]),
                    ("fsr", gen[fsr]), ("others", gen[others]),
                    ("matched_gen", gen[gen_rows]), ("matched_reco", reco[reco_rows]),
                    ("isr_matched_gen", gen[gen_rows[isr_matched]]),
                    ("isr_matched_reco", reco[reco_rows[isr_matched]]),
                    ("truth_matched_gen", gen[gen_rows[truth_matched]]),
                    ("isr_truth_matched_gen", gen[gen_rows[isr_matched & truth_matched]]),
                ):
                    record(values[kind], rows)
                totals = values["events"]
                totals["isr_energy"].append(gen[isr, 0].sum())
                totals["collinear_energy"].append(gen[collinear, 0].sum())
                totals["matched_reco_energy"].append(reco[reco_rows[isr_matched], 0].sum())
                totals["truth_reco_energy"].append(reco[reco_rows[isr_matched & truth_matched], 0].sum())
                totals["cm_energy"].append(branches["Event_cmEnergy"][event])
                if len(reco):
                    values["matching"]["nearest_angle"].extend(angles.min(axis=1))
                else:
                    n_no_reco += len(gen)
                values["matching"]["energy_response"].extend(reco[reco_rows, 0] / gen[gen_rows, 0])
                n_linked += np.count_nonzero(reco_links[reco_rows] >= 0)
    for population in values.values():
        for name in population:
            population[name] = np.asarray(population[name])
    n_events = len(values["events"]["cm_energy"])
    print(
        f"{sample}: events={n_events}, stable={len(values['stable_gen']['energy'])}, "
        f"ISR={len(values['isr']['energy'])}, collinear={len(values['collinear_isr']['energy'])}, "
        f"FSR={len(values['fsr']['energy'])}, others={len(values['others']['energy'])}, "
        f"matched photons={len(values['matched_gen']['energy'])}, matched ISR={len(values['isr_matched_gen']['energy'])}, "
        f"matched with truth link={n_linked}, gen in events with no reco={n_no_reco}", flush=True,
    )
    return values
