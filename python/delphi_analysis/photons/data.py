"""Photon origins, Sim ancestry, and geometric diagnostics for the Florian samples."""
from pathlib import Path

import numpy as np
import uproot

SAMPLES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")
TARGET_PDG = {"Zee": 11, "Zmumu": 13, "Ztautau": 15, "ZKK": 321, "Zpipi": 211}
MAX_ANGLE = 0.03  # 3D opening angle in radians; no energy requirement.
SCAN_MAX_ANGLE = 0.10
MATCH_METHODS = ("direct", "truth", "angle", "recovery")
POPULATIONS = (
    "stable_gen", "stable_gen_wo_beam", "reco", "isr", "collinear_isr",
    "noncollinear_isr", "fsr", "others",
    *(f"isr_matched_gen_{method}" for method in MATCH_METHODS),
    *(f"nonbeam_isr_matched_gen_{method}" for method in MATCH_METHODS),
)
GEN_P4 = tuple(f"GenPart_vector.fCoordinates.f{axis}" for axis in "XYZT")
RECO_P4 = tuple(f"Photon_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
SIM_P4 = tuple(f"SimPart_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
BRANCHES = (
    "GenPart_pdgId", "GenPart_status", "GenPart_parentIdx", "GenPart_simIdx", *GEN_P4, *RECO_P4, *SIM_P4,
    "Photon_partIdx", "Part_simIdx", "SimPart_genIdx", "SimPart_originVtxIdx",
    "SimVtx_incomingIdx", "Event_cmEnergy",
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


def first_gen_ancestor(sim, sim_gens, sim_vertices, vertex_incoming):
    """Resolve a saved Sim object to its first Gen ancestor, without a species cut."""
    for depth in range(len(sim_gens)):
        gen = sim_gens[sim]
        if gen >= 0:
            return gen, depth
        vertex = sim_vertices[sim]
        if vertex < 0 or vertex_incoming[vertex] < 0:
            return -1, depth
        sim = vertex_incoming[vertex]
    raise ValueError("Cycle in Sim ancestry")


def linked_gen_indices(reco_index, photon_parts, part_sims, sim_gens, sim_links, sim_depths):
    """Read direct and resolved origins for one reco Photon; missing links are -1."""
    part = photon_parts[reco_index]
    if part < 0:
        return -1, -1, -1
    sim = part_sims[part]
    if sim < 0:
        return -1, -1, -1
    return sim_gens[sim], sim_links[sim], sim_depths[sim]


def angle_matches(angles, initial, reco_allowed, cut):
    """Add closest-angle-first one-to-one pairs to existing matches."""
    matched = initial.copy()
    used_reco = set(matched.values())
    gen_rows, reco_rows = np.where(angles < cut)
    for _, gen, reco in sorted((angles[g, r], g, r) for g, r in zip(gen_rows, reco_rows, strict=True)
                               if reco_allowed[r]):
        if gen not in matched and reco not in used_reco:
            matched[gen] = reco
            used_reco.add(reco)
    return matched


def record(values, photons, sqrt_s):
    values["multiplicity"].append(len(photons))
    for column, name in enumerate(("energy", "cos_theta", "phi")):
        values[name].extend(photons[:, column])
    values["x_gamma"].extend(2 * photons[:, 0] / sqrt_s)


def read_sample(input_root: Path, sample: str):
    values = {kind: {name: [] for name in ("multiplicity", "energy", "cos_theta", "phi", "x_gamma")}
              for kind in POPULATIONS}
    values["angle_scan"] = {name: [] for name in ("truth", "same_link", "no_link", "wrong_link", "beam_no_link")}
    values["response"] = {name: [] for name in (
        "gen_energy", "reco_energy", "opening_angle", "depth", "reco_count",
        "group_gen_energy", "group_reco_energy",
    )}
    values["association_counts"] = {name: [] for name in (
        "no_sim", "unresolved_sim", "other_gen", "direct_isr", "descendant_isr",
    )}
    values["generator"] = {name: [] for name in ("sqrt_s", "metadata_cm_energy")}
    stages = {group: {name: [] for name in (
        "gen_energy", "gen_cos_theta", "gen_phi", "sim_count", "sim_root_count", "reco_count",
        "gen_multiplicity", "sim_multiplicity", "reco_multiplicity",
        "root_sim_energy", "root_sim_cos_theta", "root_sim_phi",
        "reco_energy", "reco_cos_theta", "reco_phi",
    )} for group in ("stable", "isr", "nonbeam")}
    values["sim_bookkeeping"] = {name: [] for name in (
        "missing_forward_link", "forward_reverse_disagreement", "root_p4_max_difference",
    )}
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
                all_gen_p4 = np.column_stack([branches[name][event] for name in GEN_P4])
                assert np.array_equal(pdgs[:2], [11, -11]) and np.all(parents[:2] < 0)
                sqrt_s = all_gen_p4[:2, 3].sum()
                values["generator"]["sqrt_s"].append(sqrt_s)
                values["generator"]["metadata_cm_energy"].append(branches["Event_cmEnergy"][event])
                gen_p4 = all_gen_p4[gen_indices]
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
                sim_gens = np.asarray(branches["SimPart_genIdx"][event])
                sim_links, sim_depths = np.asarray([
                    first_gen_ancestor(index, sim_gens, branches["SimPart_originVtxIdx"][event],
                                       branches["SimVtx_incomingIdx"][event])
                    for index in range(len(sim_gens))
                ], dtype=int).reshape(-1, 2).T
                direct_links, reco_links, depths = np.asarray([
                    linked_gen_indices(index, branches["Photon_partIdx"][event],
                                       branches["Part_simIdx"][event], sim_gens, sim_links, sim_depths)
                    for index in range(len(reco))
                ], dtype=int).reshape(-1, 3).T
                angles = np.arccos(np.clip(gen_unit @ reco_unit.T, -1.0, 1.0))
                gen_by_index = {index: row for row, index in enumerate(gen_indices)}
                direct = {gen_by_index[link]: row for row, link in enumerate(direct_links) if link in gen_by_index}
                truth = {gen_by_index[link]: row for row, link in enumerate(reco_links) if link in gen_by_index}
                angle = angle_matches(angles, {}, np.ones(len(reco), dtype=bool), MAX_ANGLE)
                recovery = angle_matches(angles, truth, reco_links < 0, MAX_ANGLE)
                matched = {"direct": direct, "truth": truth, "angle": angle, "recovery": recovery}
                sim_p4 = np.column_stack([branches[name][event] for name in SIM_P4])
                root_mask = np.isin(sim_gens, gen_indices)
                root_sim, _ = kinematics(sim_p4[root_mask])
                root_gen_indices = sim_gens[root_mask]
                sim_counts = np.array([np.count_nonzero(sim_links == g) for g in gen_indices])
                root_counts = np.array([np.count_nonzero(sim_gens == g) for g in gen_indices])
                reco_counts = np.array([np.count_nonzero(reco_links == g) for g in gen_indices])
                assert np.all((reco_counts == 0) | (sim_counts > 0))
                forward = np.asarray(branches["GenPart_simIdx"][event])[gen_indices]
                present = forward >= 0
                values["sim_bookkeeping"]["missing_forward_link"].append(np.count_nonzero(~present))
                values["sim_bookkeeping"]["forward_reverse_disagreement"].append(
                    np.count_nonzero(sim_gens[forward[present]] != gen_indices[present]))
                values["sim_bookkeeping"]["root_p4_max_difference"].append(
                    np.max(np.abs(sim_p4[root_mask] - all_gen_p4[root_gen_indices]), initial=0))
                for group, selected in (
                    ("stable", np.ones(len(gen), dtype=bool)),
                    ("isr", isr), ("nonbeam", isr & ~collinear),
                ):
                    stage = stages[group]
                    selected_indices = gen_indices[selected]
                    selected_roots = np.isin(root_gen_indices, selected_indices)
                    selected_reco = np.isin(reco_links, selected_indices)
                    for column, name in enumerate(("energy", "cos_theta", "phi")):
                        stage[f"gen_{name}"].extend(gen[selected, column])
                        stage[f"root_sim_{name}"].extend(root_sim[selected_roots, column])
                        stage[f"reco_{name}"].extend(reco[selected_reco, column])
                    stage["sim_count"].extend(sim_counts[selected])
                    stage["sim_root_count"].extend(root_counts[selected])
                    stage["reco_count"].extend(reco_counts[selected])
                    stage["gen_multiplicity"].append(np.count_nonzero(selected))
                    stage["sim_multiplicity"].append(np.count_nonzero(sim_counts[selected]))
                    stage["reco_multiplicity"].append(np.count_nonzero(reco_counts[selected]))
                for kind, rows in (
                    ("stable_gen", gen), ("stable_gen_wo_beam", gen[~collinear]),
                    ("reco", reco), ("isr", gen[isr]),
                    ("collinear_isr", gen[collinear]), ("noncollinear_isr", gen[isr & ~collinear]),
                    ("fsr", gen[fsr]), ("others", gen[others]),
                ):
                    record(values[kind], rows, sqrt_s)
                for method, pairs in matched.items():
                    rows = np.array(list(pairs), dtype=int)
                    record(values[f"isr_matched_gen_{method}"], gen[rows[isr[rows]]], sqrt_s)
                    record(values[f"nonbeam_isr_matched_gen_{method}"], gen[rows[isr[rows] & ~collinear[rows]]], sqrt_s)
                response = values["response"]
                for g in np.flatnonzero(isr):
                    linked_reco = np.flatnonzero(reco_links == gen_indices[g])
                    response["reco_count"].append(len(linked_reco))
                    if len(linked_reco):
                        response["group_gen_energy"].append(gen[g, 0])
                        response["group_reco_energy"].append(reco[linked_reco, 0].sum())
                    for r in linked_reco:
                        response["gen_energy"].append(gen[g, 0])
                        response["reco_energy"].append(reco[r, 0])
                        response["opening_angle"].append(angles[g, r])
                        response["depth"].append(depths[r])
                linked_isr = np.isin(reco_links, gen_indices[isr])
                for name, selection in (
                    ("no_sim", depths < 0), ("unresolved_sim", (depths >= 0) & (reco_links < 0)),
                    ("other_gen", (reco_links >= 0) & ~linked_isr),
                    ("direct_isr", linked_isr & (depths == 0)),
                    ("descendant_isr", linked_isr & (depths > 0)),
                ):
                    values["association_counts"][name].append(np.count_nonzero(selection))
                scan = values["angle_scan"]
                scan["truth"].extend(angles[gen_by_index[reco_links[r]], r] for r in np.flatnonzero(linked_isr))
                # Closest-first pairs below 0.10 rad are unchanged by a tighter cut.
                for g, r in angle_matches(angles, {}, np.ones(len(reco), dtype=bool), SCAN_MAX_ANGLE).items():
                    if not isr[g]:
                        continue
                    if reco_links[r] == gen_indices[g]:
                        name = "same_link"
                    elif reco_links[r] < 0:
                        name = "no_link"
                        if collinear[g]:
                            scan["beam_no_link"].append(angles[g, r])
                    else:
                        name = "wrong_link"
                    scan[name].append(angles[g, r])
    for population in values.values():
        for name in population:
            population[name] = np.asarray(population[name])
    values["stages"] = {
        group: {name: np.asarray(entries) for name, entries in stage.items()}
        for group, stage in stages.items()
    }
    n_events = len(values["stable_gen"]["multiplicity"])
    print(
        f"{sample}: events={n_events}, stable={len(values['stable_gen']['energy'])}, "
        f"ISR={len(values['isr']['energy'])}, collinear={len(values['collinear_isr']['energy'])}, "
        f"FSR={len(values['fsr']['energy'])}, others={len(values['others']['energy'])}, "
        "ISR matched (direct/truth/angle/recovery)="
        + "/".join(str(len(values[f'isr_matched_gen_{method}']['energy'])) for method in MATCH_METHODS),
        flush=True,
    )
    return values
