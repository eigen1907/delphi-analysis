"""Stable photon origins and stored Gen–Sim–Reco ancestry for Florian samples."""
from pathlib import Path

import numpy as np
import uproot

SAMPLES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")
POPULATIONS = ("01_gamma", "02_ISR", "03_noBeamISR")
TARGET_PDG = {"Zee": 11, "Zmumu": 13, "Ztautau": 15, "ZKK": 321, "Zpipi": 211}
SIM_GAMMA_CODE = 21  # SimPart_pdgId stores DELPHI mass codes, not PDG IDs.
GEN_P4 = tuple(f"GenPart_vector.fCoordinates.f{axis}" for axis in "XYZT")
RECO_P4 = tuple(f"Photon_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
BRANCHES = (
    "GenPart_pdgId", "GenPart_status", "GenPart_parentIdx", *GEN_P4, *RECO_P4,
    "Photon_partIdx", "Part_simIdx", "SimPart_genIdx", "SimPart_pdgId",
    "SimPart_originVtxIdx", "SimVtx_incomingIdx",
)


def unit_momenta(p4):
    momentum = np.linalg.norm(p4[:, :3], axis=1)
    if not np.all(np.isfinite(p4)) or np.any(momentum == 0) or np.any(p4[:, 3] <= 0):
        raise ValueError("Photon has nonfinite, zero-momentum, or nonpositive-energy four-vector")
    return p4[:, :3] / momentum[:, None]


def photon_origin(index, pdgs, parents, hard_parent, target_pdg):
    """Keep the existing ISR/FSR ancestry definition."""
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
    """Trace any Sim species to its first Gen anchor; return origin and depth."""
    for depth in range(len(sim_gens)):
        gen = sim_gens[sim]
        if gen >= 0:
            return gen, depth
        vertex = sim_vertices[sim]
        if vertex < 0 or vertex_incoming[vertex] < 0:
            return -1, depth
        sim = vertex_incoming[vertex]
    raise ValueError("Cycle in Sim ancestry")


def read_sample(input_root: Path, sample: str):
    # All fields below are aligned: exactly one row per stable Gen photon.
    photons = {name: [] for name in (
        "energy", "cos_theta", "isr", "beam", "sim_count", "descendant_count",
        "sim_gamma_count", "reco_count", "leading_ratio", "summed_ratio", "opening_angle",
    )}
    association = dict(reco_candidates=0, no_sim_link=0, unresolved_sim_origin=0)
    n_events = 0
    directory = input_root / f"20260828_100kTest_{sample}_photosFSR" / "final_root"
    paths = sorted(directory.glob("job_*/nanoaod.root"))
    if not paths:
        raise FileNotFoundError(directory)
    for path in paths:
        with uproot.open(path) as root_file:
            tree = root_file["Events"]
            branches = {name: tree[name].array(library="ak").to_list() for name in BRANCHES}
            n_events += tree.num_entries
            for event in range(tree.num_entries):
                pdgs = np.asarray(branches["GenPart_pdgId"][event])
                statuses = np.asarray(branches["GenPart_status"][event])
                parents = np.asarray(branches["GenPart_parentIdx"][event])
                target_status = 21 if sample == "Ztautau" else 1
                daughters = np.flatnonzero((np.abs(pdgs) == TARGET_PDG[sample]) & (statuses == target_status))
                if len(daughters) != 2 or parents[daughters[0]] != parents[daughters[1]]:
                    raise ValueError(f"Expected two daughters with a shared parent: {path}, entry {event}")
                hard_parent = parents[daughters[0]]
                gen_indices = np.flatnonzero((pdgs == 22) & (statuses == 1))
                gen_p4 = np.column_stack([branches[name][event] for name in GEN_P4])[gen_indices]
                reco_p4 = np.column_stack([branches[name][event] for name in RECO_P4])
                gen_unit, reco_unit = unit_momenta(gen_p4), unit_momenta(reco_p4)
                sim_gens = np.asarray(branches["SimPart_genIdx"][event])
                sim_codes = np.asarray(branches["SimPart_pdgId"][event])
                sim_links, sim_depths = np.asarray([
                    first_gen_ancestor(index, sim_gens, branches["SimPart_originVtxIdx"][event],
                                       branches["SimVtx_incomingIdx"][event])
                    for index in range(len(sim_gens))
                ], dtype=int).reshape(-1, 2).T
                assert np.all(sim_links < len(pdgs)), "Sim origin exceeds Gen collection"
                assert np.all(sim_codes[np.isin(sim_gens, gen_indices)] == SIM_GAMMA_CODE), "Unexpected photon mass code"
                reco_links = np.full(len(reco_p4), -1, dtype=int)
                for r, part in enumerate(branches["Photon_partIdx"][event]):
                    sim = branches["Part_simIdx"][event][part] if part >= 0 else -1
                    if sim < 0:
                        association["no_sim_link"] += 1
                    else:
                        reco_links[r] = sim_links[sim]
                        association["unresolved_sim_origin"] += int(reco_links[r] < 0)
                association["reco_candidates"] += len(reco_p4)
                for row, g in enumerate(gen_indices):
                    isr = photon_origin(g, pdgs, parents, hard_parent, TARGET_PDG[sample]) == "isr"
                    parent = parents[g]
                    beam = (isr and gen_p4[row, 0] == 0 and gen_p4[row, 1] == 0
                            and 0 <= parent < len(pdgs) and abs(pdgs[parent]) == 11
                            and parents[parent] < 0)
                    lineage = sim_links == g
                    secondary = lineage & (sim_depths > 0)  # Excludes the directly Gen-linked anchor.
                    linked_reco = np.flatnonzero(reco_links == g)
                    n_sim, n_reco = np.count_nonzero(lineage), len(linked_reco)
                    assert n_reco == 0 or n_sim > 0
                    leading_ratio, summed_ratio, opening_angle = np.nan, np.nan, np.nan
                    if n_reco:
                        leading = linked_reco[np.argmax(reco_p4[linked_reco, 3])]
                        leading_ratio = reco_p4[leading, 3] / gen_p4[row, 3]
                        summed_ratio = reco_p4[linked_reco, 3].sum() / gen_p4[row, 3]
                        opening_angle = np.arccos(np.clip(gen_unit[row] @ reco_unit[leading], -1, 1))
                    for name, value in (
                        ("energy", gen_p4[row, 3]), ("cos_theta", gen_unit[row, 2]),
                        ("isr", isr), ("beam", beam), ("sim_count", n_sim),
                        ("descendant_count", np.count_nonzero(secondary)),
                        ("sim_gamma_count", np.count_nonzero(secondary & (sim_codes == SIM_GAMMA_CODE))),
                        ("reco_count", n_reco), ("leading_ratio", leading_ratio),
                        ("summed_ratio", summed_ratio), ("opening_angle", opening_angle),
                    ):
                        photons[name].append(value)
    photons = {name: np.asarray(entries) for name, entries in photons.items()}
    # Operational topology: saved secondary particles, not a process-ID/shower claim.
    photons["topology"] = np.where(photons["sim_count"] == 0, 0,
                                    np.where(photons["descendant_count"] > 0, 2, 1))
    selections = (np.ones(len(photons["energy"]), dtype=bool), photons["isr"],
                  photons["isr"] & ~photons["beam"])
    populations = {name: {field: entries[selected] for field, entries in photons.items()}
                   for name, selected in zip(POPULATIONS, selections, strict=True)}
    for name, pop in populations.items():
        topology = np.bincount(pop["topology"], minlength=3)
        n_sim, n_reco = np.count_nonzero(pop["sim_count"]), np.count_nonzero(pop["reco_count"])
        assert topology.sum() == len(pop["energy"]) and n_reco <= n_sim <= len(pop["energy"])
        matched = pop["reco_count"] > 0
        assert np.all(pop["sim_gamma_count"] <= pop["descendant_count"])
        assert np.all(pop["summed_ratio"][matched] >= pop["leading_ratio"][matched])
        print(f"{sample} {name}: Gen={len(pop['energy'])}, no_sim/direct/shower="
              f"{'/'.join(map(str, topology))}, GS={n_sim}, SR={n_reco}/{n_sim}, GR={n_reco}", flush=True)
    return dict(events=n_events, populations=populations, association=association)
