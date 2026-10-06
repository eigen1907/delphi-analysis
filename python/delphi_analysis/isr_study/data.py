"""Stable Gen photons and truth-linked Photon or Photon+conversion candidates."""
from collections import Counter
from pathlib import Path

import numpy as np
import uproot

SAMPLES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")
TARGET_PDG = {"Zee": 11, "Zmumu": 13, "Ztautau": 15, "ZKK": 321, "Zpipi": 211}
SIM_GAMMA_CODE = 21  # Part/SimPart use DELPHI mass codes; GenPart uses PDG IDs.
GEN_P4 = tuple(f"GenPart_vector.fCoordinates.f{axis}" for axis in "XYZT")
RECO_P4 = tuple(f"Photon_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
CONV_P4 = tuple(f"PhotonConv_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
GEN_SELECTIONS = ("gen_gamma", "gen_isr", "gen_no_isr", "gen_isr_non_beam")
ANALYSIS_SELECTIONS = ("gen_gamma", "gen_isr", "gen_no_isr")
GEN_LABELS = (r"$\gamma_{\mathrm{all}}$", r"$\gamma_{\mathrm{ISR}}$", r"$\gamma_{\mathrm{no\ ISR}}$")
RECO_LABELS = {"gamma": r"$\gamma$", "gamma_plus_conversion": r"$\gamma+\gamma_{\mathrm{conv}}$"}
LINK_STATUSES = ("agreement", "forward_only", "reverse_only", "conflict", "no_origin")
BRANCHES = (
    "GenPart_pdgId", "GenPart_status", "GenPart_parentIdx", *GEN_P4,
    "Photon_partIdx", *RECO_P4, "Part_simIdx", "Part_charge", "Part_pdgId",
    "Part_originVtxIdx", "Part_decayVtxIdx", "Vtx_incomingIdx",
    "SimPart_genIdx", "SimPart_partIdx",
    "SimPart_originVtxIdx", "SimVtx_incomingIdx", *CONV_P4,
    "PhotonConv_firstDaughterIdx", "PhotonConv_secondDaughterIdx",
)


def gen_selections(category):
    isr = np.isin(category, ("beam_isr", "nonbeam_isr"))
    return dict(gen_gamma=np.ones(len(category), dtype=bool), gen_isr=isr,
                gen_no_isr=~isr, gen_isr_non_beam=category == "nonbeam_isr")


def analysis_selections(category):
    """Stages 02–04 exclude the beam-collinear generator component."""
    return dict(gen_gamma=category != "beam_isr", gen_isr=category == "nonbeam_isr",
                gen_no_isr=np.isin(category, ("fsr", "decayed")))


def unit_momenta(p4):
    momentum = np.linalg.norm(p4[:, :3], axis=1)
    if not np.all(np.isfinite(p4)) or np.any(momentum == 0) or np.any(p4[:, 3] <= 0):
        raise ValueError("Gen photon has nonfinite, zero-momentum, or nonpositive-energy four-vector")
    return p4[:, :3] / momentum[:, None]


def photon_origin(index, pdgs, parents, hard_parent, target_pdg):
    """The existing generator-ancestry ISR/FSR definition."""
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
    """Stop at the first saved Gen anchor, including non-photon anchors."""
    for _ in range(len(sim_gens)):
        gen = sim_gens[sim]
        if gen >= 0:
            return int(gen)
        vertex = sim_vertices[sim]
        if vertex < 0 or vertex_incoming[vertex] < 0:
            return -1
        sim = vertex_incoming[vertex]
    raise ValueError("Cycle in Sim ancestry")


def descendants(index, children):
    """Saved descendants, with no inference from angle or energy."""
    found = set()
    pending = list(children[index])
    while pending:
        child = pending.pop()
        if child == index or child in found:
            raise ValueError("Cycle or duplicate edge in saved ancestry")
        found.add(child)
        pending.extend(children[child])
    return found


def response(gen_idx, gen_p4, origins, object_p4, valid, raw_origins, invalid_gens):
    """Each Gen photon counts once; ambiguous associated energy remains NaN."""
    counts, ratios, defined = [], [], []
    for gen, p4 in zip(gen_idx, gen_p4, strict=True):
        selected = origins == gen
        count = int(selected.sum())
        unmatched = count == 0 and not np.any(raw_origins == gen)
        usable = unmatched or (gen not in invalid_gens and np.all(valid[selected]))
        counts.append(count)
        ratios.append(object_p4[selected, 3].sum() / p4[3] if usable else np.nan)
        defined.append(usable)
    return dict(reco_count=np.asarray(counts), energy_ratio=np.asarray(ratios),
                energy_valid=np.asarray(defined))


def link_status(forward, reverse):
    if len(forward | reverse) > 1:
        return "conflict"
    if forward and reverse:
        return "agreement"
    if forward:
        return "forward_only"
    if reverse:
        return "reverse_only"
    return "no_origin"


def analyze_event(raw, sample):
    """Resolve all Gen origins before selecting the stable-photon populations."""
    pdgs, statuses, parents = (np.asarray(raw[name], dtype=int) for name in
                               ("GenPart_pdgId", "GenPart_status", "GenPart_parentIdx"))
    daughters = np.flatnonzero((np.abs(pdgs) == TARGET_PDG[sample]) &
                              (statuses == (21 if sample == "Ztautau" else 1)))
    assert len(daughters) == 2 and parents[daughters[0]] == parents[daughters[1]]
    hard_parent = parents[daughters[0]]
    all_gamma_idx = np.flatnonzero((pdgs == 22) & (statuses == 1))
    all_gamma_p4 = np.column_stack([raw[name] for name in GEN_P4])[all_gamma_idx]
    unit_momenta(all_gamma_p4)
    all_gamma_category = []
    for gen, p4 in zip(all_gamma_idx, all_gamma_p4, strict=True):
        origin = photon_origin(gen, pdgs, parents, hard_parent, TARGET_PDG[sample])
        is_beam = (origin == "isr" and p4[0] == p4[1] == 0 and 0 <= parents[gen] < len(pdgs)
                   and abs(pdgs[parents[gen]]) == 11 and parents[parents[gen]] < 0)
        all_gamma_category.append("beam_isr" if is_beam else "nonbeam_isr" if origin == "isr"
                                  else "fsr" if origin == "fsr" else "decayed")
    all_gamma_category = np.asarray(all_gamma_category)
    selections = gen_selections(all_gamma_category)
    analysis = analysis_selections(all_gamma_category)
    stable_gamma = set(all_gamma_idx)
    sim_gens = np.asarray(raw["SimPart_genIdx"], dtype=int)
    n_sim, n_part = len(sim_gens), len(raw["Part_simIdx"])
    assert np.all((sim_gens >= -1) & (sim_gens < len(pdgs)))
    sim_children = [[] for _ in range(n_sim)]
    for sim, vertex in enumerate(raw["SimPart_originVtxIdx"]):
        incoming = raw["SimVtx_incomingIdx"][vertex] if vertex >= 0 else -1
        if incoming >= 0:
            assert incoming < n_sim
            sim_children[incoming].append(sim)

    # A is an explicit downward traversal, independent of the Part->Sim traversal B.
    reverse = [set() for _ in range(n_part)]
    downward_sim_origins = np.full(n_sim, -1, dtype=int)
    for anchor in np.flatnonzero(sim_gens >= 0):
        gen, pending, seen = int(sim_gens[anchor]), [int(anchor)], set()
        while pending:
            sim = pending.pop()
            if sim in seen:
                raise ValueError("Cycle in downward Sim ancestry")
            seen.add(sim)
            if sim_gens[sim] >= 0 and sim_gens[sim] != gen:
                continue
            assert downward_sim_origins[sim] in (-1, gen)
            downward_sim_origins[sim] = gen
            part = raw["SimPart_partIdx"][sim]
            if part >= 0:
                assert part < n_part
                reverse[part].add(gen)
            pending.extend(sim_children[sim])
    forward = np.full(n_part, -1, dtype=int)
    resolved = {}
    for part, sim in enumerate(raw["Part_simIdx"]):
        if sim >= 0:
            assert sim < n_sim
            if sim not in resolved:
                resolved[sim] = first_gen_ancestor(sim, sim_gens, raw["SimPart_originVtxIdx"],
                                                  raw["SimVtx_incomingIdx"])
                assert resolved[sim] == downward_sim_origins[sim], "Up/down Sim origin disagreement"
            forward[part] = resolved[sim]
    known = [origins | ({int(forward[part])} if forward[part] >= 0 else set())
             for part, origins in enumerate(reverse)]
    part_origins = np.asarray([next(iter(origins)) if len(origins) == 1 else -1 for origins in known], dtype=int)
    photon_parts = np.asarray(raw["Photon_partIdx"], dtype=int)
    photon_p4 = np.column_stack([raw[name] for name in RECO_P4])
    assert np.all((photon_parts >= 0) & (photon_parts < n_part))
    assert len(set(photon_parts)) == len(photon_parts)
    # Only reciprocal saved Reco vertex relations define parent/child structure.
    part_children = [[] for _ in range(n_part)]
    for part, vertex in enumerate(raw["Part_originVtxIdx"]):
        incoming = raw["Vtx_incomingIdx"][vertex] - 1 if vertex >= 0 else -1
        if 0 <= incoming < n_part and raw["Part_decayVtxIdx"][incoming] == vertex:
            part_children[incoming].append(part)
    part_descendants = {part: descendants(part, part_children) for part in range(n_part) if part_children[part]}
    bad_parts = set()
    for part, children in part_descendants.items():
        group = children | {part}
        anchors = set().union(*(known[p] for p in group))
        if len(anchors) > 1:
            bad_parts.update(group)
    conv_p4 = np.column_stack([raw[name] for name in CONV_P4])
    conv_origins = np.full(len(conv_p4), -1, dtype=int)
    conv_parents, conv_daughters, conv_known, conv_forward, conv_reverse = [], [], [], [], []
    for conv in range(len(conv_p4)):
        daughters = {raw[name][conv] for name in ("PhotonConv_firstDaughterIdx", "PhotonConv_secondDaughterIdx")
                     if raw[name][conv] >= 0}
        assert all(part < n_part for part in daughters)
        vertices = {raw["Part_originVtxIdx"][part] for part in daughters}
        parent = -1
        if len(vertices) == 1:
            vertex = next(iter(vertices))
            candidate = raw["Vtx_incomingIdx"][vertex] - 1 if vertex >= 0 else -1
            if (0 <= candidate < n_part and raw["Part_decayVtxIdx"][candidate] == vertex
                    and raw["Part_charge"][candidate] == 0 and raw["Part_pdgId"][candidate] == SIM_GAMMA_CODE):
                parent = candidate
        members = daughters | ({parent} if parent >= 0 else set())
        anchors = set().union(*(known[p] for p in members))
        if parent >= 0 and len(anchors) == 1:
            conv_origins[conv] = next(iter(anchors))
        conv_parents.append(parent)
        conv_daughters.append(daughters)
        conv_known.append(anchors)
        conv_forward.append({int(forward[p]) for p in members if forward[p] >= 0})
        conv_reverse.append(set().union(*(reverse[p] for p in members)))

    # Duplicate conversion rows are one parent representation; combine their evidence first.
    by_parent = {}
    for conv, parent in enumerate(conv_parents):
        if parent >= 0:
            by_parent.setdefault(parent, []).append(conv)
    for parent, rows in by_parent.items():
        anchors = set().union(*(conv_known[c] for c in rows))
        forward_anchors = set().union(*(conv_forward[c] for c in rows))
        reverse_anchors = set().union(*(conv_reverse[c] for c in rows))
        for conv in rows:
            conv_known[conv] = anchors
            conv_origins[conv] = next(iter(anchors)) if len(anchors) == 1 else -1
            conv_forward[conv], conv_reverse[conv] = forward_anchors, reverse_anchors

    # Explicit conversions represent their parent and saved descendants once.
    kept_convs, seen_parents, removed_photons = [], set(), set()
    photon_invalid_gens = set()
    for conv, parent in enumerate(conv_parents):
        if parent >= 0 and parent in seen_parents:
            continue
        kept_convs.append(conv)
        if parent >= 0:
            seen_parents.add(parent)
            group = {parent} | part_descendants.get(parent, set()) | conv_daughters[conv]
            anchors = set().union(*(known[p] for p in group))
            if len(anchors) > 1:
                photon_invalid_gens.update(anchors & stable_gamma)
            # Conflicting representations remain unresolved, not silently reassigned.
            if len(anchors) <= 1:
                removed_photons.update(r for r, part in enumerate(photon_parts) if part in group
                    and (part_origins[part] < 0 or part_origins[part] == conv_origins[conv]))
    # A nested conversion is also a descendant representation of its accepted ancestor conversion.
    nested_convs = set()
    for conv in kept_convs:
        parent = conv_parents[conv]
        if parent < 0 or conv_origins[conv] < 0:
            continue
        group = {parent} | part_descendants.get(parent, set())
        anchors = set().union(*(known[p] for p in group))
        if anchors == {int(conv_origins[conv])}:
            nested_convs.update(c for c in kept_convs if c != conv and conv_parents[c] in group
                               and conv_origins[c] == conv_origins[conv])
        else:
            photon_invalid_gens.update(anchors & stable_gamma)
    kept_convs = [conv for conv in kept_convs if conv not in nested_convs]
    kept_photons = np.asarray([r for r in range(len(photon_parts)) if r not in removed_photons], dtype=int)
    kept_convs = np.asarray(kept_convs, dtype=int)
    gamma_p4 = np.concatenate((photon_p4[kept_photons], conv_p4[kept_convs]))
    gamma_origins = np.concatenate((part_origins[photon_parts[kept_photons]], conv_origins[kept_convs]))
    gamma_valid = np.asarray([part not in bad_parts and len(known[part]) <= 1
                             for part in photon_parts[kept_photons]] +
                            [conv_parents[c] >= 0 and len(conv_known[c]) <= 1 and
                             conv_parents[c] not in bad_parts for c in kept_convs], dtype=bool)
    raw_gamma_origins = np.concatenate((part_origins[photon_parts], conv_origins))
    # Geometry has no Gen-origin requirement. Only overlapping saved Reco
    # representations make its energy ambiguous; an isolated conflicted Part
    # still has a measured four-vector. Invalidate the retained ancestor, so
    # a cone containing only its daughters can retain a defined energy sum.
    kept_photon_parts = set(photon_parts[kept_photons])
    conversion_groups = []
    for conv in kept_convs:
        members = conv_daughters[conv] | ({conv_parents[conv]} if conv_parents[conv] >= 0 else set())
        conversion_groups.append(members | set().union(*(part_descendants.get(part, set()) for part in members)))
    gamma_geometry_valid = []
    for part in photon_parts[kept_photons]:
        children = part_descendants.get(part, set())
        gamma_geometry_valid.append(not (children & kept_photon_parts) and
                                    not any(children & group for group in conversion_groups))
    for row, group in enumerate(conversion_groups):
        overlap = bool(group & kept_photon_parts) or any(
            group & other and not group < other for other_row, other in enumerate(conversion_groups)
            if other_row != row)
        gamma_geometry_valid.append(not overlap)
    gamma_geometry_valid = np.asarray(gamma_geometry_valid, dtype=bool)
    gamma_forward = [{int(forward[p])} if forward[p] >= 0 else set()
                     for p in photon_parts[kept_photons]] + [conv_forward[c] for c in kept_convs]
    gamma_reverse = [reverse[p] for p in photon_parts[kept_photons]] + [conv_reverse[c] for c in kept_convs]
    combined_stats = response(all_gamma_idx, all_gamma_p4, gamma_origins, gamma_p4,
                              gamma_valid, raw_gamma_origins, photon_invalid_gens)

    # Photon-only keeps every original Photon row, before conversion replacement.
    # Saved ancestor/daughter Photon rows can overlap in energy; keep their counts
    # and mark their energy sum as ambiguous instead of discarding candidates.
    photon_origins = part_origins[photon_parts]
    all_photon_parts = set(photon_parts)
    photon_geometry_valid = np.asarray([
        not (part_descendants.get(part, set()) & all_photon_parts) for part in photon_parts
    ], dtype=bool)
    photon_valid = np.asarray([
        part not in bad_parts and len(known[part]) <= 1 for part in photon_parts
    ], dtype=bool) & photon_geometry_valid
    photon_stats = response(all_gamma_idx, all_gamma_p4, photon_origins, photon_p4,
                            photon_valid, photon_origins, set())
    reco = dict(
        gamma=dict(p4=photon_p4, gen_idx=photon_origins, energy_valid=photon_valid,
                   geometric_energy_valid=photon_geometry_valid, stats=photon_stats),
        gamma_plus_conversion=dict(p4=gamma_p4, gen_idx=gamma_origins, energy_valid=gamma_valid,
                                   geometric_energy_valid=gamma_geometry_valid, stats=combined_stats),
    )
    links = dict(
        gamma=([{int(forward[p])} if forward[p] >= 0 else set() for p in photon_parts],
               [reverse[p] for p in photon_parts]),
        gamma_plus_conversion=(gamma_forward, gamma_reverse),
    )
    matching_validation = {}
    selected_origins = {name: set(all_gamma_idx[mask]) for name, mask in analysis.items()}
    for mode, (forward_links, reverse_links) in links.items():
        counts = {gen: Counter(dict.fromkeys(LINK_STATUSES, 0)) for gen in ANALYSIS_SELECTIONS}
        for first, second in zip(forward_links, reverse_links, strict=True):
            status = link_status(first, second)
            for gen_name, selected in selected_origins.items():
                if (first | second) & selected:
                    counts[gen_name][status] += 1
        matching_validation[mode] = counts
    return dict(all_gamma_idx=all_gamma_idx, all_gamma_p4=all_gamma_p4,
                all_gamma_category=all_gamma_category, selections=selections,
                analysis_selections=analysis, reco=reco,
                matching_validation=matching_validation)


def read_sample(input_root: Path, sample: str, angle_accumulator):
    """Stream files and keep only compact per-Gen observables, not full events."""
    gamma = {name: [] for name in ("energy", "cos_theta", "category")}
    results = {mode: dict(
        reco={field: [] for field in ("reco_count", "energy_ratio", "energy_valid")},
        reco_distributions={field: [] for field in ("energy", "cos_theta")},
        reco_event_counts=[],
        event_truth_energy={gen: [] for gen in ANALYSIS_SELECTIONS},
        matching_validation={gen: Counter(dict.fromkeys(LINK_STATUSES, 0)) for gen in ANALYSIS_SELECTIONS},
    ) for mode in RECO_LABELS}
    event_gen_energy = {gen: [] for gen in ANALYSIS_SELECTIONS}
    gen_event_counts = {name: [] for name in
                        ("all_gamma", "all_isr", "beam_isr", "nonbeam_isr", "fsr", "decayed")}
    events = 0
    directory = input_root / f"20260828_100kTest_{sample}_photosFSR" / "final_root"
    paths = sorted(directory.glob("job_*/nanoaod.root"))
    if not paths:
        raise FileNotFoundError(directory)
    for path in paths:
        with uproot.open(path) as root_file:
            tree = root_file["Events"]
            columns = {name: tree[name].array(library="ak").to_list() for name in BRANCHES}
            events += tree.num_entries
            for entry in range(tree.num_entries):
                event = analyze_event({name: values[entry] for name, values in columns.items()}, sample)
                gamma_p4 = event["all_gamma_p4"]
                categories = event["all_gamma_category"]
                for gen_name, selected in event["analysis_selections"].items():
                    event_gen_energy[gen_name].append(gamma_p4[selected, 3].sum())
                counts = dict(all_gamma=len(categories),
                              all_isr=int(event["selections"]["gen_isr"].sum()))
                counts.update({name: int(np.count_nonzero(categories == name))
                               for name in ("beam_isr", "nonbeam_isr", "fsr", "decayed")})
                for name, count in counts.items():
                    gen_event_counts[name].append(count)
                gamma["energy"].extend(gamma_p4[:, 3])
                gamma["cos_theta"].extend(unit_momenta(gamma_p4)[:, 2])
                gamma["category"].extend(categories)
                for mode, objects in event["reco"].items():
                    values = results[mode]
                    for field, entries in values["reco"].items():
                        entries.extend(objects["stats"][field])
                    p4 = objects["p4"]
                    for gen_name, selected in event["analysis_selections"].items():
                        matched = np.isin(objects["gen_idx"], event["all_gamma_idx"][selected])
                        truth_energy = (p4[matched, 3].sum()
                                        if objects["stats"]["energy_valid"][selected].all() else np.nan)
                        values["event_truth_energy"][gen_name].append(truth_energy)
                    momentum = np.linalg.norm(p4[:, :3], axis=1)
                    cos_theta = np.divide(p4[:, 2], momentum, out=np.full(len(p4), np.nan), where=momentum > 0)
                    values["reco_distributions"]["energy"].extend(p4[:, 3])
                    values["reco_distributions"]["cos_theta"].extend(cos_theta)
                    values["reco_event_counts"].append(len(p4))
                    for gen_name in ANALYSIS_SELECTIONS:
                        values["matching_validation"][gen_name].update(event["matching_validation"][mode][gen_name])
                angle_accumulator(event)
        print(f"{sample} {path.parent.name}: {events} events", flush=True)
    gamma = {name: np.asarray(values) for name, values in gamma.items()}
    gen_event_counts = {name: np.asarray(values) for name, values in gen_event_counts.items()}
    selections = gen_selections(gamma["category"])
    shared = dict(events=events, gamma=gamma, selections=selections,
                  analysis_selections=analysis_selections(gamma["category"]), gen_event_counts=gen_event_counts,
                  event_gen_energy={gen: np.asarray(energy) for gen, energy in event_gen_energy.items()})
    for mode, values in results.items():
        values["reco"] = {field: np.asarray(entries) for field, entries in values["reco"].items()}
        values["reco_distributions"] = {field: np.asarray(entries)
                                        for field, entries in values["reco_distributions"].items()}
        values["reco_event_counts"] = np.asarray(values["reco_event_counts"])
        values["event_truth_energy"] = {gen: np.asarray(energy)
                                        for gen, energy in values["event_truth_energy"].items()}
        values["matching_validation"] = {gen: dict(counts) for gen, counts in values["matching_validation"].items()}
        values.update(shared, reco_label=RECO_LABELS[mode])
        reco = values["reco"]
        print(f"{sample} {mode}: Gen={len(gamma['energy'])}, with Reco={np.count_nonzero(reco['reco_count'])}, "
              f"associated objects={reco['reco_count'].sum()}", flush=True)
    return results
