"""Stable Gen photons, two stored-link directions, and Reco representations."""
from collections import Counter
from pathlib import Path

import numpy as np
import uproot

SAMPLES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")
TARGET_PDG = {"Zee": 11, "Zmumu": 13, "Ztautau": 15, "ZKK": 321, "Zpipi": 211}
SIM_GAMMA_CODE = 21  # Part/SimPart use DELPHI mass codes; GenPart uses PDG IDs.
GEN_P4 = tuple(f"GenPart_vector.fCoordinates.f{axis}" for axis in "XYZT")
RECO_P4 = tuple(f"Photon_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
PART_P4 = tuple(f"Part_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
CONV_P4 = tuple(f"PhotonConv_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
CHANNELS = ("photon_conversion", "all_lineage")
GEN_SELECTIONS = ("gen_gamma_all", "gen_isr_all", "gen_isr_non_beam")
RECO_DEFINITIONS = ("reco_gamma", "reco_gamma_plus_conv", "reco_all")
LINK_STATUSES = ("agreement", "forward_only", "reverse_only", "conflict", "no_origin")
BRANCHES = (
    "GenPart_pdgId", "GenPart_status", "GenPart_parentIdx", "GenPart_simIdx", *GEN_P4,
    "Photon_partIdx", *RECO_P4, "Part_simIdx", "Part_charge", "Part_pdgId", *PART_P4,
    "Part_originVtxIdx", "Part_decayVtxIdx", "Vtx_incomingIdx",
    "SimPart_genIdx", "SimPart_partIdx", "SimPart_pdgId",
    "SimPart_originVtxIdx", "SimVtx_incomingIdx", *CONV_P4,
    "PhotonConv_simPhotonIdx", "PhotonConv_firstDaughterIdx", "PhotonConv_secondDaughterIdx",
)


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
    for depth in range(len(sim_gens)):
        gen = sim_gens[sim]
        if gen >= 0:
            return int(gen), depth
        vertex = sim_vertices[sim]
        if vertex < 0 or vertex_incoming[vertex] < 0:
            return -1, depth
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


def response(gen_idx, gen_p4, origins, object_p4, valid, raw_origins, raw_p4, invalid_gens):
    """Unmatched sums are zero; structurally ambiguous associated energy is NaN."""
    fields = {name: [] for name in ("reco_count", "raw_count", "energy_count", "energy_ratio",
                                    "leading_ratio", "raw_energy_ratio", "energy_valid")}
    for gen, p4 in zip(gen_idx, gen_p4, strict=True):
        selected, raw_selected = origins == gen, raw_origins == gen
        n, raw_n = int(selected.sum()), int(raw_selected.sum())
        usable = (n == 0 and raw_n == 0) or (gen not in invalid_gens and bool(np.all(valid[selected])))
        energy = object_p4[selected, 3]
        values = (n, raw_n, n, float(energy.sum() / p4[3]) if usable else np.nan,
                  float(energy.max() / p4[3]) if n and usable else np.nan,
                  float(raw_p4[raw_selected, 3].sum() / p4[3]), usable)
        for name, value in zip(fields, values, strict=True):
            fields[name].append(value)
    return {name: np.asarray(values) for name, values in fields.items()}


def link_status(forward, reverse):
    origins = forward | reverse
    return ("conflict" if len(origins) > 1 else "no_origin" if not origins else
            "agreement" if forward and reverse else "forward_only" if forward else "reverse_only")


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
    selections = dict(gen_gamma_all=np.ones(len(all_gamma_idx), dtype=bool),
                      gen_isr_all=np.isin(all_gamma_category, ("beam_isr", "nonbeam_isr")),
                      gen_isr_non_beam=all_gamma_category == "nonbeam_isr")
    gen_idx, gen_p4 = all_gamma_idx[selections["gen_isr_all"]], all_gamma_p4[selections["gen_isr_all"]]
    beam = all_gamma_category[selections["gen_isr_all"]] == "beam_isr"
    stable_gamma = set(all_gamma_idx)
    isr = set(gen_idx)
    nonbeam_isr = set(gen_idx[~beam])
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
                                                  raw["SimVtx_incomingIdx"])[0]
                assert resolved[sim] == downward_sim_origins[sim], "Up/down Sim origin disagreement"
            forward[part] = resolved[sim]
    known = [origins | ({int(forward[part])} if forward[part] >= 0 else set())
             for part, origins in enumerate(reverse)]
    part_origins = np.asarray([next(iter(origins)) if len(origins) == 1 else -1 for origins in known], dtype=int)
    validation = {name: Counter() for name in ("all_parts", "ISR_involved", "noBeamISR_involved", "direct_anchors")}
    for part, origins in enumerate(known):
        kind = ("conflict" if len(origins) > 1 else "no_origin" if not origins else
                "agreement" if reverse[part] and forward[part] >= 0 else
                "forward_only" if forward[part] >= 0 else "reverse_only")
        validation["all_parts"][kind] += 1
        if origins & isr:
            validation["ISR_involved"][kind] += 1
        if origins & nonbeam_isr:
            validation["noBeamISR_involved"][kind] += 1
    for gen, sim in enumerate(raw["GenPart_simIdx"]):
        kind = ("missing" if sim < 0 else "consistent" if sim < n_sim and sim_gens[sim] == gen
                else "different_or_invalid_anchor")
        validation["direct_anchors"][kind] += 1
        if gen in isr:
            validation["direct_anchors"][f"ISR_{kind}"] += 1

    part_p4 = np.column_stack([raw[name] for name in PART_P4])
    photon_parts = np.asarray(raw["Photon_partIdx"], dtype=int)
    photon_p4 = np.column_stack([raw[name] for name in RECO_P4])
    assert np.all((photon_parts >= 0) & (photon_parts < n_part))
    assert len(set(photon_parts)) == len(photon_parts)
    assert np.array_equal(photon_p4, part_p4[photon_parts])
    # Only reciprocal saved Reco vertex relations define parent/child structure.
    part_children = [[] for _ in range(n_part)]
    for part, vertex in enumerate(raw["Part_originVtxIdx"]):
        incoming = raw["Vtx_incomingIdx"][vertex] - 1 if vertex >= 0 else -1
        if 0 <= incoming < n_part and raw["Part_decayVtxIdx"][incoming] == vertex:
            part_children[incoming].append(part)
    part_descendants = {part: descendants(part, part_children) for part in range(n_part) if part_children[part]}
    diagnostics = Counter()
    invalid_gens, bad_parts = set(), set()
    for part, children in part_descendants.items():
        group = children | {part}
        anchors = set().union(*(known[p] for p in group))
        if len(anchors) > 1:
            bad_parts.update(group)
            invalid_gens.update(anchors & stable_gamma)
            diagnostics["conflicting_Reco_structural_groups"] += 1
    diagnostics["missing_saved_children"] += sum(vertex >= 0 and not part_children[part]
        for part, vertex in enumerate(raw["Part_decayVtxIdx"]))
    removed_parts = set()
    for part, children in part_descendants.items():
        anchors = set().union(*(known[p] for p in children | {part}))
        if part_origins[part] >= 0 and anchors == {int(part_origins[part])}:
            removed_parts.update(children)
    canonical_parts = np.asarray([part for part in range(n_part) if part not in removed_parts], dtype=int)
    diagnostics["all_lineage_removed_descendant_Parts"] = len(removed_parts)
    all_valid = np.asarray([part not in bad_parts and len(known[part]) <= 1 for part in canonical_parts], dtype=bool)
    # An unresolved retained parent may overlap its retained daughters geometrically.
    # Do not invent its Gen origin: nominal child response stays valid, but a cone
    # containing that parent has no unambiguous total reconstructed energy.
    all_geometry_valid = all_valid.copy()
    retained_parts = set(canonical_parts)
    for row, part in enumerate(canonical_parts):
        if part_origins[part] < 0 and part_descendants.get(part, set()) & retained_parts:
            all_geometry_valid[row] = False
            diagnostics["all_lineage_unresolved_parent_geometry_overlap"] += 1

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
        diagnostics["conversion_rows"] += 1
        diagnostics["conversion_direct_sim_index_nonnegative"] += raw["PhotonConv_simPhotonIdx"][conv] >= 0
        diagnostics["conversion_missing_parent_topology"] += parent < 0
        diagnostics["conversion_conflicting_origins"] += len(anchors) > 1
        diagnostics["conversion_conflict_involving_ISR"] += len(anchors) > 1 and bool(anchors & isr)
        diagnostics["conversion_unique_ISR"] += conv_origins[conv] in isr

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
        diagnostics["duplicate_conversion_parent_conflicting_origins"] += len(rows) > 1 and len(anchors) > 1

    # Explicit conversions represent their parent and saved descendants once.
    kept_convs, seen_parents, removed_photons = [], set(), set()
    photon_invalid_gens = set()
    for conv, parent in enumerate(conv_parents):
        if parent >= 0 and parent in seen_parents:
            diagnostics["duplicate_conversion_rows_same_parent"] += 1
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
    diagnostics["removed_nested_conversion_representations"] = len(nested_convs)
    kept_photons = np.asarray([r for r in range(len(photon_parts)) if r not in removed_photons], dtype=int)
    kept_convs = np.asarray(kept_convs, dtype=int)
    gamma_p4 = np.concatenate((photon_p4[kept_photons], conv_p4[kept_convs]))
    gamma_origins = np.concatenate((part_origins[photon_parts[kept_photons]], conv_origins[kept_convs]))
    gamma_valid = np.asarray([part not in bad_parts and len(known[part]) <= 1
                             for part in photon_parts[kept_photons]] +
                            [conv_parents[c] >= 0 and len(conv_known[c]) <= 1 and
                             conv_parents[c] not in bad_parts for c in kept_convs], dtype=bool)
    diagnostics["photon_conversion_removed_Photon_representations"] = len(removed_photons)
    raw_gamma_p4 = np.concatenate((photon_p4, conv_p4))
    raw_gamma_origins = np.concatenate((part_origins[photon_parts], conv_origins))
    part_forward = [{int(gen)} if gen >= 0 else set() for gen in forward]
    photon_valid = np.asarray([part not in bad_parts and len(known[part]) <= 1
                               for part in photon_parts], dtype=bool)
    gamma_forward = [part_forward[p] for p in photon_parts[kept_photons]] + [conv_forward[c] for c in kept_convs]
    gamma_reverse = [reverse[p] for p in photon_parts[kept_photons]] + [conv_reverse[c] for c in kept_convs]
    gamma_known = [known[p] for p in photon_parts[kept_photons]] + [conv_known[c] for c in kept_convs]
    reco = {}
    for name, p4, origins, valid, raw_p4, raw_origins, ambiguous, forward_sets, reverse_sets, known_sets, kind in (
        ("reco_gamma", photon_p4, part_origins[photon_parts], photon_valid, photon_p4,
         part_origins[photon_parts], invalid_gens, [part_forward[p] for p in photon_parts],
         [reverse[p] for p in photon_parts], [known[p] for p in photon_parts],
         np.full(len(photon_parts), "gamma")),
        ("reco_gamma_plus_conv", gamma_p4, gamma_origins, gamma_valid, raw_gamma_p4, raw_gamma_origins,
         photon_invalid_gens, gamma_forward, gamma_reverse, gamma_known,
         np.asarray(["gamma"] * len(kept_photons) + ["gamma_conv"] * len(kept_convs))),
        ("reco_all", part_p4[canonical_parts], part_origins[canonical_parts], all_valid,
         part_p4, part_origins, invalid_gens, [part_forward[p] for p in canonical_parts],
         [reverse[p] for p in canonical_parts], [known[p] for p in canonical_parts],
         np.full(len(canonical_parts), "part")),
    ):
        stats = response(all_gamma_idx, all_gamma_p4, origins, p4, valid, raw_origins, raw_p4, ambiguous)
        if name == "reco_all":
            stats["reco_count"] = stats["raw_count"].copy()  # Any stored footprint; energy uses canonical objects.
        geometric_valid = all_geometry_valid if name == "reco_all" else valid
        reco[name] = dict(p4=p4, gen_idx=origins, energy_valid=geometric_valid, stats=stats,
                          forward_origin_sets=forward_sets, reverse_origin_sets=reverse_sets,
                          origin_sets=known_sets, kind=kind)
    # All-Part validation counts raw stored footprints; its energy uses canonical objects.
    reco["reco_all"].update(link_forward_sets=part_forward, link_reverse_sets=reverse, link_origin_sets=known)
    matching_validation = {gen: {name: Counter(dict.fromkeys(LINK_STATUSES, 0)) for name in RECO_DEFINITIONS}
                           for gen in GEN_SELECTIONS}
    unmatched_validation = {name: Counter(dict.fromkeys(LINK_STATUSES, 0)) for name in RECO_DEFINITIONS}
    selected_origins = {name: set(all_gamma_idx[mask]) for name, mask in selections.items()}
    for name, objects in reco.items():
        forward_sets = objects.get("link_forward_sets", objects["forward_origin_sets"])
        reverse_sets = objects.get("link_reverse_sets", objects["reverse_origin_sets"])
        for first, second in zip(forward_sets, reverse_sets, strict=True):
            status = link_status(first, second)
            unmatched_validation[name][status] += 1
            for gen_name, selected in selected_origins.items():
                if (first | second) & selected:
                    matching_validation[gen_name][name][status] += 1
    # Legacy ISR fields remain exact slices for the independent manual-tree audit.
    channels = {}
    for old_name, name in (("photon_conversion", "reco_gamma_plus_conv"), ("all_lineage", "reco_all")):
        objects = reco[name]
        stats = {field: values[selections["gen_isr_all"]] for field, values in objects["stats"].items()}
        diagnostics[f"{old_name}_Gen_with_undefined_energy"] += int(np.count_nonzero(~stats["energy_valid"]))
        channels[old_name] = {field: objects[field] for field in ("p4", "gen_idx", "energy_valid")}
        channels[old_name]["stats"] = stats

    return dict(all_gamma_idx=all_gamma_idx, all_gamma_p4=all_gamma_p4, all_gamma_category=all_gamma_category,
                selections=selections, reco=reco, matching_validation=matching_validation,
                unmatched_validation=unmatched_validation,
                gen_idx=gen_idx, gen_p4=gen_p4, beam=beam, part_origins=part_origins,
                part_origin_sets=known, forward_origins=forward, reverse_origins=reverse,
                sim_origins=downward_sim_origins,
                conversion_origins=conv_origins, conversion_parents=conv_parents,
                conversion_daughters=conv_daughters, conversion_origin_sets=conv_known, channels=channels,
                link_validation=validation, diagnostics=diagnostics)


def read_sample(input_root: Path, sample: str, angle_accumulator=None):
    """Stream files and keep only compact per-Gen observables, not full events."""
    gamma = {name: [] for name in ("energy", "cos_theta", "category")}
    reco = {name: {} for name in RECO_DEFINITIONS}
    reco_distributions = {name: {field: [] for field in ("energy", "cos_theta", "kind", "origin_category")}
                          for name in RECO_DEFINITIONS}
    reco_event_counts = {name: [] for name in RECO_DEFINITIONS}
    matching_validation = {gen: {name: Counter(dict.fromkeys(LINK_STATUSES, 0)) for name in RECO_DEFINITIONS}
                           for gen in GEN_SELECTIONS}
    unmatched_validation = {name: Counter(dict.fromkeys(LINK_STATUSES, 0)) for name in RECO_DEFINITIONS}
    validation = {name: Counter() for name in ("all_parts", "ISR_involved", "noBeamISR_involved", "direct_anchors")}
    diagnostics, events = Counter(), 0
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
                for name, values in (("energy", gamma_p4[:, 3]),
                                     ("cos_theta", unit_momenta(gamma_p4)[:, 2]),
                                     ("category", event["all_gamma_category"])):
                    gamma[name].extend(values)
                origin_categories = dict(zip(event["all_gamma_idx"], event["all_gamma_category"], strict=True))
                for name, objects in event["reco"].items():
                    for field, values in objects["stats"].items():
                        reco[name].setdefault(field, []).extend(values)
                    p4 = objects["p4"]
                    momentum = np.linalg.norm(p4[:, :3], axis=1)
                    cos_theta = np.divide(p4[:, 2], momentum, out=np.full(len(p4), np.nan), where=momentum > 0)
                    category = [origin_categories.get(gen, "other_anchor") if gen >= 0 else "unresolved"
                                for gen in objects["gen_idx"]]
                    for field, values in (("energy", p4[:, 3]), ("cos_theta", cos_theta),
                                          ("kind", objects["kind"]), ("origin_category", category)):
                        reco_distributions[name][field].extend(values)
                    reco_event_counts[name].append(len(event["part_origins"]) if name == "reco_all" else len(p4))
                    unmatched_validation[name].update(event["unmatched_validation"][name])
                    for gen_name in GEN_SELECTIONS:
                        matching_validation[gen_name][name].update(event["matching_validation"][gen_name][name])
                for name, values in event["link_validation"].items():
                    validation[name].update(values)
                diagnostics.update(event["diagnostics"])
                if angle_accumulator is not None:
                    angle_accumulator(event)
        print(f"{sample} {path.parent.name}: {events} events", flush=True)
    gamma = {name: np.asarray(values) for name, values in gamma.items()}
    reco = {name: {field: np.asarray(values) for field, values in fields.items()}
            for name, fields in reco.items()}
    reco_distributions = {name: {field: np.asarray(values) for field, values in fields.items()}
                          for name, fields in reco_distributions.items()}
    reco_event_counts = {name: np.asarray(values) for name, values in reco_event_counts.items()}
    for name, fields in reco.items():
        unmatched = fields["reco_count"] == 0
        assert np.all(fields["energy_ratio"][unmatched & fields["energy_valid"]] == 0)
        assert np.all(fields["energy_ratio"][~fields["energy_valid"]] != fields["energy_ratio"][~fields["energy_valid"]])
        print(f"{sample} {name}: Gen={len(gamma['energy'])}, with Reco={np.count_nonzero(~unmatched)}, "
              f"raw objects={fields['raw_count'].sum()}, canonical objects={fields['energy_count'].sum()}", flush=True)
    selections = dict(gen_gamma_all=np.ones(len(gamma["energy"]), dtype=bool),
                      gen_isr_all=np.isin(gamma["category"], ("beam_isr", "nonbeam_isr")),
                      gen_isr_non_beam=gamma["category"] == "nonbeam_isr")
    return dict(events=events, gamma=gamma, selections=selections, reco=reco,
                reco_distributions=reco_distributions, reco_event_counts=reco_event_counts,
                matching_validation={gen: {name: dict(counts) for name, counts in definitions.items()}
                                     for gen, definitions in matching_validation.items()},
                unmatched_validation={name: dict(counts) for name, counts in unmatched_validation.items()},
                link_validation={name: dict(values) for name, values in validation.items()}, diagnostics=dict(diagnostics))
