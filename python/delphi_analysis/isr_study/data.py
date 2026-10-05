"""One stable Gen ISR photon; two stored-link directions and Reco representations."""
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


def analyze_event(raw, sample):
    """Resolve all known Gen origins before selecting the ISR analysis units."""
    pdgs, statuses, parents = (np.asarray(raw[name], dtype=int) for name in
                               ("GenPart_pdgId", "GenPart_status", "GenPart_parentIdx"))
    daughters = np.flatnonzero((np.abs(pdgs) == TARGET_PDG[sample]) &
                              (statuses == (21 if sample == "Ztautau" else 1)))
    assert len(daughters) == 2 and parents[daughters[0]] == parents[daughters[1]]
    hard_parent = parents[daughters[0]]
    gen_idx = np.asarray([gen for gen in np.flatnonzero((pdgs == 22) & (statuses == 1))
                          if photon_origin(gen, pdgs, parents, hard_parent, TARGET_PDG[sample]) == "isr"], dtype=int)
    gen_p4 = np.column_stack([raw[name] for name in GEN_P4])[gen_idx]
    gen_unit = unit_momenta(gen_p4)
    beam = np.asarray([p4[0] == 0 and p4[1] == 0 and 0 <= parents[gen] < len(pdgs)
                       and abs(pdgs[parents[gen]]) == 11 and parents[parents[gen]] < 0
                       for gen, p4 in zip(gen_idx, gen_p4, strict=True)], dtype=bool)
    isr = set(gen_idx)
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
    validation = {name: Counter() for name in ("all_parts", "ISR_involved", "direct_anchors")}
    for part, origins in enumerate(known):
        kind = ("conflict" if len(origins) > 1 else "no_origin" if not origins else
                "agreement" if reverse[part] and forward[part] >= 0 else
                "forward_only" if forward[part] >= 0 else "reverse_only")
        validation["all_parts"][kind] += 1
        if origins & isr:
            validation["ISR_involved"][kind] += 1
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
            invalid_gens.update(anchors & isr)
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
    conv_parents, conv_daughters, conv_known = [], [], []
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
        anchors = set().union(*(known[p] for p in daughters | ({parent} if parent >= 0 else set())))
        if parent >= 0 and len(anchors) == 1:
            conv_origins[conv] = next(iter(anchors))
        conv_parents.append(parent)
        conv_daughters.append(daughters)
        conv_known.append(anchors)
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
        for conv in rows:
            conv_known[conv] = anchors
            conv_origins[conv] = next(iter(anchors)) if len(anchors) == 1 else -1
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
                photon_invalid_gens.update(anchors & isr)
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
            photon_invalid_gens.update(anchors & isr)
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
    channels = {}
    for name, p4, origins, valid, raw_p4, raw_origins, ambiguous in (
        ("photon_conversion", gamma_p4, gamma_origins, gamma_valid, raw_gamma_p4, raw_gamma_origins, photon_invalid_gens),
        ("all_lineage", part_p4[canonical_parts], part_origins[canonical_parts], all_valid,
         part_p4, part_origins, invalid_gens),
    ):
        stats = response(gen_idx, gen_p4, origins, p4, valid, raw_origins, raw_p4, ambiguous)
        if name == "all_lineage":
            stats["reco_count"] = stats["raw_count"].copy()  # Any stored footprint; energy uses canonical objects.
        diagnostics[f"{name}_Gen_with_undefined_energy"] += int(np.count_nonzero(~stats["energy_valid"]))
        geometric_valid = all_geometry_valid if name == "all_lineage" else valid
        channels[name] = dict(p4=p4, gen_idx=origins, energy_valid=geometric_valid, stats=stats)

    gen_rows = {int(gen): row for row, gen in enumerate(gen_idx)}
    composition = {name: [] for name in ("gen_pt", "gen_cos_theta", "beam", "category")}
    for part, gen in enumerate(part_origins):
        if gen not in isr:
            continue
        row = gen_rows[gen]
        charged, code = raw["Part_charge"][part] != 0, raw["Part_pdgId"][part]
        category = ("charged_e" if charged and abs(code) == 2 else "other_charged" if charged
                    else "neutral_gamma" if code == SIM_GAMMA_CODE else "other_neutral")
        for name, value in (("gen_pt", np.hypot(*gen_p4[row, :2])), ("gen_cos_theta", gen_unit[row, 2]),
                            ("beam", beam[row]), ("category", category)):
            composition[name].append(value)
    return dict(gen_idx=gen_idx, gen_p4=gen_p4, beam=beam, part_origins=part_origins,
                part_origin_sets=known, forward_origins=forward, reverse_origins=reverse,
                sim_origins=downward_sim_origins,
                conversion_origins=conv_origins, conversion_parents=conv_parents,
                conversion_daughters=conv_daughters, conversion_origin_sets=conv_known, channels=channels,
                composition=composition, link_validation=validation, diagnostics=diagnostics)


def read_sample(input_root: Path, sample: str, angle_accumulator=None):
    """Stream files and keep only compact per-Gen observables, not full events."""
    photons = {name: [] for name in ("energy", "pt", "cos_theta", "beam")}
    channels = {name: {} for name in CHANNELS}
    composition = {name: [] for name in ("gen_pt", "gen_cos_theta", "beam", "category")}
    validation = {name: Counter() for name in ("all_parts", "ISR_involved", "direct_anchors")}
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
                p4 = event["gen_p4"]
                for name, values in (("energy", p4[:, 3]), ("pt", np.hypot(p4[:, 0], p4[:, 1])),
                                     ("cos_theta", unit_momenta(p4)[:, 2]), ("beam", event["beam"])):
                    photons[name].extend(values)
                for channel in CHANNELS:
                    for name, values in event["channels"][channel]["stats"].items():
                        channels[channel].setdefault(name, []).extend(values)
                for name, values in event["composition"].items():
                    composition[name].extend(values)
                for name, values in event["link_validation"].items():
                    validation[name].update(values)
                diagnostics.update(event["diagnostics"])
                if angle_accumulator is not None:
                    angle_accumulator(event)
        print(f"{sample} {path.parent.name}: {events} events", flush=True)
    photons = {name: np.asarray(values) for name, values in photons.items()}
    channels = {channel: {name: np.asarray(values) for name, values in fields.items()}
                for channel, fields in channels.items()}
    composition = {name: np.asarray(values) for name, values in composition.items()}
    for name, fields in channels.items():
        unmatched = fields["reco_count"] == 0
        assert np.all(fields["energy_ratio"][unmatched & fields["energy_valid"]] == 0)
        assert np.all(fields["energy_ratio"][~fields["energy_valid"]] != fields["energy_ratio"][~fields["energy_valid"]])
        print(f"{sample} {name}: Gen={len(photons['energy'])}, with Reco={np.count_nonzero(~unmatched)}, "
              f"raw objects={fields['raw_count'].sum()}, canonical objects={fields['energy_count'].sum()}", flush=True)
    return dict(events=events, photons=photons, channels=channels, composition=composition,
                link_validation={name: dict(values) for name, values in validation.items()}, diagnostics=dict(diagnostics))
