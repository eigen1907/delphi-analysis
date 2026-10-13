"""Stable Gen photons and truth-linked Reco photon candidates."""
from collections import Counter
from pathlib import Path

import numpy as np
import uproot

from common.matching import descendants, link_status, resolve_part_origins

SAMPLES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")
TARGET_PDG = {"Zee": 11, "Zmumu": 13, "Ztautau": 15, "ZKK": 321, "Zpipi": 211}
SIM_GAMMA_CODE = 21  # Part/SimPart use DELPHI mass codes; GenPart uses PDG IDs.
GEN_P4 = tuple(f"GenPart_vector.fCoordinates.f{axis}" for axis in "XYZT")
RECO_P4 = tuple(f"Photon_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
PART_P4 = tuple(f"Part_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
CONV_P4 = tuple(f"PhotonConv_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT")
GEN_SELECTIONS = ("gen_gamma", "gen_isr", "gen_no_isr", "gen_isr_non_beam")
ANALYSIS_SELECTIONS = ("gen_gamma", "gen_isr", "gen_no_isr")
GEN_LABELS = (r"$\gamma_{\mathrm{all}}$", r"$\gamma_{\mathrm{ISR}}$", r"$\gamma_{\mathrm{no\ ISR}}$")
GEN_CUTS = (("cut_no", 0, False), ("cut_energy_0p1", 0.1, False),
            ("cut_fiducial", 0, True), ("cut_fiducial_energy_0p1", 0.1, True))
EM_REGIONS = {"HPC": ((44, 88), (92, 136)),
              "FEMC": ((11, 36), (144, 169)),
              "STIC": ((2, 10), (170, 178))}
RECO_LABELS = {"Photon": r"$\gamma$", "Photon_PhotonConv": r"$\gamma+\gamma_{\mathrm{conv}}$",
               "Part_pdgId_21": r"$\gamma_{\mathrm{Part}}$"}
LINK_STATUSES = ("agreement", "forward_only", "reverse_only", "conflict", "no_origin")
BRANCHES = (
    "GenPart_pdgId", "GenPart_status", "GenPart_parentIdx", *GEN_P4,
    "Photon_partIdx", *RECO_P4, "Part_simIdx", "Part_charge", "Part_pdgId", *PART_P4,
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
    """Matching studies exclude the beam-collinear generator component."""
    return dict(gen_gamma=category != "beam_isr", gen_isr=category == "nonbeam_isr",
                gen_no_isr=np.isin(category, ("fsr", "decayed")))


def em_regions(theta):
    """Conservative Gen acceptance, with theta in degrees and open bounds."""
    theta = np.asarray(theta)
    return {name: np.logical_or.reduce([(theta > low) & (theta < high)
                                       for low, high in intervals])
            for name, intervals in EM_REGIONS.items()}


def em_acceptance(cos_theta):
    """Gen acceptance from the union of the same theta regions used in plots."""
    theta = np.rad2deg(np.arccos(np.clip(cos_theta, -1, 1)))
    return np.logical_or.reduce(list(em_regions(theta).values()))


def matching_masks(energy, cos_theta, category):
    """All four cuts apply to Gen only, after the Beam ISR exclusion."""
    populations = analysis_selections(category)
    accepted = em_acceptance(cos_theta)
    return {name: {population: selected & (energy >= minimum) & (accepted if fiducial else True)
                   for population, selected in populations.items()}
            for name, minimum, fiducial in GEN_CUTS}


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


def gen_event(raw, sample, p4_branches=GEN_P4):
    """The same stable-photon populations in both NanoAOD files."""
    pdgs, statuses, parents = (np.asarray(raw[name], dtype=int) for name in
                               ("GenPart_pdgId", "GenPart_status", "GenPart_parentIdx"))
    daughters = np.flatnonzero((np.abs(pdgs) == TARGET_PDG[sample]) &
                              (statuses == (21 if sample == "Ztautau" else 1)))
    assert len(daughters) == 2 and parents[daughters[0]] == parents[daughters[1]]
    hard_parent = parents[daughters[0]]
    all_gamma_idx = np.flatnonzero((pdgs == 22) & (statuses == 1))
    all_gamma_p4 = np.column_stack([raw[name] for name in p4_branches])[all_gamma_idx]
    all_gamma_category = []
    for gen, p4 in zip(all_gamma_idx, all_gamma_p4, strict=True):
        origin = photon_origin(gen, pdgs, parents, hard_parent, TARGET_PDG[sample])
        is_beam = (origin == "isr" and p4[0] == p4[1] == 0 and 0 <= parents[gen] < len(pdgs)
                   and abs(pdgs[parents[gen]]) == 11 and parents[parents[gen]] < 0)
        all_gamma_category.append("beam_isr" if is_beam else "nonbeam_isr" if origin == "isr"
                                  else "fsr" if origin == "fsr" else "decayed")
    all_gamma_category = np.asarray(all_gamma_category)
    return dict(all_gamma_idx=all_gamma_idx, all_gamma_p4=all_gamma_p4,
                all_gamma_category=all_gamma_category, selections=gen_selections(all_gamma_category),
                analysis_selections=analysis_selections(all_gamma_category),
                cuts=matching_masks(all_gamma_p4[:, 3], unit_momenta(all_gamma_p4)[:, 2], all_gamma_category))


def analyze_event(raw, sample):
    """Resolve all Gen origins before selecting the stable-photon populations."""
    event = gen_event(raw, sample)
    all_gamma_idx, all_gamma_p4 = event['all_gamma_idx'], event['all_gamma_p4']
    stable_gamma = set(all_gamma_idx)
    n_part = len(raw["Part_simIdx"])
    truth = resolve_part_origins(raw)
    forward, reverse, known = truth["forward"], truth["reverse"], truth["known"]
    part_origins = truth["origin"]
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

    # This selection starts from all Parts, independent of the Photon collection.
    # Part_pdgId stores the DELPHI mass code, not the standard PDG ID.
    # Neutral Parts without an identification can also receive the default code 21.
    pid_parts = np.flatnonzero(np.asarray(raw["Part_pdgId"]) == SIM_GAMMA_CODE)
    pid_p4 = np.column_stack([raw[name] for name in PART_P4])[pid_parts]
    pid_origins = part_origins[pid_parts]
    all_pid_parts = set(pid_parts)
    pid_geometry_valid = np.asarray([
        not (part_descendants.get(part, set()) & all_pid_parts) for part in pid_parts
    ], dtype=bool)
    pid_valid = np.asarray([
        part not in bad_parts and len(known[part]) <= 1 for part in pid_parts
    ], dtype=bool) & pid_geometry_valid
    pid_stats = response(all_gamma_idx, all_gamma_p4, pid_origins, pid_p4,
                         pid_valid, pid_origins, set())
    reco = dict(
        Photon=dict(p4=photon_p4, gen_idx=photon_origins, energy_valid=photon_valid,
                   geometric_energy_valid=photon_geometry_valid, stats=photon_stats),
        Photon_PhotonConv=dict(p4=gamma_p4, gen_idx=gamma_origins, energy_valid=gamma_valid,
                                   geometric_energy_valid=gamma_geometry_valid, stats=combined_stats),
        Part_pdgId_21=dict(p4=pid_p4, gen_idx=pid_origins, energy_valid=pid_valid,
                       geometric_energy_valid=pid_geometry_valid, stats=pid_stats),
    )
    links = dict(
        Photon=([{int(forward[p])} if forward[p] >= 0 else set() for p in photon_parts],
               [reverse[p] for p in photon_parts]),
        Photon_PhotonConv=(gamma_forward, gamma_reverse),
        Part_pdgId_21=([{int(forward[p])} if forward[p] >= 0 else set() for p in pid_parts],
                   [reverse[p] for p in pid_parts]),
    )
    cuts = event['cuts']
    matching_validation, link_counts = {}, {}
    selected_origins = {cut: {name: set(all_gamma_idx[mask]) for name, mask in masks.items()}
                        for cut, masks in cuts.items()}
    for mode, (forward_links, reverse_links) in links.items():
        counts = {cut: {gen: Counter(dict.fromkeys(LINK_STATUSES, 0)) for gen in ANALYSIS_SELECTIONS}
                  for cut, _, _ in GEN_CUTS}
        link_counts[mode] = Counter(dict.fromkeys(LINK_STATUSES, 0))
        for first, second in zip(forward_links, reverse_links, strict=True):
            status = link_status(first, second)
            link_counts[mode][status] += 1
            for cut, origins in selected_origins.items():
                for gen_name, selected in origins.items():
                    if (first | second) & selected:
                        counts[cut][gen_name][status] += 1
        matching_validation[mode] = counts
    event.update(reco=reco, matching_validation=matching_validation, link_counts=link_counts)
    return event


def read_sample(input_root: Path, sample: str, angle_accumulator, source='nanoaod'):
    """Stream files and retain compact observables; all cuts share the same Reco pool."""
    from hashlib import sha1
    from . import raw_data

    native = source == 'nanoaod'
    labels = RECO_LABELS if native else raw_data.RECO_LABELS
    branches = BRANCHES if native else raw_data.BRANCHES
    gen_p4 = GEN_P4 if native else raw_data.GEN_P4
    analyze = analyze_event if native else raw_data.analyze_event
    gamma = {name: [] for name in ("energy", "cos_theta", "category")}
    results = {mode: dict(
        reco={field: [] for field in ("reco_count", "energy_ratio", "energy_valid")},
        reco_distributions={field: [] for field in ("energy", "cos_theta")},
        reco_event_counts=[], link_counts=Counter(dict.fromkeys(LINK_STATUSES, 0)),
        matching_validation={cut: {gen: Counter(dict.fromkeys(LINK_STATUSES, 0))
                                  for gen in ANALYSIS_SELECTIONS} for cut, _, _ in GEN_CUTS},
        branch_counts=Counter(),
    ) for mode in labels}
    gen_event_counts = {name: [] for name in
                        ("all_gamma", "all_isr", "beam_isr", "nonbeam_isr", "fsr", "decayed")}
    events, metadata_excluded = 0, 0
    event_keys, gen_signatures = [], []
    directory = input_root / f"20260828_100kTest_{sample}_photosFSR" / "final_root"
    paths = sorted(directory.glob(f"job_*/{source}.root"))
    if not paths:
        raise FileNotFoundError(directory)
    for path in paths:
        with uproot.open(path) as root_file:
            tree = root_file["Events"]
            event_number = 'Event_evtNumber' if native else 'Event_eventNumber'
            names = (*branches, 'Event_runNumber', event_number)
            columns = {name: tree[name].array(library="ak").to_list() for name in names}
            for entry in range(tree.num_entries):
                raw = {name: values[entry] for name, values in columns.items()}
                # Raw SDST has one empty initialization record per job, not a physics event.
                if not raw['GenPart_pdgId']:
                    metadata_excluded += 1
                    continue
                events += 1
                event_keys.append((path.parent.name, raw['Event_runNumber'], raw[event_number]))
                signature = sha1()
                for name in ('GenPart_pdgId', 'GenPart_status', 'GenPart_parentIdx', *gen_p4):
                    signature.update(np.asarray(raw[name], dtype=np.float64).tobytes())
                gen_signatures.append(signature.digest())
                event = analyze(raw, sample)
                p4, categories = event["all_gamma_p4"], event["all_gamma_category"]
                counts = dict(all_gamma=len(categories), all_isr=int(event["selections"]["gen_isr"].sum()))
                counts.update({name: int(np.count_nonzero(categories == name))
                               for name in ("beam_isr", "nonbeam_isr", "fsr", "decayed")})
                for name, count in counts.items():
                    gen_event_counts[name].append(count)
                gamma["energy"].extend(p4[:, 3])
                gamma["cos_theta"].extend(unit_momenta(p4)[:, 2])
                gamma["category"].extend(categories)
                for mode, objects in event["reco"].items():
                    values = results[mode]
                    if native:
                        for field, entries in values["reco"].items():
                            entries.extend(objects["stats"][field])
                    reco_p4 = objects["p4"]
                    momentum = np.linalg.norm(reco_p4[:, :3], axis=1)
                    cos_theta = np.divide(reco_p4[:, 2], momentum, out=np.full(len(reco_p4), np.nan), where=momentum > 0)
                    values["reco_distributions"]["energy"].extend(reco_p4[:, 3])
                    values["reco_distributions"]["cos_theta"].extend(cos_theta)
                    values["reco_event_counts"].append(len(reco_p4))
                    if native:
                        values["link_counts"].update(event["link_counts"][mode])
                        for cut, _, _ in GEN_CUTS:
                            for population in ANALYSIS_SELECTIONS:
                                values["matching_validation"][cut][population].update(
                                    event["matching_validation"][mode][cut][population])
                    else:
                        values['branch_counts'].update(objects['branch_counts'])
                angle_accumulator(event)
        print(f"{sample} {source} {path.parent.name}: {events} events, empty Gen records excluded={metadata_excluded}", flush=True)
    gamma = {name: np.asarray(entries) for name, entries in gamma.items()}
    gen_event_counts = {name: np.asarray(entries) for name, entries in gen_event_counts.items()}
    shared = dict(events=events, source=source, event_keys=event_keys, gen_signatures=gen_signatures,
                  metadata_excluded=metadata_excluded, gamma=gamma, selections=gen_selections(gamma["category"]),
                  analysis_selections=analysis_selections(gamma["category"]),
                  cuts=matching_masks(gamma["energy"], gamma["cos_theta"], gamma["category"]),
                  gen_event_counts=gen_event_counts,
                  gen_event_index=np.repeat(np.arange(events), gen_event_counts["all_gamma"]))
    for mode, values in results.items():
        for field in ("reco", "reco_distributions"):
            values[field] = {name: np.asarray(entries) for name, entries in values[field].items()}
        values["reco_event_counts"] = np.asarray(values["reco_event_counts"])
        values.update(shared, collection=mode, reco_label=labels[mode])
        print(f"{sample} {source}/{mode}: Gen={len(gamma['energy'])}, Reco={len(values['reco_distributions']['energy'])}", flush=True)
    return results
