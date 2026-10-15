"""Track-associated saved RICH fields, without interpreting PID quality bits."""
from collections import Counter

import numpy as np
import uproot

from common.dataset import SOURCE_TREES
from common.matching import angle_matrix, match_tracks
from ..track import data as track_data

ANGLE_CUT, SCAN_CUTS = track_data.ANGLE_CUT, track_data.SCAN_CUTS
read_events, truth_track_pairs = track_data.read_events, track_data.truth_track_pairs

SOURCES = {"nanoaod_raw_sdst": "HaidRaw", "nanoaod": "Haid_Rich"}
FIELDS = {
    "nanoaod_raw_sdst": {
        "quality": "HaidRaw_richQuality",
        **{f"{radiator}_{field}": f"HaidRaw_{branch}{suffix}"
           for radiator, suffix in (("gas", "Gas"), ("liquid", "Liq"))
           for field, branch in (("angle", "theta"), ("sigma", "sigma"),
                                 ("photons", "nph"), ("expected", "nep"), ("flag", "flag"))},
    },
    "nanoaod": {
        "quality": "Haid_richQuality",
        **{f"{radiator}_{field}": f"Rich_{branch}{suffix}"
           for radiator, suffix in (("gas", "g"), ("liquid", "l"))
           for field, branch in (("angle", "the"), ("sigma", "sig"),
                                 ("photons", "nph"), ("expected", "nep"), ("flag", "flag"))},
    },
}
LINK_LABELS = ("Agree", "Forward only", "Reverse only", "Conflict", "No origin")


def saved_angle(angle):
    """A positive normal stored angle; this alone is not a PID decision."""
    angle = np.asarray(angle)
    return np.isfinite(angle) & (angle >= np.finfo(np.float32).tiny)


def corrected_schema(files, source):
    """Do not use the old ABI or VECP-as-PA indices for RICH coverage."""
    for path in files:
        with uproot.open(path) as root:
            tree = root[SOURCE_TREES[source]]
            if source == "nanoaod_raw_sdst":
                # The paired native schema records the corrected producer even
                # when a rewritten ROOT file has no original branch descriptions.
                if "Part_paIdx" not in root[SOURCE_TREES["nanoaod"]]:
                    return False
            elif "Part_paIdx" not in tree:
                return False  # Added with the corrected RICH producer.
    return True


def track_fields(event, source, supported=True):
    """Join raw HaidRaw and TracRaw through their zero-based PA indices.

    In nanoaod, Haid/Rich rows use Part order, reached through Trac_partIdx. Duplicate raw
    PA IDs remain unresolved. Old raw VECP indices cannot be used as PA indices.
    """
    raw, tracks = event["raw"], event["reco_idx"]
    if not supported:
        rows = np.full(len(tracks), -1, dtype=int)
    elif source == "nanoaod":
        rows = event["part_idx"]
    else:
        ids = np.asarray(raw["HaidRaw_paIdx"], dtype=int)
        unique, counts = np.unique(ids, return_counts=True)
        usable = set(unique[(counts == 1) & (unique >= 0)])
        lookup = {part: row for row, part in enumerate(ids) if part in usable}
        rows = np.asarray([lookup.get(raw["TracRaw_paIdx"][track], -1) for track in tracks], dtype=int)
    result = {}
    for name, branch in FIELDS[source].items():
        values = np.asarray(raw[branch])
        selected = np.full(len(tracks), np.nan)
        valid = (rows >= 0) & (rows < len(values))
        selected[valid] = values[rows[valid]]
        result[name] = selected
    result["linked"] = (rows >= 0) & (rows < len(raw[FIELDS[source]["quality"]]))
    return result


def read_sample(files, source):
    """Read once and retain the small arrays needed for the four plot stages."""
    supported = corrected_schema(files, source)
    gen = {name: [] for name in ("p", "cos_theta", "pdg")}
    reco = {name: [] for name in ("p", "cos_theta")}
    stored = {name: [] for name in FIELDS[source]}
    methods = ("angular", "truth") if source == "nanoaod" else ("angular",)
    results = {method: {name: [] for name in ("matched", "gas", "liquid", "either", "angle",
                                               "reco_count", "saved_count", "track_p", "track_cos_theta",
                                               "track_gas", "track_liquid", "track_either")}
               for method in methods}
    scan = {name: [] for name in ("track", "gas", "liquid", "either")}
    multiplicity = {name: [] for name in ("gen", "reco", "gas", "liquid", "either")}
    links, counts = Counter(dict.fromkeys(LINK_LABELS, 0)), Counter()
    branches = list(FIELDS[source].values())
    if source == "nanoaod_raw_sdst":
        branches += ["HaidRaw_paIdx", "nHaidRaw"]
    for event in read_events(files, source, branches):
        fields = track_fields(event, source, supported)
        raw = event["raw"]
        bank_ids = raw["HaidRaw_paIdx"] if source == "nanoaod_raw_sdst" else range(len(raw["Part_charge"]))
        counts["bank_length_mismatch_events"] += int(any(len(raw[branch]) != len(bank_ids)
                                                        for branch in FIELDS[source].values()))
        usable_ids = [index for index in bank_ids if index >= 0]
        counts["duplicate_bank_id_events"] += int(len(set(usable_ids)) != len(usable_ids))
        if source == "nanoaod_raw_sdst":
            counts["reported_count_mismatch_events"] += int(raw["nHaidRaw"] != len(bank_ids))
        gen_xyz, reco_xyz = event["gen_xyz"], event["reco_xyz"]
        gen_p, reco_p = np.linalg.norm(gen_xyz, axis=1), np.linalg.norm(reco_xyz, axis=1)
        for target, xyz, momentum in ((gen, gen_xyz, gen_p), (reco, reco_xyz, reco_p)):
            target["p"].append(momentum)
            target["cos_theta"].append(xyz[:, 2] / momentum)
        gen["pdg"].append(event["gen_pdg"])
        for name, branch in FIELDS[source].items():
            stored[name].append(np.asarray(raw[branch]))
        coverage = {name: saved_angle(fields[f"{name}_angle"]) for name in ("gas", "liquid")}
        coverage["either"] = coverage["gas"] | coverage["liquid"]
        multiplicity["gen"].append(len(gen_p))
        multiplicity["reco"].append(len(reco_p))
        for name in coverage:
            multiplicity[name].append(int(coverage[name].sum()))
        counts.update(events=1, gen=len(gen_p), tracks=len(reco_p))
        if supported:
            counts["unmapped_rich_rows"] += int((~fields["linked"]).sum())
        else:
            counts["association_not_evaluated"] += len(reco_p)
        angles = angle_matrix(gen_xyz, reco_xyz)
        gi, ri = match_tracks(gen_xyz, reco_xyz, event["gen_charge"], event["reco_charge"], SCAN_CUTS[-1])
        for name in scan:
            value = np.full(len(gen_p), np.inf)
            selected = np.ones(len(ri), dtype=bool) if name == "track" else coverage[name][ri]
            value[gi[selected]] = angles[gi[selected], ri[selected]]
            scan[name].append(value)
        selected = angles[gi, ri] < ANGLE_CUT
        pairs = {"angular": (gi[selected], ri[selected])}
        if source == "nanoaod":
            truth = event["truth"]
            pairs["truth"] = truth_track_pairs(event)
            for part in event["part_idx"]:
                forward, reverse, known = truth["forward"][part], truth["reverse"][part], truth["known"][part]
                label = ("Conflict" if len(known) > 1 else "No origin" if not known else
                         "Agree" if forward >= 0 and reverse else "Forward only" if forward >= 0 else "Reverse only")
                links[label] += 1
        for method, (gi, ri) in pairs.items():
            result = results[method]
            for name in ("matched", "gas", "liquid", "either"):
                value = np.zeros(len(gen_p), dtype=bool)
                selected = np.ones(len(ri), dtype=bool) if name == "matched" else coverage[name][ri]
                value[gi[selected]] = True
                result[name].append(value)
            result["reco_count"].append(np.bincount(gi, minlength=len(gen_p)))
            result["saved_count"].append(np.bincount(gi[coverage["either"][ri]], minlength=len(gen_p)))
            result["angle"].append(angles[gi, ri])
            result["track_p"].append(reco_p[ri])
            result["track_cos_theta"].append(reco_xyz[ri, 2] / reco_p[ri])
            for name in coverage:
                result[f"track_{name}"].append(coverage[name][ri])
    return dict(gen={key: np.concatenate(value) for key, value in gen.items()},
                reco={key: np.concatenate(value) for key, value in reco.items()},
                stored={key: np.concatenate(value) for key, value in stored.items()},
                methods={method: {key: np.concatenate(value) for key, value in result.items()}
                         for method, result in results.items()},
                scan={key: np.concatenate(value) for key, value in scan.items()},
                multiplicity={key: np.asarray(value) for key, value in multiplicity.items()},
                links=dict(links), counts=dict(counts), association_supported=supported)
