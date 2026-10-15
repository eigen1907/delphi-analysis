"""Read stable charged Gen particles and the two saved track collections."""
from collections import Counter
from hashlib import sha256

import awkward as ak
import numpy as np
import uproot

from common.dataset import SOURCE_TREES
from common.matching import angle_matrix, match_tracks, resolve_part_origins

SOURCES = {"nanoaod_raw_sdst": "TracRaw", "nanoaod": "Trac"}
CHARGED_PDGS = (11, 13, 211, 321, 2212)
ANGLE_CUT = 0.05  # radians; same-charge, greedy one-to-one angular matching.
SCAN_CUTS = np.linspace(0, 0.20, 201)
TRUTH_BRANCHES = ("Part_simIdx", "SimPart_genIdx", "SimPart_partIdx",
                  "SimPart_originVtxIdx", "SimVtx_incomingIdx")
LINK_LABELS = ("Agree", "Forward only", "Reverse only", "Conflict", "No origin")


def truth_track_pairs(event):
    """All direct/descendant truth associations with the same charge as Gen.

    Opposite-charge shower tracks do not reconstruct the generated charged track.
    No angle or PID requirement is added to the unique first-Gen-anchor link.
    """
    origin = event["truth"]["origin"][event["part_idx"]]
    eligible = ((event["gen_idx"][:, None] == origin[None, :]) &
                (event["gen_charge"][:, None] == event["reco_charge"][None, :]))
    return np.argwhere(eligible).T


def read_events(files, source, extra_branches=()):
    """Yield events; extra branches and original row indices let RICH reuse tracks.

    Raw track momentum is B/|invR|. Normal Trac momentum and charge are read
    through Trac_partIdx -> Part. No track-quality, PID or fiducial cut is added.
    Both Gen and Reco require finite momenta and nonzero transverse momentum;
    Gen also requires finite positive energy.
    """
    gen_name = "GenPart_fourMomentum" if source == "nanoaod_raw_sdst" else "GenPart_vector"
    gen_p4 = tuple(f"{gen_name}.fCoordinates.f{axis}" for axis in "XYZT")
    event_number = "Event_eventNumber" if source == "nanoaod_raw_sdst" else "Event_evtNumber"
    branches = ["job_id", "GenPart_status", "GenPart_pdgId", "GenPart_parentIdx", "Event_runNumber", event_number, *gen_p4]
    if source == "nanoaod_raw_sdst":
        branches += ["Event_bFieldGevCm", "TracRaw_theta", "TracRaw_phi",
                     "TracRaw_invR", "TracRaw_charge", "TracRaw_paIdx"]
    else:
        branches += ["Trac_partIdx", "Part_charge", *TRUTH_BRANCHES]
        branches += [f"Part_fourMomentum.fCoordinates.f{axis}" for axis in "XYZT"]
    for path in files:
        with uproot.open(path) as root:
            tree = root[SOURCE_TREES[source]]
            excluded = 0
            for start in range(0, tree.num_entries, 5000):
                # Read split four-vector leaves directly; tree.iterate omits some
                # dotted leaves in these merged ROOT files.
                columns = {name: ak.to_list(tree[name].array(entry_start=start, entry_stop=start + 5000, library="ak"))
                           for name in dict.fromkeys([*branches, *extra_branches])}
                for entry in range(len(next(iter(columns.values())))):
                    raw = {name: values[entry] for name, values in columns.items()}
                    pdg = np.asarray(raw["GenPart_pdgId"], dtype=int)
                    if len(pdg) == 0:
                        # Raw SDST saves an empty setup record before the MC events.
                        # This removes metadata, not events with zero selected tracks.
                        excluded += 1
                        continue
                    p4 = np.column_stack([raw[name] for name in gen_p4])
                    signature = sha256()
                    for name in ("GenPart_pdgId", "GenPart_status", "GenPart_parentIdx"):
                        signature.update(np.asarray(raw[name], dtype=np.int64).tobytes())
                    signature.update(np.asarray(p4, dtype=np.float64).tobytes())
                    pt = np.linalg.norm(p4[:, :2], axis=1)
                    selected = ((np.asarray(raw["GenPart_status"]) == 1) & np.isin(abs(pdg), CHARGED_PDGS)
                                & np.isfinite(p4).all(axis=1) & (p4[:, 3] > 0) & (pt > 0))
                    gen_idx = np.flatnonzero(selected)
                    gen_charge = np.sign(pdg[selected]) * np.where(np.isin(abs(pdg[selected]), (11, 13)), -1, 1)
                    truth = None
                    if source == "nanoaod_raw_sdst":
                        theta, phi, inv_r = (np.asarray(raw[f"TracRaw_{name}"]) for name in ("theta", "phi", "invR"))
                        charge = np.asarray(raw["TracRaw_charge"], dtype=int)
                        raw_count = len(theta)
                        track_pt = np.divide(raw["Event_bFieldGevCm"], abs(inv_r),
                                             out=np.full(raw_count, np.nan), where=inv_r != 0)
                        momentum = np.divide(track_pt, np.sin(theta), out=np.full(raw_count, np.nan),
                                             where=np.sin(theta) > 0)
                        xyz = np.column_stack((track_pt * np.cos(phi), track_pt * np.sin(phi),
                                               momentum * np.cos(theta)))
                        part_idx = np.full(raw_count, -1, dtype=int)
                    else:
                        part_idx = np.asarray(raw["Trac_partIdx"], dtype=int)
                        raw_count = len(part_idx)
                        valid_part = (part_idx >= 0) & (part_idx < len(raw["Part_charge"]))
                        xyz = np.full((raw_count, 3), np.nan)
                        charge = np.zeros(raw_count, dtype=int)
                        part_p4 = np.column_stack([raw[f"Part_fourMomentum.fCoordinates.f{axis}"] for axis in "XYZT"])
                        xyz[valid_part] = part_p4[part_idx[valid_part], :3]
                        charge[valid_part] = np.asarray(raw["Part_charge"])[part_idx[valid_part]]
                        truth = resolve_part_origins(raw)
                    valid = np.isfinite(xyz).all(axis=1) & (np.linalg.norm(xyz[:, :2], axis=1) > 0) & (charge != 0)
                    reco_idx = np.flatnonzero(valid)
                    yield dict(raw=raw, gen_idx=gen_idx, gen_xyz=p4[selected, :3],
                               gen_energy=p4[selected, 3], gen_pdg=pdg[selected], gen_charge=gen_charge,
                               reco_idx=reco_idx, reco_xyz=xyz[valid], reco_charge=charge[valid],
                               gen_p=np.linalg.norm(p4[selected, :3], axis=1), gen_pt=pt[selected],
                               reco_p=np.linalg.norm(xyz[valid], axis=1),
                               reco_pt=np.linalg.norm(xyz[valid, :2], axis=1),
                               part_idx=part_idx[valid], truth=truth, raw_reco_count=raw_count,
                               event_key=(int(raw["job_id"]), int(raw["Event_runNumber"]), int(raw[event_number])),
                               gen_signature=signature.digest())
            if excluded:
                print(f"{source} {path.parent.name}: metadata_excluded={excluded} (empty Gen record)", flush=True)


def read_sample(files, source):
    """Keep compact distributions and per-Gen matched track quantities."""
    gen = {name: [] for name in ("p", "pt", "cos_theta", "phi")}
    reco = {name: [] for name in gen}
    gen["pdg"] = []
    methods = ("angular", "truth") if source == "nanoaod" else ("angular",)
    results = {method: {name: [] for name in ("matched", "p_response", "pt_response", "angle", "theta_residual", "phi_residual", "reco_count", "momentum_ratio")}
               for method in methods}
    gen_counts, reco_counts, scan_angles = [], [], []
    links = Counter(dict.fromkeys(LINK_LABELS, 0))
    counts = Counter()
    records, duplicates = {}, set()
    for path in files:
        with uproot.open(path) as root:
            counts["metadata_excluded"] += int(np.count_nonzero(root[SOURCE_TREES[source]]["nGenPart"].array(library="np") == 0))
    for event in read_events(files, source):
        gen_xyz, reco_xyz = event["gen_xyz"], event["reco_xyz"]
        gen_p, reco_p = np.linalg.norm(gen_xyz, axis=1), np.linalg.norm(reco_xyz, axis=1)
        gen_pt, reco_pt = np.linalg.norm(gen_xyz[:, :2], axis=1), np.linalg.norm(reco_xyz[:, :2], axis=1)
        for target, xyz, momentum, transverse in ((gen, gen_xyz, gen_p, gen_pt), (reco, reco_xyz, reco_p, reco_pt)):
            target["p"].append(momentum)
            target["pt"].append(transverse)
            target["cos_theta"].append(xyz[:, 2] / momentum)
            target["phi"].append(np.arctan2(xyz[:, 1], xyz[:, 0]))
        gen["pdg"].append(event["gen_pdg"])
        if event["event_key"] in records:
            duplicates.add(event["event_key"])
        records[event["event_key"]] = event["gen_signature"]
        gen_counts.append(len(gen_p))
        reco_counts.append(len(reco_p))
        counts.update(events=1, stable_charged=len(gen_p), raw_tracks=event["raw_reco_count"], valid_tracks=len(reco_p))
        angles = angle_matrix(gen_xyz, reco_xyz)
        gen_match, reco_match = match_tracks(gen_xyz, reco_xyz, event["gen_charge"], event["reco_charge"], SCAN_CUTS[-1])
        scan = np.full(len(gen_p), np.inf)
        scan[gen_match] = angles[gen_match, reco_match]
        scan_angles.append(scan)
        pairs = {"angular": (gen_match[scan[gen_match] < ANGLE_CUT], reco_match[scan[gen_match] < ANGLE_CUT])}
        multiplicity = {"angular": np.bincount(pairs["angular"][0], minlength=len(gen_p))}
        if source == "nanoaod":
            truth = event["truth"]
            gen_rows, reco_rows = truth_track_pairs(event)
            multiplicity["truth"] = np.bincount(gen_rows, minlength=len(gen_p))
            truth_gen = np.unique(gen_rows)
            # A Gen track counts once; the leading-p associated track defines response.
            truth_reco = [reco_rows[gen_rows == row][np.argmax(reco_p[reco_rows[gen_rows == row]])]
                          for row in truth_gen]
            pairs["truth"] = (truth_gen, np.asarray(truth_reco, dtype=int))
            for part in event["part_idx"]:
                forward, reverse, known = truth["forward"][part], truth["reverse"][part], truth["known"][part]
                label = ("Conflict" if len(known) > 1 else "No origin" if not known else
                         "Agree" if forward >= 0 and reverse else "Forward only" if forward >= 0 else "Reverse only")
                links[label] += 1
        for method, (gen_match, reco_match) in pairs.items():
            matched = np.zeros(len(gen_p), dtype=bool)
            matched[gen_match] = True
            p_response, pt_response, angle = (np.full(len(gen_p), np.nan) for _ in range(3))
            p_response[gen_match] = (reco_p[reco_match] - gen_p[gen_match]) / gen_p[gen_match]
            pt_response[gen_match] = (reco_pt[reco_match] - gen_pt[gen_match]) / gen_pt[gen_match]
            angle[gen_match] = angles[gen_match, reco_match]
            theta_residual, phi_residual = (np.full(len(gen_p), np.nan) for _ in range(2))
            theta_residual[gen_match] = (np.arctan2(reco_pt[reco_match], reco_xyz[reco_match, 2])
                                       - np.arctan2(gen_pt[gen_match], gen_xyz[gen_match, 2]))
            delta_phi = (np.arctan2(reco_xyz[reco_match, 1], reco_xyz[reco_match, 0])
                         - np.arctan2(gen_xyz[gen_match, 1], gen_xyz[gen_match, 0]))
            phi_residual[gen_match] = (delta_phi + np.pi) % (2 * np.pi) - np.pi
            momentum_ratio = np.zeros(len(gen_p))
            momentum_ratio[gen_match] = reco_p[reco_match] / gen_p[gen_match]
            for key, value in zip(results[method], (matched, p_response, pt_response, angle, theta_residual, phi_residual,
                                                   multiplicity[method], momentum_ratio), strict=True):
                results[method][key].append(value)
    return dict(gen={key: np.concatenate(value) for key, value in gen.items()},
                reco={key: np.concatenate(value) for key, value in reco.items()},
                methods={name: {key: np.concatenate(value) for key, value in fields.items()}
                         for name, fields in results.items()}, scan_angle=np.concatenate(scan_angles),
                gen_counts=np.asarray(gen_counts), reco_counts=np.asarray(reco_counts),
                counts=dict(counts), links=dict(links), event_records=records, duplicate_keys=duplicates)
