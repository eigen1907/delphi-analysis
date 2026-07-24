from __future__ import annotations

import awkward as ak
import numpy as np


SOURCE_FILE_NAME = "nanoaod_raw_sdst.root"
TREE_NAME = "Events"
N_OBJECTS = 4

EVENT_FEATURES = (
    "n_track",
    "n_em_shower",
    "n_had_shower",
    "n_muon",
    "n_electron",
    "n_hadron_pid",
    "n_vertex",
    "em_energy_sum",
    "had_energy_sum",
    "calorimeter_energy_sum",
    "track_p_sum",
    "track_momentum_imbalance",
)
TRACK_FEATURES = (
    "p",
    "abs_cos_theta",
    "charge",
    "impact_rphi",
    "impact_z",
    "is_muon",
    "is_electron",
    "em_energy_over_p",
    "had_energy_over_p",
    "kaon_combined",
    "kaon_dedx",
    "kaon_rich",
    "pion_rich",
    "proton_combined",
    "proton_dedx",
    "proton_rich",
    "rich_quality",
    "tpc_dedx80",
    "tpc_dedx80_sigma",
    "tpc_dedx65",
    "tpc_dedx65_sigma",
    "tpc_dedx80_integrated",
)
SHOWER_FEATURES = (
    *(f"em_shower{rank}_energy" for rank in range(1, N_OBJECTS + 1)),
    *(
        f"had_shower{rank}_{feature}"
        for rank in range(1, N_OBJECTS + 1)
        for feature in ("energy", "n_hits", "n_hit_rows")
    ),
)
PAIR_FEATURES = (
    "track12_acolinearity",
    "track12_momentum_ratio",
    "track12_charge_product",
    "track12_tpc_dedx80_min",
    "track12_tpc_dedx80_max",
    "track12_tpc_dedx80_abs_diff",
    "track12_kaon_dedx_min",
    "track12_kaon_combined_min",
    "track12_both_have_tpc",
    "track12_both_have_had_pid",
)
FEATURE_NAMES = (
    *EVENT_FEATURES,
    *(
        f"track{rank}_{feature}"
        for rank in range(1, N_OBJECTS + 1)
        for feature in TRACK_FEATURES
    ),
    *SHOWER_FEATURES,
    *PAIR_FEATURES,
)

INPUT_BRANCHES = (
    "Event_bFieldGevCm",
    "nTracRaw",
    "TracRaw_paIdx",
    "TracRaw_impactRPhi",
    "TracRaw_impactZ",
    "TracRaw_theta",
    "TracRaw_phi",
    "TracRaw_invR",
    "TracRaw_charge",
    "nEmShower",
    "EmShower_paIdx",
    "EmShower_energy",
    "nHadShower",
    "HadShower_paIdx",
    "HadShower_energy",
    "HadShower_nHits",
    "HadShower_nHitRows",
    "nMuidRaw",
    "MuidRaw_paIdx",
    "nElidRaw",
    "ElidRaw_paIdx",
    "nHaidRaw",
    "HaidRaw_paIdx",
    "HaidRaw_kaonCombined",
    "HaidRaw_kaonDedx",
    "HaidRaw_kaonRich",
    "HaidRaw_pionRich",
    "HaidRaw_protonCombined",
    "HaidRaw_protonDedx",
    "HaidRaw_protonRich",
    "HaidRaw_richQuality",
    "MtpcRaw_tracRawIdx",
    "MtpcRaw_dEdx80Max",
    "MtpcRaw_dEdx80Sigma",
    "MtpcRaw_dEdx65Max",
    "MtpcRaw_dEdx65Sigma",
    "MtpcRaw_dEdx80Integrated",
    "nVtx",
)


def _read(tree, branch: str):
    if branch not in tree:
        raise KeyError(f"Missing BDT input branch: {branch}")
    return tree[branch].array(library="ak")


def _to_float32(array) -> np.ndarray:
    return np.asarray(ak.to_numpy(array), dtype=np.float32)


def _top_k(array, order=None) -> np.ndarray:
    if order is not None:
        array = array[order]
    padded = ak.pad_none(array, N_OBJECTS, axis=1, clip=True)
    return _to_float32(ak.fill_none(padded, np.nan))


def _top_k_descending(array) -> np.ndarray:
    order = ak.argsort(array, axis=1, ascending=False)
    return _top_k(array, order)


def _event_sum(array) -> np.ndarray:
    return _to_float32(ak.sum(array, axis=1))


def _match_to_top_tracks(
    track_keys,
    order,
    object_keys,
    object_values,
) -> np.ndarray:
    output = np.full((len(track_keys), N_OBJECTS), np.nan, dtype=np.float32)

    event_values = zip(
        ak.to_list(track_keys),
        ak.to_list(order),
        ak.to_list(object_keys),
        ak.to_list(object_values),
        strict=True,
    )
    for event, (tracks, top_order, keys, values) in enumerate(event_values):
        lookup = {}
        for key, value in zip(keys, values, strict=True):
            key = int(key)
            value = float(value)
            if key < 0 or not np.isfinite(value):
                continue
            if key not in lookup or value > lookup[key]:
                lookup[key] = value

        for rank, track_index in enumerate(top_order[:N_OBJECTS]):
            key = int(tracks[int(track_index)])
            if key in lookup:
                output[event, rank] = lookup[key]

    return output


def _matched_flag(
    track_keys,
    order,
    object_keys,
    track_exists: np.ndarray,
) -> np.ndarray:
    values = ak.ones_like(object_keys)
    matched = _match_to_top_tracks(track_keys, order, object_keys, values)
    return np.where(track_exists, np.nan_to_num(matched, nan=0.0), np.nan).astype(np.float32)


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0),
    )


def _pair_min(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    valid = np.isfinite(left) & np.isfinite(right)
    return np.where(valid, np.minimum(left, right), np.nan).astype(np.float32)


def _pair_max(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    valid = np.isfinite(left) & np.isfinite(right)
    return np.where(valid, np.maximum(left, right), np.nan).astype(np.float32)


def extract_features(tree) -> dict[str, np.ndarray]:
    arrays = {branch: _read(tree, branch) for branch in INPUT_BRANCHES}

    theta = arrays["TracRaw_theta"]
    phi = arrays["TracRaw_phi"]
    inv_r = arrays["TracRaw_invR"]
    bfield = ak.broadcast_arrays(arrays["Event_bFieldGevCm"], inv_r)[0]

    denominator = np.abs(inv_r * np.sin(theta))
    valid_momentum = denominator > 0
    safe_denominator = ak.where(valid_momentum, denominator, 1.0)
    momentum = ak.where(valid_momentum, np.abs(bfield) / safe_denominator, np.nan)
    sort_values = ak.where(np.isfinite(momentum), momentum, -np.inf)
    order = ak.argsort(sort_values, axis=1, ascending=False)

    track_index = ak.local_index(arrays["TracRaw_paIdx"], axis=1)
    top_track_index = _top_k(track_index, order)
    track_exists = np.isfinite(top_track_index)

    track_p = _top_k(momentum, order)
    track_theta = _top_k(theta, order)
    track_phi = _top_k(phi, order)
    track_charge = _top_k(arrays["TracRaw_charge"], order)

    finite_momentum = ak.where(np.isfinite(momentum), momentum, 0.0)
    track_px = finite_momentum * np.sin(theta) * np.cos(phi)
    track_py = finite_momentum * np.sin(theta) * np.sin(phi)
    track_pz = finite_momentum * np.cos(theta)
    track_p_sum = _event_sum(finite_momentum)
    track_momentum_imbalance = np.sqrt(
        _event_sum(track_px) ** 2
        + _event_sum(track_py) ** 2
        + _event_sum(track_pz) ** 2
    ).astype(np.float32)

    matched_em_energy = _match_to_top_tracks(
        arrays["TracRaw_paIdx"],
        order,
        arrays["EmShower_paIdx"],
        arrays["EmShower_energy"],
    )
    matched_had_energy = _match_to_top_tracks(
        arrays["TracRaw_paIdx"],
        order,
        arrays["HadShower_paIdx"],
        arrays["HadShower_energy"],
    )

    track_values = {
        "p": track_p,
        "abs_cos_theta": np.abs(np.cos(track_theta)).astype(np.float32),
        "charge": track_charge,
        "impact_rphi": _top_k(arrays["TracRaw_impactRPhi"], order),
        "impact_z": _top_k(arrays["TracRaw_impactZ"], order),
        "is_muon": _matched_flag(
            arrays["TracRaw_paIdx"],
            order,
            arrays["MuidRaw_paIdx"],
            track_exists,
        ),
        "is_electron": _matched_flag(
            arrays["TracRaw_paIdx"],
            order,
            arrays["ElidRaw_paIdx"],
            track_exists,
        ),
        "em_energy_over_p": _safe_ratio(matched_em_energy, track_p),
        "had_energy_over_p": _safe_ratio(matched_had_energy, track_p),
        "kaon_combined": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_kaonCombined"],
        ),
        "kaon_dedx": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_kaonDedx"],
        ),
        "kaon_rich": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_kaonRich"],
        ),
        "pion_rich": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_pionRich"],
        ),
        "proton_combined": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_protonCombined"],
        ),
        "proton_dedx": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_protonDedx"],
        ),
        "proton_rich": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_protonRich"],
        ),
        "rich_quality": _match_to_top_tracks(
            arrays["TracRaw_paIdx"],
            order,
            arrays["HaidRaw_paIdx"],
            arrays["HaidRaw_richQuality"],
        ),
        "tpc_dedx80": _match_to_top_tracks(
            track_index,
            order,
            arrays["MtpcRaw_tracRawIdx"],
            arrays["MtpcRaw_dEdx80Max"],
        ),
        "tpc_dedx80_sigma": _match_to_top_tracks(
            track_index,
            order,
            arrays["MtpcRaw_tracRawIdx"],
            arrays["MtpcRaw_dEdx80Sigma"],
        ),
        "tpc_dedx65": _match_to_top_tracks(
            track_index,
            order,
            arrays["MtpcRaw_tracRawIdx"],
            arrays["MtpcRaw_dEdx65Max"],
        ),
        "tpc_dedx65_sigma": _match_to_top_tracks(
            track_index,
            order,
            arrays["MtpcRaw_tracRawIdx"],
            arrays["MtpcRaw_dEdx65Sigma"],
        ),
        "tpc_dedx80_integrated": _match_to_top_tracks(
            track_index,
            order,
            arrays["MtpcRaw_tracRawIdx"],
            arrays["MtpcRaw_dEdx80Integrated"],
        ),
    }

    cos_opening = (
        np.sin(track_theta[:, 0])
        * np.sin(track_theta[:, 1])
        * np.cos(track_phi[:, 0] - track_phi[:, 1])
        + np.cos(track_theta[:, 0]) * np.cos(track_theta[:, 1])
    )
    opening_angle = np.arccos(np.clip(cos_opening, -1.0, 1.0))

    em_energy_sum = _event_sum(arrays["EmShower_energy"])
    had_energy_sum = _event_sum(arrays["HadShower_energy"])
    features = {
        "n_track": _to_float32(arrays["nTracRaw"]),
        "n_em_shower": _to_float32(arrays["nEmShower"]),
        "n_had_shower": _to_float32(arrays["nHadShower"]),
        "n_muon": _to_float32(arrays["nMuidRaw"]),
        "n_electron": _to_float32(arrays["nElidRaw"]),
        "n_hadron_pid": _to_float32(arrays["nHaidRaw"]),
        "n_vertex": _to_float32(arrays["nVtx"]),
        "em_energy_sum": em_energy_sum,
        "had_energy_sum": had_energy_sum,
        "calorimeter_energy_sum": (em_energy_sum + had_energy_sum).astype(np.float32),
        "track_p_sum": track_p_sum,
        "track_momentum_imbalance": track_momentum_imbalance,
    }
    for rank in range(N_OBJECTS):
        for name in TRACK_FEATURES:
            features[f"track{rank + 1}_{name}"] = track_values[name][:, rank]

    em_showers = _top_k_descending(arrays["EmShower_energy"])
    had_shower_order = ak.argsort(arrays["HadShower_energy"], axis=1, ascending=False)
    had_showers = _top_k(arrays["HadShower_energy"], had_shower_order)
    had_shower_n_hits = _top_k(arrays["HadShower_nHits"], had_shower_order)
    had_shower_n_hit_rows = _top_k(arrays["HadShower_nHitRows"], had_shower_order)
    for rank in range(N_OBJECTS):
        features[f"em_shower{rank + 1}_energy"] = em_showers[:, rank]
    for rank in range(N_OBJECTS):
        features[f"had_shower{rank + 1}_energy"] = had_showers[:, rank]
        features[f"had_shower{rank + 1}_n_hits"] = had_shower_n_hits[:, rank]
        features[f"had_shower{rank + 1}_n_hit_rows"] = had_shower_n_hit_rows[:, rank]

    tpc_dedx80_1 = track_values["tpc_dedx80"][:, 0]
    tpc_dedx80_2 = track_values["tpc_dedx80"][:, 1]
    kaon_dedx_1 = track_values["kaon_dedx"][:, 0]
    kaon_dedx_2 = track_values["kaon_dedx"][:, 1]
    kaon_combined_1 = track_values["kaon_combined"][:, 0]
    kaon_combined_2 = track_values["kaon_combined"][:, 1]
    two_tracks = track_exists[:, 0] & track_exists[:, 1]
    features.update(
        {
            "track12_acolinearity": (np.pi - opening_angle).astype(np.float32),
            "track12_momentum_ratio": _safe_ratio(track_p[:, 1], track_p[:, 0]),
            "track12_charge_product": (track_charge[:, 0] * track_charge[:, 1]).astype(
                np.float32
            ),
            "track12_tpc_dedx80_min": _pair_min(tpc_dedx80_1, tpc_dedx80_2),
            "track12_tpc_dedx80_max": _pair_max(tpc_dedx80_1, tpc_dedx80_2),
            "track12_tpc_dedx80_abs_diff": np.where(
                np.isfinite(tpc_dedx80_1) & np.isfinite(tpc_dedx80_2),
                np.abs(tpc_dedx80_1 - tpc_dedx80_2),
                np.nan,
            ).astype(np.float32),
            "track12_kaon_dedx_min": _pair_min(kaon_dedx_1, kaon_dedx_2),
            "track12_kaon_combined_min": _pair_min(
                kaon_combined_1,
                kaon_combined_2,
            ),
            "track12_both_have_tpc": np.where(
                two_tracks,
                np.isfinite(tpc_dedx80_1) & np.isfinite(tpc_dedx80_2),
                np.nan,
            ).astype(np.float32),
            "track12_both_have_had_pid": np.where(
                two_tracks,
                np.isfinite(kaon_combined_1) & np.isfinite(kaon_combined_2),
                np.nan,
            ).astype(np.float32),
        }
    )

    if tuple(features) != FEATURE_NAMES:
        raise RuntimeError("FEATURE_NAMES and extracted feature order differ")

    n_events = int(tree.num_entries)
    for name, values in features.items():
        if len(values) != n_events:
            raise RuntimeError(f"Feature {name} has {len(values)} rows, expected {n_events}")

    return features
