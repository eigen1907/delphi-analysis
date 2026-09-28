from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path
from typing import Any

import awkward as ak
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "bdt" / "features" / "pid.json"
SOURCE_FILE_NAME = "nanoaod_raw_sdst.root"
TREE_NAME = "Events"
N_OBJECTS = 10
N_COVARIANCE_TRACKS = 4
COVARIANCE_SIZE = 15

COLLECTION_ORDER = (
    "Event",
    "TracRaw",
    "EmShower",
    "EmLayer",
    "HadShower",
    "HadHit",
    "Stic",
    "MuidRaw",
    "ElidRaw",
    "HaidRaw",
    "VdAssocHit",
    "VdUnassocHit",
    "MtpcRaw",
    "Vtx",
)

MATCH_RULES = {
    "EmShower": ("EmShower_paIdx", "particle", "EmShower_energy"),
    "HadShower": ("HadShower_paIdx", "particle", "HadShower_energy"),
    "MuidRaw": ("MuidRaw_paIdx", "particle", None),
    "ElidRaw": ("ElidRaw_paIdx", "particle", None),
    "HaidRaw": ("HaidRaw_paIdx", "particle", None),
    "VdAssocHit": (
        "VdAssocHit_tracRawIdx",
        "track",
        "VdAssocHit_signalToNoise",
    ),
    "MtpcRaw": ("MtpcRaw_tracRawIdx", "track", None),
}

ORDER_RULES = {
    "EmLayer": "EmLayer_energy",
    "HadHit": "HadHit_energy",
    "Stic": "Stic_energyFromMain",
    "VdUnassocHit": "VdUnassocHit_signalToNoise",
    "Vtx": "Vtx_nOutgoing",
}


def load_config(config: Path | str | dict[str, Any] | None = None) -> dict[str, Any]:
    if isinstance(config, dict):
        loaded = config
    else:
        path = Path(config) if config is not None else DEFAULT_CONFIG_PATH
        with path.open() as input_file:
            loaded = json.load(input_file)

    if not isinstance(loaded.get("name"), str):
        raise ValueError("Feature config needs a string name")
    if not isinstance(loaded.get("features"), dict):
        raise ValueError("Feature config needs a features object")

    unknown = set(loaded["features"]) - set(COLLECTION_ORDER)
    if unknown:
        raise ValueError(f"Unknown ROOT collections: {sorted(unknown)}")
    if "TracRaw" not in loaded["features"]:
        raise ValueError("Feature config needs a TracRaw collection")
    for collection, branches in loaded["features"].items():
        if not isinstance(branches, list) or not branches:
            raise ValueError(
                f"Feature collection {collection} needs a non-empty branch list"
            )
    return loaded


def _is_scalar(branch: str) -> bool:
    return branch.startswith("Event_") or branch.startswith("n")


def feature_names(
    config: Path | str | dict[str, Any] | None = None,
) -> tuple[str, ...]:
    features = load_config(config)["features"]
    names: list[str] = []
    for collection in COLLECTION_ORDER:
        for branch in features.get(collection, ()):
            if _is_scalar(branch):
                names.append(branch)
            elif branch == "TracRaw_weightMatrix":
                names.extend(
                f"{branch}__{track}__{element}"
                for track in range(N_COVARIANCE_TRACKS)
                for element in range(COVARIANCE_SIZE)
                )
            else:
                names.extend(
                    f"{branch}__{index}" for index in range(N_OBJECTS)
                )
    return tuple(names)


def input_branches(
    config: Path | str | dict[str, Any] | None = None,
) -> tuple[str, ...]:
    features = load_config(config)["features"]
    branches: list[str] = ["TracRaw_invR"]
    for collection in COLLECTION_ORDER:
        if collection not in features:
            continue
        branches.extend(features[collection])
        if collection in MATCH_RULES:
            source_id, target, source_sort = MATCH_RULES[collection]
            branches.append(source_id)
            if target == "particle":
                branches.append("TracRaw_paIdx")
            if source_sort is not None:
                branches.append(source_sort)
        elif collection in ORDER_RULES:
            branches.append(ORDER_RULES[collection])
    return tuple(dict.fromkeys(branches))


def _padded_matrix(array: ak.Array, width: int) -> np.ndarray:
    values = ak.values_astype(array, np.float32)
    padded = ak.pad_none(values, width, axis=1, clip=True)
    return np.asarray(ak.to_numpy(ak.fill_none(padded, np.nan)), dtype=np.float32)


def _padded_ids(array: ak.Array, width: int) -> np.ndarray:
    padded = ak.pad_none(array, width, axis=1, clip=True)
    matrix = ak.to_numpy(padded, allow_missing=True)
    return np.asarray(np.ma.filled(matrix, -1), dtype=np.int64)


def _descending_order(array: ak.Array) -> ak.Array:
    return ak.argsort(array, axis=1, ascending=False, stable=True)


def _ordered_matrices(
    arrays: dict[str, ak.Array],
    branches: list[str],
    order: ak.Array,
    width: int = N_OBJECTS,
) -> dict[str, np.ndarray]:
    return {
        branch: _padded_matrix(arrays[branch][order], width)
        for branch in branches
    }


def _matched_matrices(
    arrays: dict[str, ak.Array],
    branches: list[str],
    source_id_branch: str,
    target_ids: np.ndarray,
    source_sort_branch: str | None = None,
) -> dict[str, np.ndarray]:
    source_order = None
    if source_sort_branch is not None:
        source_order = ak.argsort(
            arrays[source_sort_branch],
            axis=1,
            ascending=True,
            stable=True,
        )

    source_ids = arrays[source_id_branch]
    if source_id_branch == "HaidRaw_paIdx":
        source_ids = source_ids - 1
    if source_order is not None:
        source_ids = source_ids[source_order]

    source_width = int(ak.max(ak.num(source_ids, axis=1), initial=0))
    empty = np.full(target_ids.shape, np.nan, dtype=np.float32)
    if source_width == 0:
        return {branch: empty.copy() for branch in branches}

    source_id_matrix = _padded_ids(source_ids, source_width)
    valid_ids = source_id_matrix[source_id_matrix >= 0]
    if len(valid_ids) == 0:
        return {branch: empty.copy() for branch in branches}

    n_events = len(target_ids)
    dense_width = int(valid_ids.max()) + 1
    rows = np.arange(n_events)
    row_matrix = np.broadcast_to(rows[:, np.newaxis], target_ids.shape)
    source_positions = np.full((n_events, dense_width), -1, dtype=np.int32)
    for index in range(source_width):
        source_index = source_id_matrix[:, index]
        valid = (source_index >= 0) & (source_index < dense_width)
        source_positions[rows[valid], source_index[valid]] = index

    positions = np.full(target_ids.shape, -1, dtype=np.int32)
    valid_target = (target_ids >= 0) & (target_ids < dense_width)
    positions[valid_target] = source_positions[
        row_matrix[valid_target],
        target_ids[valid_target],
    ]
    valid_match = positions >= 0
    matrices = {}
    for branch in branches:
        values = arrays[branch]
        if source_order is not None:
            values = values[source_order]
        value_matrix = _padded_matrix(values, source_width)
        matched = empty.copy()
        matched[valid_match] = value_matrix[
            row_matrix[valid_match],
            positions[valid_match],
        ]
        matrices[branch] = matched
    return matrices


def extract_features(
    tree: Any,
    config: Path | str | dict[str, Any] | None = None,
) -> OrderedDict[str, np.ndarray]:
    loaded = load_config(config)
    features = loaded["features"]
    branches = input_branches(loaded)
    missing = [branch for branch in branches if branch not in tree]
    if missing:
        raise KeyError(f"Missing RAW-SDST branches: {missing}")

    arrays = {branch: tree[branch].array(library="ak") for branch in branches}
    matrices: dict[str, np.ndarray] = {}

    track_order = ak.argsort(
        abs(arrays["TracRaw_invR"]),
        axis=1,
        ascending=True,
        stable=True,
    )
    track_branches = [
        branch
        for branch in features["TracRaw"]
        if not _is_scalar(branch) and branch != "TracRaw_weightMatrix"
    ]
    matrices.update(
        _ordered_matrices(arrays, track_branches, track_order)
    )

    covariance = None
    if "TracRaw_weightMatrix" in features["TracRaw"]:
        values = arrays["TracRaw_weightMatrix"][track_order]
        padded = ak.pad_none(
            values,
            N_COVARIANCE_TRACKS,
            axis=1,
            clip=True,
        )
        covariance = np.asarray(
            np.ma.filled(ak.to_numpy(padded, allow_missing=True), np.nan),
            dtype=np.float32,
        )
        if covariance.shape[2] != COVARIANCE_SIZE:
            raise ValueError("TracRaw_weightMatrix does not contain 15 elements")

    leading_track_indices = _padded_ids(
        ak.local_index(arrays["TracRaw_invR"], axis=1)[track_order],
        N_OBJECTS,
    )
    needs_particle_ids = any(
        collection in features and rule[1] == "particle"
        for collection, rule in MATCH_RULES.items()
    )
    leading_particle_ids = None
    if needs_particle_ids:
        leading_particle_ids = _padded_ids(
            arrays["TracRaw_paIdx"][track_order],
            N_OBJECTS,
        )

    targets = {
        "track": leading_track_indices,
        "particle": leading_particle_ids,
    }
    for collection in COLLECTION_ORDER:
        if collection in {"Event", "TracRaw"} or collection not in features:
            continue
        vector_branches = [
            branch
            for branch in features[collection]
            if not _is_scalar(branch)
        ]
        if not vector_branches:
            continue
        if collection in MATCH_RULES:
            source_id, target, source_sort = MATCH_RULES[collection]
            collection_matrices = _matched_matrices(
                arrays,
                vector_branches,
                source_id,
                targets[target],
                source_sort,
            )
        else:
            order = _descending_order(arrays[ORDER_RULES[collection]])
            collection_matrices = _ordered_matrices(
                arrays,
                vector_branches,
                order,
            )
        matrices.update(collection_matrices)

    output: OrderedDict[str, np.ndarray] = OrderedDict()
    for collection in COLLECTION_ORDER:
        for branch in features.get(collection, ()):
            if _is_scalar(branch):
                output[branch] = np.asarray(arrays[branch], dtype=np.float32)
            elif branch == "TracRaw_weightMatrix":
                for track in range(N_COVARIANCE_TRACKS):
                    for element in range(COVARIANCE_SIZE):
                        output[f"{branch}__{track}__{element}"] = covariance[
                            :, track, element
                        ]
            else:
                for index in range(N_OBJECTS):
                    output[f"{branch}__{index}"] = matrices[branch][:, index]

    if tuple(output) != feature_names(loaded):
        raise RuntimeError("Extracted feature order does not match the feature config")
    return output
