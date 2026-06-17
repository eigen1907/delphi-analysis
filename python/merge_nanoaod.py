from __future__ import annotations

import fnmatch
from collections import Counter
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from utils import META_TREE_NAME, TREE_NAME


RUN = "Event_runNumber"
EVENT = "Event_eventNumber"
EVENT_ALIASES = ("Event_eventNumber", "Event_evtNumber")
RAW_DROP = ("Event_*", "nGenPart", "GenPart", "GenPart_*")
P4_CANDIDATES = ("GenPart_fourMomentum", "GenPart_vector")
PREFIXES = {"sdst": "SDST_", "raw_sdst": "RAWSDST_", "raw_fadana": "RAWFADANA_"}
ROUND_DECIMALS = 6


def unique(items: list[str]) -> list[str]:
    return list(dict.fromkeys(items))


def parent_branch(name: str, keys: set[str]) -> str:
    parts = name.split(".")
    for stop in range(1, len(parts)):
        parent = ".".join(parts[:stop])
        if parent in keys:
            return parent
    return name


def selected_branches(path: Path, raw: bool = False) -> tuple[list[str], list[str]]:
    with uproot.open(path) as root_file:
        keys = list(root_file[TREE_NAME].keys())

    selected = []
    key_set = set(keys)
    for key in keys:
        if key in {RUN, *EVENT_ALIASES}:
            continue
        if raw and any(fnmatch.fnmatchcase(key, pattern) for pattern in RAW_DROP):
            continue
        selected.append(parent_branch(key, key_set))
    return unique(selected), keys


def event_branch(keys: list[str]) -> str:
    for name in EVENT_ALIASES:
        if name in keys:
            return name
    raise RuntimeError(f"Missing event number branch. Tried {EVENT_ALIASES}")


def p4_branch(keys: list[str]) -> str:
    for name in P4_CANDIDATES:
        if name in keys:
            return name
    raise RuntimeError(f"Missing GenPart four-momentum branch. Tried {P4_CANDIDATES}")


def read_sample(path: Path, branches: list[str], keys: list[str]) -> dict[str, object]:
    event_name = event_branch(keys)
    p4_name = p4_branch(keys)
    read_names = [event_name, RUN, "nGenPart", f"{p4_name}.*"]

    for branch in branches:
        read_names.append(f"{branch}.*" if any(key.startswith(f"{branch}.") for key in keys) else branch)

    with uproot.open(path) as root_file:
        arrays = root_file[TREE_NAME].arrays(unique(read_names), library="ak", how=dict)

    p4 = arrays[p4_name]
    return {
        "arrays": {branch: arrays[branch] for branch in branches if branch in arrays},
        "branches": [branch for branch in branches if branch in arrays],
        "run": np.asarray(ak.to_numpy(arrays[RUN]), dtype=np.int64),
        "event": np.asarray(ak.to_numpy(arrays[event_name]), dtype=np.int64),
        "n_genpart": np.asarray(ak.to_numpy(arrays["nGenPart"]), dtype=np.int64),
        "px_sum": np.round(ak.to_numpy(ak.sum(p4.fCoordinates.fX, axis=1)), ROUND_DECIMALS),
        "py_sum": np.round(ak.to_numpy(ak.sum(p4.fCoordinates.fY, axis=1)), ROUND_DECIMALS),
        "pz_sum": np.round(ak.to_numpy(ak.sum(p4.fCoordinates.fZ, axis=1)), ROUND_DECIMALS),
        "e_sum": np.round(ak.to_numpy(ak.sum(p4.fCoordinates.fT, axis=1)), ROUND_DECIMALS),
    }


def event_keys(sample: dict[str, object]) -> list[tuple[int, int, int, float, float, float, float]]:
    return [
        (int(run), int(event), int(n_genpart), float(px), float(py), float(pz), float(e))
        for run, event, n_genpart, px, py, pz, e in zip(
            sample["run"],
            sample["event"],
            sample["n_genpart"],
            sample["px_sum"],
            sample["py_sum"],
            sample["pz_sum"],
            sample["e_sum"],
            strict=True,
        )
    ]


def valid_mask(sample: dict[str, object]) -> np.ndarray:
    keys = event_keys(sample)
    has_gen = sample["n_genpart"] > 0
    counts = Counter(key for key, keep in zip(keys, has_gen, strict=True) if keep)
    return np.asarray([keep and counts[key] == 1 for key, keep in zip(keys, has_gen, strict=True)], dtype=bool)


def key_to_index(sample: dict[str, object], mask: np.ndarray) -> dict[tuple, int]:
    return {key: index for index, key in enumerate(event_keys(sample)) if mask[index]}


def take(sample: dict[str, object], indexer: np.ndarray) -> dict[str, object]:
    return {
        "arrays": {name: array[indexer] for name, array in sample["arrays"].items()},
        "branches": sample["branches"],
        "run": sample["run"][indexer],
        "event": sample["event"][indexer],
        "n_genpart": sample["n_genpart"][indexer],
    }


def add_prefixed(output: dict[str, object], sample: dict[str, object], prefix: str) -> None:
    for branch in sample["branches"]:
        output[f"{prefix}{branch}"] = sample["arrays"][branch]


def meta_scalar(value: int) -> np.ndarray:
    return np.asarray([value], dtype=np.int64)


def meta_vector(values) -> ak.Array:
    return ak.Array([values])


def duplicate_count(sample: dict[str, object]) -> int:
    has_gen = sample["n_genpart"] > 0
    counts = Counter(key for key, keep in zip(event_keys(sample), has_gen, strict=True) if keep)
    return sum(count - 1 for count in counts.values() if count > 1)


def merge_nanoaods(sdst_path: Path, raw_sdst_path: Path, raw_fadana_path: Path, output_path: Path) -> dict[str, int]:
    sdst_branches, sdst_keys = selected_branches(sdst_path)
    raw_sdst_branches, raw_sdst_keys = selected_branches(raw_sdst_path, raw=True)
    raw_fadana_branches, raw_fadana_keys = selected_branches(raw_fadana_path, raw=True)

    sdst = read_sample(sdst_path, sdst_branches, sdst_keys)
    raw_sdst = read_sample(raw_sdst_path, raw_sdst_branches, raw_sdst_keys)
    raw_fadana = read_sample(raw_fadana_path, raw_fadana_branches, raw_fadana_keys)

    sdst_valid = valid_mask(sdst)
    raw_sdst_valid = valid_mask(raw_sdst)
    raw_fadana_valid = valid_mask(raw_fadana)
    raw_sdst_index = key_to_index(raw_sdst, raw_sdst_valid)
    raw_fadana_index = key_to_index(raw_fadana, raw_fadana_valid)

    sdst_keys_all = event_keys(sdst)
    sdst_merged = np.asarray(
        [valid and key in raw_sdst_index and key in raw_fadana_index for key, valid in zip(sdst_keys_all, sdst_valid, strict=True)],
        dtype=bool,
    )
    sdst_index = np.flatnonzero(sdst_merged)
    output_keys = [sdst_keys_all[index] for index in sdst_index]
    raw_sdst_indexer = np.asarray([raw_sdst_index[key] for key in output_keys], dtype=np.int64)
    raw_fadana_indexer = np.asarray([raw_fadana_index[key] for key in output_keys], dtype=np.int64)

    raw_sdst_merged = np.zeros(len(raw_sdst["event"]), dtype=bool)
    raw_fadana_merged = np.zeros(len(raw_fadana["event"]), dtype=bool)
    raw_sdst_merged[raw_sdst_indexer] = True
    raw_fadana_merged[raw_fadana_indexer] = True

    sdst_out = take(sdst, sdst_index)
    raw_sdst_out = take(raw_sdst, raw_sdst_indexer)
    raw_fadana_out = take(raw_fadana, raw_fadana_indexer)

    events = {RUN: sdst_out["run"], EVENT: sdst_out["event"]}
    add_prefixed(events, sdst_out, PREFIXES["sdst"])
    add_prefixed(events, raw_sdst_out, PREFIXES["raw_sdst"])
    add_prefixed(events, raw_fadana_out, PREFIXES["raw_fadana"])

    meta = {
        "n_sdst_input": meta_scalar(len(sdst["event"])),
        "n_sdst_invalid_genpart0": meta_scalar(np.count_nonzero(sdst["n_genpart"] == 0)),
        "n_sdst_duplicate_after_gen_filter": meta_scalar(duplicate_count(sdst)),
        "n_sdst_valid": meta_scalar(np.count_nonzero(sdst_valid)),
        "n_sdst_merged": meta_scalar(np.count_nonzero(sdst_merged)),
        "n_raw_sdst_input": meta_scalar(len(raw_sdst["event"])),
        "n_raw_sdst_invalid_genpart0": meta_scalar(np.count_nonzero(raw_sdst["n_genpart"] == 0)),
        "n_raw_sdst_duplicate_after_gen_filter": meta_scalar(duplicate_count(raw_sdst)),
        "n_raw_sdst_valid": meta_scalar(np.count_nonzero(raw_sdst_valid)),
        "n_raw_sdst_merged": meta_scalar(np.count_nonzero(raw_sdst_merged)),
        "n_raw_fadana_input": meta_scalar(len(raw_fadana["event"])),
        "n_raw_fadana_invalid_genpart0": meta_scalar(np.count_nonzero(raw_fadana["n_genpart"] == 0)),
        "n_raw_fadana_duplicate_after_gen_filter": meta_scalar(duplicate_count(raw_fadana)),
        "n_raw_fadana_valid": meta_scalar(np.count_nonzero(raw_fadana_valid)),
        "n_raw_fadana_merged": meta_scalar(np.count_nonzero(raw_fadana_merged)),
        "sdst_run_number": meta_vector(sdst["run"]),
        "sdst_event_number": meta_vector(sdst["event"]),
        "sdst_n_genpart": meta_vector(sdst["n_genpart"]),
        "sdst_event_valid": meta_vector(sdst_valid),
        "sdst_event_merged": meta_vector(sdst_merged),
        "raw_sdst_run_number": meta_vector(raw_sdst["run"]),
        "raw_sdst_event_number": meta_vector(raw_sdst["event"]),
        "raw_sdst_n_genpart": meta_vector(raw_sdst["n_genpart"]),
        "raw_sdst_event_valid": meta_vector(raw_sdst_valid),
        "raw_sdst_event_merged": meta_vector(raw_sdst_merged),
        "raw_fadana_run_number": meta_vector(raw_fadana["run"]),
        "raw_fadana_event_number": meta_vector(raw_fadana["event"]),
        "raw_fadana_n_genpart": meta_vector(raw_fadana["n_genpart"]),
        "raw_fadana_event_valid": meta_vector(raw_fadana_valid),
        "raw_fadana_event_merged": meta_vector(raw_fadana_merged),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as root_file:
        root_file[TREE_NAME] = events
        root_file[META_TREE_NAME] = meta

    return {
        "events": len(sdst_out["event"]),
        "sdst_valid": int(np.count_nonzero(sdst_valid)),
        "sdst_merged": int(np.count_nonzero(sdst_merged)),
        "raw_sdst_valid": int(np.count_nonzero(raw_sdst_valid)),
        "raw_sdst_merged": int(np.count_nonzero(raw_sdst_merged)),
        "raw_fadana_valid": int(np.count_nonzero(raw_fadana_valid)),
        "raw_fadana_merged": int(np.count_nonzero(raw_fadana_merged)),
    }
