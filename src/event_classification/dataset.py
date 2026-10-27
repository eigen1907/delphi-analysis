"""Load aligned analysis-level and RAW event datasets for classification."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import uproot

from common.samples import sample_name
from event_classification.bdt import (
    extract_paired_features,
    paired_feature_config,
    paired_feature_names,
)


LABELS = {"Zee": 0, "Zmumu": 1, "ZKK": 2, "Zpipi": 3, "Ztautau": 4}


def load_metadata(dataset_dir: Path) -> dict[str, Any]:
    with (dataset_dir / "metadata.json").open() as metadata_file:
        metadata = json.load(metadata_file)
    if metadata.get("format_version") not in (5, 6):
        raise ValueError(
            "BDT dataset uses an older feature format; rerun "
            "scripts/01_sample_preparation/prepare_dataset.sh"
        )
    return metadata


def feature_metadata(metadata: dict, config=None) -> dict:
    """Describe model inputs independently of the reusable paired dataset."""
    if metadata["format_version"] == 5:
        return metadata
    loaded = paired_feature_config(config)
    return {
        **metadata,
        "feature_set": loaded["name"],
        "feature_config": loaded,
        "features": list(paired_feature_names(loaded)),
        "labels": LABELS,
    }


def load_dataset_split(dataset_dir: Path, split: str, metadata: dict):
    """Read one split, joining aligned RAW and analysis-level feature columns."""
    if metadata["format_version"] == 5:
        return load_split(
            dataset_dir / f"{split}.root", metadata["features"], metadata["tree"]
        )

    paths = []
    for sample_dir, sample in metadata["samples"].items():
        name = sample.get("class", sample.get("sample", sample_name(sample_dir)))
        if name in LABELS:
            paths.append((name, dataset_dir / sample["splits"][split]["file"]))
    if {name for name, _ in paths} != set(LABELS):
        raise ValueError("Paired dataset must contain all configured BDT classes")
    paths.sort(key=lambda item: item[0])

    n_events = 0
    for _, path in paths:
        with uproot.open(path) as root_file:
            n_events += root_file[metadata["trees"]["raw_sdst"]].num_entries
    features = np.empty((n_events, len(metadata["features"])), dtype=np.float32)
    labels = np.empty(n_events, dtype=np.int8)
    offset = 0
    for name, path in paths:
        with uproot.open(path) as root_file:
            extracted = extract_paired_features(
                root_file[metadata["trees"]["raw_sdst"]],
                root_file[metadata["trees"]["analysisLevel"]],
                metadata["feature_config"],
            )
            count = len(next(iter(extracted.values())))
            for column, feature in enumerate(metadata["features"]):
                features[offset : offset + count, column] = extracted[feature]
            labels[offset : offset + count] = LABELS[name]
            offset += count
        print(f"{split}: {name}: {count} events", flush=True)
    return features, labels


def load_split(
    path: Path,
    selected_features: list[str],
    tree_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    with uproot.open(path) as root_file:
        tree = root_file[tree_name]
        features = [
            np.asarray(tree[name].array(library="np"), dtype=np.float32)
            for name in selected_features
        ]
        labels = np.asarray(tree["label"].array(library="np"), dtype=np.int8)
    return np.column_stack(features), labels
