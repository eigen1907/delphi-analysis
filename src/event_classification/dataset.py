"""Build job-separated training, validation, and test feature datasets."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import uproot

from common.samples import sample_name
from event_classification.bdt import (
    SOURCE_FILE_NAME,
    TREE_NAME,
    extract_features,
    feature_names,
    load_config,
)


LABELS = {
    "Zee": 0,
    "Zmumu": 1,
    "ZKK": 2,
    "Zpipi": 3,
}
SPLIT_NAMES = ("train", "val", "test")


def class_name(sample_dir: Path) -> str:
    return sample_name(sample_dir)


def job_number(job_dir: Path) -> int:
    match = re.search(r"(\d+)$", job_dir.name)
    if match is None:
        raise ValueError(f"Could not read job number from {job_dir}")
    return int(match.group(1))


def discover_jobs(input_root: Path, source_file: str) -> dict[str, list[Path]]:
    jobs_by_class = {}
    for sample_dir in sorted(path for path in input_root.iterdir() if path.is_dir()):
        sample = class_name(sample_dir)
        if sample not in LABELS:
            continue
        if sample in jobs_by_class:
            raise ValueError(f"Multiple input directories resolve to class {sample}")

        jobs = []
        for job_dir in sorted((sample_dir / "final_root").glob("job_*"), key=job_number):
            source_path = job_dir / source_file
            if source_path.exists():
                jobs.append(source_path)
        if not jobs:
            raise FileNotFoundError(f"No {source_file} chunks found under {sample_dir}")
        jobs_by_class[sample] = jobs

    missing = set(LABELS) - set(jobs_by_class)
    if missing:
        raise FileNotFoundError(f"Missing BDT classes under {input_root}: {sorted(missing)}")
    return jobs_by_class


def split_jobs(
    jobs_by_class: dict[str, list[Path]],
    fractions: tuple[float, float, float],
    seed: int,
) -> dict[str, list[tuple[str, Path]]]:
    rng = np.random.default_rng(seed)
    split_paths = {name: [] for name in SPLIT_NAMES}

    for sample in LABELS:
        jobs = list(jobs_by_class[sample])
        rng.shuffle(jobs)
        n_jobs = len(jobs)
        n_train = int(n_jobs * fractions[0])
        n_val = int(n_jobs * fractions[1])
        boundaries = (0, n_train, n_train + n_val, n_jobs)

        for index, split in enumerate(SPLIT_NAMES):
            selected = jobs[boundaries[index] : boundaries[index + 1]]
            if not selected:
                raise ValueError(f"Split {split} has no {sample} jobs")
            split_paths[split].extend((sample, path) for path in selected)

    for paths in split_paths.values():
        rng.shuffle(paths)
    return split_paths


def output_branch_types(selected_features: tuple[str, ...]) -> dict[str, object]:
    return {
        "run_number": np.int32,
        "event_number": np.int32,
        "job_id": np.int16,
        "label": np.int8,
        **{feature: np.float32 for feature in selected_features},
    }


def write_split(
    output_path: Path,
    paths: list[tuple[str, Path]],
    selected_features: tuple[str, ...],
    feature_config: dict[str, Any],
) -> tuple[dict[str, int], dict[str, list[int]]]:
    event_counts = {sample: 0 for sample in LABELS}
    jobs = {sample: [] for sample in LABELS}

    with uproot.recreate(output_path) as output_file:
        output_tree = output_file.mktree(
            TREE_NAME,
            output_branch_types(selected_features),
        )

        for sample, source_path in paths:
            job_id = job_number(source_path.parent)
            with uproot.open(source_path) as source_file:
                if TREE_NAME not in source_file:
                    raise KeyError(f"Missing tree {TREE_NAME} in {source_path}")
                tree = source_file[TREE_NAME]
                features = extract_features(tree, feature_config)
                n_events = int(tree.num_entries)
                data = {
                    "run_number": np.asarray(tree["Event_runNumber"].array(library="np"), dtype=np.int32),
                    "event_number": np.asarray(tree["Event_eventNumber"].array(library="np"), dtype=np.int32),
                    "job_id": np.full(n_events, job_id, dtype=np.int16),
                    "label": np.full(n_events, LABELS[sample], dtype=np.int8),
                    **features,
                }
                output_tree.extend(data)

            event_counts[sample] += n_events
            jobs[sample].append(job_id)

    for sample in jobs:
        jobs[sample].sort()
    return event_counts, jobs


def prepare_dataset(
    input_root: Path,
    output_root: Path,
    feature_config_path: Path,
    fractions: tuple[float, float, float] = (0.6, 0.2, 0.2),
    seed: int = 1907,
) -> dict[str, Any]:
    """Build feature datasets with disjoint jobs for training, validation, and test."""
    input_root = Path(input_root)
    output_root = Path(output_root)
    if len(fractions) != 3 or any(value <= 0 for value in fractions):
        raise ValueError("Three positive split fractions are required")
    if not np.isclose(sum(fractions), 1.0):
        raise ValueError("Split fractions must sum to 1")

    feature_config = load_config(feature_config_path)
    jobs_by_class = discover_jobs(input_root, SOURCE_FILE_NAME)
    split_paths = split_jobs(jobs_by_class, fractions, seed)
    selected_features = feature_names(feature_config)

    output_root.mkdir(parents=True, exist_ok=True)
    split_metadata = {}
    for split in SPLIT_NAMES:
        output_path = output_root / f"{split}.root"
        event_counts, jobs = write_split(
            output_path,
            split_paths[split],
            selected_features,
            feature_config,
        )
        split_metadata[split] = {
            "file": output_path.name,
            "event_counts": event_counts,
            "jobs": jobs,
        }
        print(f"{split}: {sum(event_counts.values())} events -> {output_path}")

    metadata = {
        "format_version": 5,
        "source_file": SOURCE_FILE_NAME,
        "tree": TREE_NAME,
        "feature_set": feature_config["name"],
        "feature_config": feature_config,
        "seed": seed,
        "fractions": dict(zip(SPLIT_NAMES, fractions, strict=True)),
        "labels": LABELS,
        "features": list(selected_features),
        "splits": split_metadata,
    }
    metadata_path = output_root / "metadata.json"
    with metadata_path.open("w") as output_file:
        json.dump(metadata, output_file, indent=2, sort_keys=True)
        output_file.write("\n")
    print(f"metadata: {metadata_path}")
    return metadata


def load_metadata(dataset_dir: Path) -> dict[str, Any]:
    with (dataset_dir / "metadata.json").open() as metadata_file:
        metadata = json.load(metadata_file)
    if metadata.get("format_version") != 5:
        raise ValueError(
            "BDT dataset uses an older feature format; rerun "
            "scripts/02_sample_preparation/prepare_dataset.py"
        )
    return metadata


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
