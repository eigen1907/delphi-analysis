#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import uproot


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[2]))

from bdt import (
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
    return sample_dir.name.split("_")[-1]


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path, help="filtered chunk root")
    parser.add_argument("-o", "--output", required=True, type=Path, help="BDT dataset output directory")
    parser.add_argument(
        "--fractions",
        nargs=3,
        type=float,
        default=(0.6, 0.2, 0.2),
        metavar=("TRAIN", "VAL", "TEST"),
        help="job split fractions (default: 0.6 0.2 0.2)",
    )
    parser.add_argument("--seed", type=int, default=1907, help="split seed (default: 1907)")
    parser.add_argument(
        "--feature-config",
        type=Path,
        default=PROJECT_ROOT / "config" / "bdt" / "features" / "pid.json",
        help="RAW feature JSON config (default: config/bdt/features/pid.json)",
    )
    args = parser.parse_args()

    feature_config = load_config(args.feature_config)
    feature_set = feature_config["name"]
    fractions = tuple(args.fractions)
    if any(value <= 0 for value in fractions):
        parser.error("split fractions must be positive")
    if not np.isclose(sum(fractions), 1.0):
        parser.error("split fractions must sum to 1")
    jobs_by_class = discover_jobs(args.input, SOURCE_FILE_NAME)
    split_paths = split_jobs(jobs_by_class, fractions, args.seed)
    selected_features = feature_names(feature_config)

    args.output.mkdir(parents=True, exist_ok=True)
    split_metadata = {}
    for split in SPLIT_NAMES:
        output_path = args.output / f"{split}.root"
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
        "format_version": 4,
        "source_file": SOURCE_FILE_NAME,
        "tree": TREE_NAME,
        "feature_set": feature_set,
        "feature_config": feature_config,
        "seed": args.seed,
        "fractions": dict(zip(SPLIT_NAMES, fractions, strict=True)),
        "labels": LABELS,
        "features": list(selected_features),
        "splits": split_metadata,
    }
    metadata_path = args.output / "metadata.json"
    with metadata_path.open("w") as output_file:
        json.dump(metadata, output_file, indent=2, sort_keys=True)
        output_file.write("\n")
    print(f"metadata: {metadata_path}")


if __name__ == "__main__":
    main()
