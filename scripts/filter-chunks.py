#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
from collections import Counter
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))
SAMPLE_SET = "100kTest"
INPUT_ROOT = PROJECT_ROOT / "output" / "chunks_raw" / SAMPLE_SET
OUTPUT_ROOT = PROJECT_ROOT / "output" / "chunks" / SAMPLE_SET
CHECK_PATH = PROJECT_ROOT / "data" / "check" / SAMPLE_SET / "event-filter.csv"

ROOT_FILES = (
    "nanoaod.root",
    "nanoaod_raw_sdst.root",
    "nanoaod_raw_fadana.root",
    "nanoaod_ttree.root",
    "pythiastandalone.root",
)
EVENT_BRANCHES = ("Event_eventNumber", "Event_evtNumber", "EventNo")
RUN_BRANCH = "Event_runNumber"


def sample_name(sample_dir: Path) -> str:
    return sample_dir.name.split("_")[-1]


def job_number(job_dir: Path) -> int:
    try:
        return int(job_dir.name.removeprefix("job_"))
    except ValueError:
        return 10**9


def is_tree(obj) -> bool:
    return hasattr(obj, "arrays") and hasattr(obj, "keys") and hasattr(obj, "num_entries")


def find_event_branch(tree) -> str | None:
    branches = set(tree.keys())
    for branch in EVENT_BRANCHES:
        if branch in branches:
            return branch
    return None


def object_names(root_file) -> list[str]:
    names = list(root_file.keys(cycle=False))
    if "Events" in names:
        names.remove("Events")
        names.insert(0, "Events")
    return names


def event_mask(events, valid_events: set[int]) -> tuple[np.ndarray, int, int]:
    keep = np.zeros(len(events), dtype=bool)
    seen = set()
    n_invalid = 0
    n_duplicated = 0

    for i, event_number in enumerate(events):
        event_number = int(event_number)
        if event_number not in valid_events:
            n_invalid += 1
            continue
        if event_number in seen:
            n_duplicated += 1
            continue

        seen.add(event_number)
        keep[i] = True

    return keep, n_invalid, n_duplicated


def event_keys(arrays, event_branch: str) -> list[tuple[int, int, int]]:
    return [
        (int(run), int(event), int(n_gen))
        for run, event, n_gen in zip(
            ak.to_numpy(arrays[RUN_BRANCH]),
            ak.to_numpy(arrays[event_branch]),
            ak.to_numpy(arrays["nGenPart"]),
            strict=True,
        )
    ]


def valid_mask(arrays, event_branch: str) -> tuple[np.ndarray, int, int]:
    keys = event_keys(arrays, event_branch)
    has_gen = np.asarray(ak.to_numpy(arrays["nGenPart"]) > 0, dtype=bool)
    counts = Counter(key for key, keep in zip(keys, has_gen, strict=True) if keep)
    keep = np.asarray([has_gen[i] and counts[key] == 1 for i, key in enumerate(keys)], dtype=bool)
    n_invalid = int(np.count_nonzero(~has_gen))
    n_duplicated = sum(count - 1 for count in counts.values() if count > 1)
    return keep, n_invalid, n_duplicated


def filter_tree(tree, sample: str, job: str, file_name: str, tree_name: str, valid_events, rows):
    event_branch = find_event_branch(tree)
    if event_branch is None:
        return tree.arrays(library="ak")

    arrays = tree.arrays(library="ak")

    if RUN_BRANCH in tree.keys() and "nGenPart" in tree.keys():
        keep, n_invalid, n_duplicated = valid_mask(arrays, event_branch)
        valid_events = set(int(event) for event in ak.to_numpy(arrays[event_branch][keep]))
    elif valid_events is not None:
        events = np.asarray(ak.to_numpy(arrays[event_branch]))
        keep, n_invalid, n_duplicated = event_mask(events, valid_events)
    else:
        return arrays

    rows.append((sample, job, file_name, tree_name, int(tree.num_entries), n_invalid, n_duplicated))

    return arrays[keep], valid_events


def filter_root_file(src: Path, dst: Path, sample: str, job: str, file_name: str, valid_events, rows) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)

    with uproot.open(src) as root_file, uproot.recreate(dst) as output_file:
        for name in object_names(root_file):
            obj = root_file[name]
            if not is_tree(obj):
                output_file[name] = obj
                continue

            result = filter_tree(
                obj,
                sample=sample,
                job=job,
                file_name=file_name,
                tree_name=name,
                valid_events=valid_events,
                rows=rows,
            )
            if isinstance(result, tuple):
                arrays, valid_events = result
            else:
                arrays = result
            output_file[name] = arrays


def main() -> None:
    rows = []

    for sample_dir in sorted(path for path in INPUT_ROOT.iterdir() if path.is_dir()):
        sample = sample_name(sample_dir)
        job_dirs = sorted((sample_dir / "final_root").glob("job_*"), key=job_number)

        for job_dir in job_dirs:
            for file_name in ROOT_FILES:
                src = job_dir / file_name
                if not src.exists():
                    continue

                dst = OUTPUT_ROOT / sample_dir.name / "final_root" / job_dir.name / file_name
                filter_root_file(src, dst, sample, job_dir.name, file_name, None, rows)

    CHECK_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CHECK_PATH.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("sample", "job", "file", "tree", "n_event", "n_invalid_event", "n_duplicated_event"))
        writer.writerows(rows)

    print(f"filtered chunks: {OUTPUT_ROOT}")
    print(f"check summary: {CHECK_PATH}")


if __name__ == "__main__":
    main()
