"""Align native and raw SDST events and create per-sample ROOT datasets."""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from common.paths import DATA_ROOT, DATASET_ROOT
from common.samples import sample_name


SPLITS = ("train", "val", "test")
SOURCES = {
    "analysisLevel": ("nanoaod.root", "Event_evtNumber"),
    "raw_sdst": ("nanoaod_raw_sdst.root", "Event_eventNumber"),
}
TREES = {source: f"{source}/Events" for source in SOURCES}


def event_indices(tree, event_branch: str) -> tuple[dict, dict]:
    """Find unique MC event keys, excluding bookkeeping and invalid events."""
    branches = ["Event_runNumber", event_branch, "nGenPart"]
    if "Event_isMC" in tree:
        branches.append("Event_isMC")
    arrays = tree.arrays(branches, library="np")
    valid = arrays["nGenPart"] > 0
    if "Event_isMC" in arrays:
        valid &= arrays["Event_isMC"] == 1
    rows = np.flatnonzero(valid)
    keys = [(int(arrays["Event_runNumber"][i]), int(arrays[event_branch][i])) for i in rows]
    counts = Counter(keys)
    indices = {key: int(i) for key, i in zip(keys, rows, strict=True) if counts[key] == 1}
    return indices, {
        "input_rows": int(tree.num_entries),
        "invalid_rows": int(np.count_nonzero(~valid)),
        "duplicate_rows": sum(count for count in counts.values() if count > 1),
    }


def inspect_job(job: Path) -> tuple[dict[str, np.ndarray], dict]:
    indices, audit = {}, {}
    for source, (filename, event_branch) in SOURCES.items():
        with uproot.open(job / filename) as root_file:
            indices[source], audit[source] = event_indices(root_file["Events"], event_branch)
    common = indices["analysisLevel"].keys() & indices["raw_sdst"].keys()
    # Native source order is retained in both output trees.
    keys = [key for key in indices["analysisLevel"] if key in common]
    paired = {source: np.asarray([lookup[key] for key in keys], dtype=np.int64)
              for source, lookup in indices.items()}
    audit["paired_events"] = len(keys)
    for source in SOURCES:
        audit[source]["unmatched_events"] = len(indices[source]) - len(keys)
    return paired, audit


def discover_samples(input_root: Path, samples: list[str] | None = None) -> list[Path]:
    requested = set(samples or ())
    paths = [path for path in sorted(input_root.glob("*photosFSR")) if path.is_dir()]
    if requested:
        paths = [path for path in paths if sample_name(path) in requested or path.name in requested
                 or f"{sample_name(path)}_photosFSR" in requested]
    if not paths:
        raise FileNotFoundError(f"No requested photosFSR sample directories in {input_root}")
    return paths


def prepare_sample(sample_dir: Path, output_root: Path, counts: tuple[int, int, int],
                   seed: int) -> dict:
    """Read one source job at a time and write the same split to both trees."""
    sample = f"{sample_name(sample_dir)}_photosFSR"
    output_dir = output_root / sample
    jobs, audit, excluded = [], {}, []
    has_production = (sample_dir / "production").is_dir()
    for job in sorted((sample_dir / "final_root").glob("job_*"), key=lambda path: int(path.name[4:])):
        if has_production and not (sample_dir / "production" / job.name / ".done").exists():
            excluded.append({"job": job.name, "reason": "production incomplete"})
            continue
        if any(not (job / filename).exists() for filename, _ in SOURCES.values()):
            excluded.append({"job": job.name, "reason": "source file missing"})
            continue
        paired, audit[job.name] = inspect_job(job)
        jobs.append((job, paired))

    available = sum(len(paired["analysisLevel"]) for _, paired in jobs)
    total = sum(counts)
    if available < total:
        raise ValueError(f"{sample}: requested {total} events, only {available} valid pairs available")
    if any((output_dir / f"{split}.root").exists() for split in SPLITS):
        raise FileExistsError(f"Dataset already exists under {output_dir}; choose another output directory")

    # One permutation controls native and raw events together. Events are never reused.
    order = np.random.default_rng(seed).permutation(available)
    assignment = np.full(available, -1, dtype=np.int8)
    start = 0
    for split_id, count in enumerate(counts):
        assignment[order[start:start + count]] = split_id
        start += count
    output_dir.mkdir(parents=True, exist_ok=True)
    written = {split: 0 for split in SPLITS}
    split_jobs = {split: [] for split in SPLITS}
    with ExitStack() as stack:
        files = {split: stack.enter_context(uproot.recreate(output_dir / f"{split}.root.part"))
                 for split in SPLITS}
        trees = {}
        offset = 0
        for job, paired in jobs:
            n_events = len(paired["analysisLevel"])
            job_assignment = assignment[offset:offset + n_events]
            offset += n_events
            job_id = int(job.name[4:])
            for source, (filename, _) in SOURCES.items():
                with uproot.open(job / filename) as root_file:
                    arrays = root_file["Events"].arrays(library="ak")
                arrays = arrays[paired[source]]
                arrays = ak.with_field(arrays, np.full(n_events, job_id, dtype=np.int32), "job_id")
                for split_id, split in enumerate(SPLITS):
                    selected = arrays[job_assignment == split_id]
                    if (split, source) not in trees:
                        trees[split, source] = files[split].mkrntuple(TREES[source], selected.layout.form)
                    if len(selected):
                        trees[split, source].extend(selected)
            for split_id, split in enumerate(SPLITS):
                count = int(np.count_nonzero(job_assignment == split_id))
                written[split] += count
                if count:
                    split_jobs[split].append(job_id)
            print(f"  {sample}/{job.name}: {n_events} paired events", flush=True)

    for split, expected in zip(SPLITS, counts, strict=True):
        if written[split] != expected:
            raise RuntimeError(f"{sample}/{split}: wrote {written[split]}, expected {expected}")
        source_path = output_dir / f"{split}.root.part"
        with uproot.open(source_path) as root_file:
            native = root_file[TREES["analysisLevel"]]
            raw = root_file[TREES["raw_sdst"]]
            if native.num_entries != expected or raw.num_entries != expected:
                raise RuntimeError(f"{sample}/{split}: output entry count mismatch")
            for native_key, raw_key in (("job_id", "job_id"), ("Event_runNumber", "Event_runNumber"),
                                        ("Event_evtNumber", "Event_eventNumber")):
                if not np.array_equal(native[native_key].array(library="np"), raw[raw_key].array(library="np")):
                    raise RuntimeError(f"{sample}/{split}: event alignment mismatch")
        source_path.replace(output_dir / f"{split}.root")
    print(f"saved {sample}: {dict(zip(SPLITS, counts, strict=True))}", flush=True)
    return {
        "class": sample_name(sample_dir), "source_directory": str(sample_dir.resolve()),
        "available_events": available, "unused_events": available - total,
        "splits": {split: {"file": f"{sample}/{split}.root", "n_events": written[split],
                           "jobs": split_jobs[split]} for split in SPLITS},
        "audit": audit, "excluded_jobs": excluded,
    }


def prepare_dataset(input_root: Path, output_root: Path, *, counts=(60000, 15000, 15000),
                    seed: int = 1907, samples: list[str] | None = None) -> dict:
    input_root, output_root = Path(input_root), Path(output_root)
    if len(counts) != 3 or any(count <= 0 for count in counts):
        raise ValueError("Three positive counts are required: train, validation, test")
    output_root.mkdir(parents=True, exist_ok=True)
    metadata_path = output_root / "metadata.json"
    metadata = {
        "format_version": 6, "campaign": input_root.name, "input_root": str(input_root.resolve()),
        "trees": TREES, "counts": dict(zip(SPLITS, counts, strict=True)),
        "seed": seed, "split_unit": "event", "samples": {},
    }
    if metadata_path.exists():
        existing = json.loads(metadata_path.read_text())
        for key in ("format_version", "campaign", "input_root", "trees", "seed", "split_unit"):
            if existing.get(key) != metadata[key]:
                raise ValueError(f"Existing metadata differs in {key}: {metadata_path}")
        metadata["samples"] = existing["samples"]
        metadata["counts"] = existing["counts"]
    for sample_dir in discover_samples(input_root, samples):
        name = f"{sample_name(sample_dir)}_photosFSR"
        metadata["samples"][name] = prepare_sample(sample_dir, output_root, tuple(counts), seed)
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, default=DATA_ROOT / "20261015_jongwon")
    parser.add_argument("-o", "--output", type=Path, help="default: output/dataset/<campaign>")
    parser.add_argument("--counts", nargs=3, type=int, default=(60000, 15000, 15000),
                        metavar=("TRAIN", "VAL", "TEST"))
    parser.add_argument("--samples", nargs="+", help="photosFSR sample classes or directory names")
    parser.add_argument("--seed", type=int, default=1907)
    args = parser.parse_args()
    prepare_dataset(args.input, args.output or DATASET_ROOT / args.input.name,
                    counts=args.counts, seed=args.seed, samples=args.samples)


if __name__ == "__main__":
    main()
