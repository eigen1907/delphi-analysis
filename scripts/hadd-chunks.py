#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import awkward as ak
import uproot


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))
SAMPLE_SET = "100kTest"
INPUT_ROOT = PROJECT_ROOT / "output" / "chunks" / SAMPLE_SET
OUTPUT_ROOT = PROJECT_ROOT / "output" / "dataset" / SAMPLE_SET

ROOT_FILES = (
    "nanoaod.root",
    "nanoaod_raw_sdst.root",
    "nanoaod_raw_fadana.root",
    "nanoaod_ttree.root",
    "pythiastandalone.root",
)


def sample_name(sample_dir: Path) -> str:
    return sample_dir.name.split("_")[-1]


def job_number(job_dir: Path) -> int:
    try:
        return int(job_dir.name.removeprefix("job_"))
    except ValueError:
        return 10**9


def is_tree(obj) -> bool:
    return hasattr(obj, "arrays") and hasattr(obj, "keys") and hasattr(obj, "num_entries")


def object_names(paths: list[Path]) -> list[str]:
    names = []
    seen = set()
    for path in paths:
        with uproot.open(path) as root_file:
            for name in root_file.keys(cycle=False):
                if name not in seen:
                    names.append(name)
                    seen.add(name)
    return names


def concat_tree(paths: list[Path], tree_name: str):
    arrays = []
    for path in paths:
        with uproot.open(path) as root_file:
            if tree_name in root_file:
                arrays.append(root_file[tree_name].arrays(library="ak"))

    if len(arrays) == 1:
        return arrays[0]
    return ak.concatenate(arrays, axis=0)


def sum_histograms(paths: list[Path], name: str):
    total = None
    for path in paths:
        with uproot.open(path) as root_file:
            if name not in root_file or not hasattr(root_file[name], "to_hist"):
                continue

            hist = root_file[name].to_hist()
            total = hist if total is None else total + hist
    return total


def build_file(paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with uproot.recreate(output_path) as output_file:
        for name in object_names(paths):
            first_is_tree = False
            for path in paths:
                with uproot.open(path) as first_file:
                    if name in first_file:
                        first_is_tree = is_tree(first_file[name])
                        break

            if first_is_tree:
                output_file[name] = concat_tree(paths, name)
                continue

            hist = sum_histograms(paths, name)
            if hist is not None:
                output_file[name] = hist
                continue

            for path in paths:
                with uproot.open(path) as root_file:
                    if name in root_file:
                        output_file[name] = root_file[name]
                        break


def main() -> None:
    for sample_dir in sorted(path for path in INPUT_ROOT.iterdir() if path.is_dir()):
        sample = sample_name(sample_dir)
        job_dirs = sorted((sample_dir / "final_root").glob("job_*"), key=job_number)

        for file_name in ROOT_FILES:
            paths = [job_dir / file_name for job_dir in job_dirs if (job_dir / file_name).exists()]
            if not paths:
                continue

            output_path = OUTPUT_ROOT / sample / file_name
            build_file(paths, output_path)
            print(f"{sample}/{file_name}: {len(paths)} chunks -> {output_path}")


if __name__ == "__main__":
    main()
