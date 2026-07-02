#!/usr/bin/env python3
from __future__ import annotations

import os
import argparse
import shutil
import subprocess
from pathlib import Path

import uproot


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))

from plot_utils import add_samples_argument

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


def write_tree(paths: list[Path], tree_name: str, output_file) -> None:
    output_tree = None

    for path in paths:
        with uproot.open(path) as root_file:
            if tree_name not in root_file:
                continue

            tree = root_file[tree_name]
            if tree.num_entries == 0:
                batches = [tree.arrays(library="ak")]
            else:
                batches = tree.iterate(step_size="100 MB", library="ak")

            for arrays in batches:
                if output_tree is None:
                    output_file[tree_name] = arrays
                    output_tree = output_file[tree_name]
                else:
                    output_tree.extend(arrays)


def sum_histograms(paths: list[Path], name: str):
    total = None
    for path in paths:
        with uproot.open(path) as root_file:
            if name not in root_file or not hasattr(root_file[name], "to_hist"):
                continue

            hist = root_file[name].to_hist()
            total = hist if total is None else total + hist
    return total


def build_file_with_hadd(paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    hadd = shutil.which("hadd")
    if hadd is None:
        raise RuntimeError("ROOT hadd was not found in PATH")

    subprocess.run(
        [hadd, "-fk", "-v", "0", str(output_path), *(str(path) for path in paths)],
        check=True,
    )


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
                write_tree(paths, name, output_file)
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
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", required=True, type=Path, help="input chunks root")
    parser.add_argument("-o", "--output", required=True, type=Path, help="merged dataset output root")
    parser.add_argument(
        "--backend",
        choices=("uproot", "hadd"),
        default="uproot",
        help="merge backend (default: uproot)",
    )
    args = parser.parse_args()

    input_root = args.input
    output_root = args.output
    requested_samples = set(args.samples or ())

    for sample_dir in sorted(path for path in input_root.iterdir() if path.is_dir()):
        sample = sample_name(sample_dir)
        if requested_samples and sample not in requested_samples and sample_dir.name not in requested_samples:
            continue
        job_dirs = sorted((sample_dir / "final_root").glob("job_*"), key=job_number)

        for file_name in ROOT_FILES:
            paths = [job_dir / file_name for job_dir in job_dirs if (job_dir / file_name).exists()]
            if not paths:
                continue

            output_path = output_root / sample / file_name
            if args.backend == "hadd":
                build_file_with_hadd(paths, output_path)
            else:
                build_file(paths, output_path)
            print(f"{sample}/{file_name}: {len(paths)} chunks -> {output_path}")


if __name__ == "__main__":
    main()
