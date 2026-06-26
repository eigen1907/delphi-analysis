#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import uproot


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))

from plot_utils import FILE_GROUPS, SAMPLE_LABELS, TREE_NAMES_BY_SOURCE, resolve_file


def relative_to_project(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=Path, default=PROJECT_ROOT / "output" / "dataset" / "100kTest")
    parser.add_argument("-o", "--output", type=Path, default=PROJECT_ROOT / "data" / "branch")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    for sample in SAMPLE_LABELS:
        summary = {"sample": sample, "trees": {}}
        sample_dir = args.input / sample

        for source, file_names in FILE_GROUPS:
            path = resolve_file(sample_dir, file_names)

            with uproot.open(path) as root_file:
                for tree_name in TREE_NAMES_BY_SOURCE.get(source, ("Events",)):
                    if tree_name not in root_file:
                        continue

                    tree = root_file[tree_name]
                    tree_summary = summary["trees"].setdefault(tree_name, {"branches": {}})
                    for branch in sorted(tree.typenames()):
                        branch_summary = tree_summary["branches"].setdefault(branch, {"source": []})
                        branch_summary["source"].append(source)

        output = args.output / f"{sample}.json"
        output.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"[{sample}] wrote {relative_to_project(output)}")


if __name__ == "__main__":
    main()
