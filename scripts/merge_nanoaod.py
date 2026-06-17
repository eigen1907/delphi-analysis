#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from merge_nanoaod import merge_nanoaods
from utils import job_dirs, project_root, sample_label


PROJECT_ROOT = project_root()
SAMPLE_SET = "100kTest"
INPUT_ROOT = PROJECT_ROOT / "output" / "samples" / SAMPLE_SET
OUTPUT_ROOT = PROJECT_ROOT / "output" / "merged" / SAMPLE_SET


def main() -> None:
    for sample_dir in sorted(path for path in INPUT_ROOT.iterdir() if path.is_dir()):
        label = sample_label(sample_dir)
        for job_dir in job_dirs(sample_dir):
            output = OUTPUT_ROOT / label / f"{job_dir.name}.root"
            summary = merge_nanoaods(
                job_dir / "nanoaod.root",
                job_dir / "nanoaod_raw_sdst.root",
                job_dir / "nanoaod_raw_fadana.root",
                output,
            )
            print(
                f"{label}/{job_dir.name}: wrote {summary['events']} events "
                f"to {output.relative_to(PROJECT_ROOT)}"
            )


if __name__ == "__main__":
    main()
