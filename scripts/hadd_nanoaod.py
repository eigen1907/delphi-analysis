#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from hadd_nanoaod import hadd_nanoaods
from utils import project_root


PROJECT_ROOT = project_root()
SAMPLE_SET = "100kTest"
INPUT_ROOT = PROJECT_ROOT / "output" / "merged" / SAMPLE_SET
OUTPUT_ROOT = PROJECT_ROOT / "output" / "hadd" / SAMPLE_SET


def main() -> None:
    for sample_dir in sorted(path for path in INPUT_ROOT.iterdir() if path.is_dir()):
        inputs = sorted(sample_dir.glob("job_*.root"))
        output = OUTPUT_ROOT / f"{sample_dir.name}.root"
        summary = hadd_nanoaods(inputs, output)
        print(
            f"{sample_dir.name}: hadd {summary['files']} files, "
            f"{summary['events']} events -> {output.relative_to(PROJECT_ROOT)}"
        )


if __name__ == "__main__":
    main()
