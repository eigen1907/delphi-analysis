#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from plot_variables import plot_variables
from utils import project_root


PROJECT_ROOT = project_root()
SAMPLE_SET = "100kTest"
INPUT_SAMPLES = [
    ("ZKK", PROJECT_ROOT / "output" / "hadd" / SAMPLE_SET / "ZKK.root"),
    ("Zpipi", PROJECT_ROOT / "output" / "hadd" / SAMPLE_SET / "Zpipi.root"),
    ("Zee", PROJECT_ROOT / "output" / "hadd" / SAMPLE_SET / "Zee.root"),
    ("Zmumu", PROJECT_ROOT / "output" / "hadd" / SAMPLE_SET / "Zmumu.root"),
]
OUTPUT_ROOT = PROJECT_ROOT / "plots" / SAMPLE_SET


def main() -> None:
    result = plot_variables(INPUT_SAMPLES, OUTPUT_ROOT)
    print(result["summary"].to_string(index=False))


if __name__ == "__main__":
    main()
