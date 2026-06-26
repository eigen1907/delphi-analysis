#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=Path, default=PROJECT_ROOT / "output" / "dataset" / "100kTest")
    parser.add_argument("-o", "--output", type=Path, default=PROJECT_ROOT / "plots" / "100kTest")
    args = parser.parse_args()

    from plot_gen_check import plot_gen_check

    plot_gen_check(args.input, args.output)


if __name__ == "__main__":
    main()
