#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))

from plot_utils import add_samples_argument, default_plot_root


def main() -> None:
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", required=True, type=Path, help="input dataset root")
    parser.add_argument("-o", "--output", type=Path, help="plot output root (default: plots/<input directory>)")
    args = parser.parse_args()

    from plot_branches_all import plot_branches_all

    input_root = args.input
    output_root = args.output or default_plot_root(PROJECT_ROOT, input_root)
    plot_branches_all(input_root, output_root, args.samples)


if __name__ == "__main__":
    main()
