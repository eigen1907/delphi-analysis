#!/usr/bin/env python3
from __future__ import annotations

import argparse
from importlib import import_module
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

from common.plot_utils import add_samples_argument, default_plot_root


def main() -> None:
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", type=Path,
                        default=PROJECT_ROOT / "data/202606xx_jongwon/dataset/20260606_100kTest",
                        help="input dataset root")
    parser.add_argument("-o", "--output", type=Path, help="plot output root (default: plots/<input directory>)")
    args = parser.parse_args()

    plot_branches_all = import_module("00_data_validation.00_branch_distribution").plot_branches_all

    input_root = args.input
    output_root = args.output or default_plot_root(PROJECT_ROOT, input_root)
    plot_branches_all(input_root, output_root, args.samples)


if __name__ == "__main__":
    main()
