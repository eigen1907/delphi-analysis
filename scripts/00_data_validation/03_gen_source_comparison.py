#!/usr/bin/env python3
from __future__ import annotations

import argparse
from importlib import import_module
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from common.plot_utils import add_samples_argument, default_check_root, default_plot_root


def main() -> None:
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", type=Path,
                        default=PROJECT_ROOT / "data/202606xx_jongwon/dataset/20260606_100kTest",
                        help="input dataset root")
    parser.add_argument("-o", "--output", type=Path, help="plot output root (default: plots/<input directory>)")
    parser.add_argument("--data-root", type=Path, help="diagnostic data root (default: data/check/<input directory>)")
    args = parser.parse_args()

    plot_gen_compare = import_module("00_data_validation.03_gen_source_comparison").plot_gen_compare

    input_root = args.input
    output_root = args.output or default_plot_root(PROJECT_ROOT, input_root)
    data_root = args.data_root or default_check_root(PROJECT_ROOT, input_root)
    plot_gen_compare(input_root, output_root, data_root, args.samples)


if __name__ == "__main__":
    main()
