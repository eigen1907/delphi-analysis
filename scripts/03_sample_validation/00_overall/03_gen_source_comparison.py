#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from common.paths import DATASET_ROOT, PLOT_ROOT, PROJECT_ROOT
from common.plot_utils import add_samples_argument, default_check_root
from sample_validation.overall.gen_source_comparison import plot_gen_compare


def main() -> None:
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", type=Path,
                        default=DATASET_ROOT / "20261015_jongwon/merged",
                        help="input dataset root")
    parser.add_argument("-o", "--output", type=Path, default=PLOT_ROOT,
                        help="plot output root (default: output/plots)")
    parser.add_argument("--data-root", type=Path, help="diagnostic data root")
    args = parser.parse_args()

    input_root = args.input
    output_root = args.output
    data_root = args.data_root or default_check_root(PROJECT_ROOT, input_root)
    plot_gen_compare(input_root, output_root, data_root, args.samples)


if __name__ == "__main__":
    main()
