#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from common.paths import DATASET_ROOT, PLOT_ROOT
from common.plot_utils import add_color_scale_argument, add_samples_argument
from sample_validation.overall.gen_distribution import plot_gen_check


def main() -> None:
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", type=Path,
                        default=DATASET_ROOT / "20261015_jongwon/merged",
                        help="input dataset root")
    parser.add_argument("-o", "--output", type=Path, default=PLOT_ROOT,
                        help="plot output root (default: output/plots)")
    add_color_scale_argument(parser)
    args = parser.parse_args()

    input_root = args.input
    output_root = args.output
    plot_gen_check(input_root, output_root, args.samples, args.color_scale)


if __name__ == "__main__":
    main()
