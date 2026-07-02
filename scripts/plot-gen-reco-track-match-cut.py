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
    parser.add_argument("--alpha-cut", type=float, default=0.05)
    parser.add_argument("--fiducial-cos-max", type=float, default=1.00)
    parser.add_argument("--pt-min", type=float, default=0.00)
    args = parser.parse_args()

    from plot_gen_reco_track_match_cut import plot_gen_reco_track_match_cut

    input_root = args.input
    output_root = args.output or default_plot_root(PROJECT_ROOT, input_root)
    plot_gen_reco_track_match_cut(input_root, output_root, args.alpha_cut, args.fiducial_cos_max, args.pt_min, args.samples)


if __name__ == "__main__":
    main()
