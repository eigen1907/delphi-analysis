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
    parser.add_argument("--alpha-cut", type=float, default=0.05)
    parser.add_argument("--fiducial-cos-max", type=float, default=1.00)
    parser.add_argument("--pt-min", type=float, default=0.00)
    args = parser.parse_args()

    from plot_gen_reco_track_match_cut import plot_gen_reco_track_match_cut

    plot_gen_reco_track_match_cut(args.input, args.output, args.alpha_cut, args.fiducial_cos_max, args.pt_min)


if __name__ == "__main__":
    main()
