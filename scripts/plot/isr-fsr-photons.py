#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from plot_utils import default_plot_root


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path, help="Florian study directory")
    parser.add_argument("-o", "--output", type=Path, help="plot output root (default: plots/<input directory>)")
    args = parser.parse_args()

    from plot_isr_fsr_photons import plot_isr_fsr_photons

    output_root = args.output or default_plot_root(PROJECT_ROOT, args.input)
    plot_isr_fsr_photons(args.input, output_root)


if __name__ == "__main__":
    main()
