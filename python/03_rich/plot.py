"""Read each RICH sample once and run the numbered analysis stages."""
import argparse
from importlib import import_module
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from common.plot_utils import sample_files, set_style
from .data import SOURCES, read_sample

PROJECT_ROOT = Path(__file__).resolve().parents[2]
STAGES = ("00_data_validation", "01_sample_distribution", "02_angular_matching", "03_truth_link_matching")


def run(input_root, output_root, stage=None):
    set_style()
    print("RICH outputs measure saved-angle coverage, not valid RICH/PID response. "
          "Gas/liquid ABI and sigma are defective; no quality bit definition is inferred.")
    modules = [import_module("." + name, __package__) for name in STAGES if stage in (None, name)]
    for source in SOURCES:
        if stage == "03_truth_link_matching" and source == "nanoaod_raw_sdst":
            continue
        for sample, files in sample_files(input_root, f"{source}.root"):
            if sample.endswith("_pythiaFSR"):
                continue
            print(f"Reading RICH {sample} {source}", flush=True)
            values = read_sample(files, source)
            for module in modules:
                module.run(values, files, source, output_root, sample)


def main(stage=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, default=PROJECT_ROOT / "plots")
    args = parser.parse_args()
    run(args.input, args.output, stage)


if __name__ == "__main__":
    main()
