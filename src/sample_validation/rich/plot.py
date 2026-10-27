"""Read each RICH sample once and run the numbered analysis stages."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from common.paths import DATASET_ROOT, PROJECT_ROOT
from common.dataset import dataset_files
from common.plot_utils import add_samples_argument, default_plot_root, set_style
from . import angular_matching, data_validation, sample_distribution, truth_link_matching
from .data import SOURCES, read_sample

STAGES = {
    "00_data_validation": data_validation,
    "01_sample_distribution": sample_distribution,
    "02_angular_matching": angular_matching,
    "03_truth_link_matching": truth_link_matching,
}


def run(input_root, output_root, stage=None, split="train", samples=None):
    set_style()
    print("RICH outputs measure saved-angle coverage, not PID efficiency. "
          "No quality bit definition is inferred.")
    modules = [module for name, module in STAGES.items() if stage in (None, name)]
    for source in SOURCES:
        if stage == "03_truth_link_matching" and source == "nanoaod_raw_sdst":
            continue
        for sample, files in dataset_files(input_root, split, samples):
            print(f"Reading RICH {sample} {source}", flush=True)
            values = read_sample(files, source)
            if not values["association_supported"]:
                print("Old RICH schema: ABI/index fixes are absent; only storage checks are usable.", flush=True)
            for module in modules:
                if not values["association_supported"] and module is not data_validation:
                    continue
                module.run(values, files, source, output_root, sample)


def main(stage=None):
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", type=Path, default=DATASET_ROOT / "20260828_florian")
    parser.add_argument("-o", "--output", type=Path, default=None)
    parser.add_argument(
        "--stage", choices=("all", "00_data_validation", "01_sample_distribution",
                            "02_angular_matching", "03_truth_link_matching"),
        default=stage or "all", help="plotting stage (default: all)",
    )
    parser.add_argument("--split", choices=("train", "val", "test"), default="train")
    args = parser.parse_args()
    output = args.output or default_plot_root(PROJECT_ROOT, args.input, args.split)
    run(args.input, output, None if args.stage == "all" else args.stage, args.split, args.samples)


if __name__ == "__main__":
    main()
