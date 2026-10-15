"""Read each RICH sample once and run the numbered analysis stages."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from common.paths import DATA_ROOT, PLOT_ROOT
from common.plot_utils import sample_files, set_style
from . import angular_matching, data_validation, sample_distribution, truth_link_matching
from .data import SOURCES, read_sample

STAGES = {
    "00_data_validation": data_validation,
    "01_sample_distribution": sample_distribution,
    "02_angular_matching": angular_matching,
    "03_truth_link_matching": truth_link_matching,
}


def run(input_root, output_root, stage=None):
    set_style()
    print("RICH outputs measure saved-angle coverage, not PID efficiency. "
          "No quality bit definition is inferred.")
    modules = [module for name, module in STAGES.items() if stage in (None, name)]
    for source in SOURCES:
        if stage == "03_truth_link_matching" and source == "nanoaod_raw_sdst":
            continue
        for sample, files in sample_files(input_root, f"{source}.root"):
            if sample.endswith("_pythiaFSR"):
                continue
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
    parser.add_argument("-i", "--input", type=Path, default=DATA_ROOT / "20260828_florian")
    parser.add_argument("-o", "--output", type=Path, default=PLOT_ROOT)
    args = parser.parse_args()
    run(args.input, args.output, stage)


if __name__ == "__main__":
    main()
