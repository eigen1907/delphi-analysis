"""Read each collection once and run the selected numbered plotting stages."""
from common.paths import DATASET_ROOT, PROJECT_ROOT
from common.dataset import dataset_files
from common.plot_utils import add_samples_argument, default_plot_root, set_style
from . import angular_matching, data_validation, sample_distribution, truth_link_matching
from .data import SOURCES, read_sample
from .results import STAGES

STAGE_MODULES = {
    STAGES[0]: data_validation,
    STAGES[1]: sample_distribution,
    STAGES[2]: angular_matching,
    STAGES[3]: truth_link_matching,
}


def run(input_root, output_root, stage=None, split="train", samples=None):
    set_style()
    comparisons = {}
    for source, collection in SOURCES.items():
        if stage == STAGES[3] and source == 'nanoaod_raw_sdst':
            continue  # The raw SDST source has no Sim truth links.
        for dataset, files in dataset_files(input_root, split, samples):
            values = read_sample(files, source)
            for name, module in STAGE_MODULES.items():
                if name == STAGES[3] and source == 'nanoaod_raw_sdst':
                    continue
                if stage in (None, name):
                    module.plot(output_root, dataset, source, collection, values)
            if stage in (None, STAGES[0]):
                if source == 'nanoaod_raw_sdst':
                    # Keep only what the cross-source comparison needs.
                    comparisons[dataset] = {key: values[key] for key in
                                            ('counts', 'event_records', 'duplicate_keys', 'reco')}
                elif dataset in comparisons:
                    data_validation.source_comparison(
                        output_root, dataset, comparisons.pop(dataset), values)
            print(f'Track {source}/{collection}/{dataset}: {values["counts"]}', flush=True)


def main(stage=None):
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument('-i', '--input', type=Path,
                        default=DATASET_ROOT / '20260828_florian')
    parser.add_argument('-o', '--output', type=Path, default=None)
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
