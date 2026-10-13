"""Read each collection once and run the selected numbered plotting stages."""
from importlib import import_module

from common.plot_utils import sample_files, set_style
from .data import SOURCES, read_sample
from .results import STAGES


def run(input_root, output_root, stage=None):
    set_style()
    comparisons = {}
    for source, (filename, collection) in SOURCES.items():
        if stage == STAGES[3] and source == 'nanoaod_raw_sdst':
            continue  # The raw SDST source has no Sim truth links.
        for dataset, files in sample_files(input_root, filename):
            if input_root.name == '20260828_florian' and not dataset.endswith('_photosFSR'):
                continue  # Use exactly the five PHOTOS samples from this study.
            values = read_sample(files, source)
            for name in STAGES:
                if name == STAGES[3] and source == 'nanoaod_raw_sdst':
                    continue
                if stage in (None, name):
                    import_module(f'{__package__}.{name}').plot(output_root, dataset, source, collection, values)
            if stage in (None, STAGES[0]):
                if source == 'nanoaod_raw_sdst':
                    # Keep only what the cross-source comparison needs.
                    comparisons[dataset] = {key: values[key] for key in
                                            ('counts', 'event_records', 'duplicate_keys', 'reco')}
                elif dataset in comparisons:
                    import_module(f'{__package__}.00_data_validation').source_comparison(
                        output_root, dataset, comparisons.pop(dataset), values)
            print(f'Track {source}/{collection}/{dataset}: {values["counts"]}', flush=True)


def main(stage=None):
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, required=True)
    parser.add_argument('-o', '--output', type=Path, default=Path(__file__).resolve().parents[2] / 'plots')
    args = parser.parse_args()
    run(args.input, args.output, stage)
