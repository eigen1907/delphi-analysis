"""Read both NanoAOD sources and execute the explicitly numbered photon stages."""
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import numpy as np

from common.paths import DATASET_ROOT, PROJECT_ROOT
from common.dataset import dataset_files
from common.plot_utils import default_plot_root, set_style
from . import angular_matching as angular
from . import data_validation as validation
from . import sample_distribution as distributions
from . import truth_link_matching as truth
from .data import RECO_LABELS, SAMPLES, read_sample
from .raw_data import RECO_LABELS as RAW_RECO_LABELS


def _read_dataset(job):
    dataset, files, source, need_angles = job
    labels = RECO_LABELS if source == 'nanoaod' else RAW_RECO_LABELS
    sample = next(name for name in SAMPLES if f'_{name}_' in f'_{dataset}_')
    angles = {mode: angular.new_accumulator() for mode in labels}

    def scan(event):
        if need_angles:
            for mode, objects in event['reco'].items():
                angular.accumulate(angles[mode], event, objects)

    diagnostics = StringIO()
    with redirect_stdout(diagnostics):
        samples = read_sample(files, sample, scan, source)
    return dataset, samples, angles, diagnostics.getvalue()


def run(input_root: Path, output_root: Path, stage=None, split="train", ncpu=5):
    """Reco definitions share fiducial Gen populations; only NanoAOD has truth links."""
    set_style()
    datasets = list(dataset_files(input_root, split))
    need_angles = stage in (None, '02_angular_matching')
    native_gen = {}
    for source, labels in (('nanoaod', RECO_LABELS), ('nanoaod_raw_sdst', RAW_RECO_LABELS)):
        if source == 'nanoaod_raw_sdst' and stage == '03_truth_link_matching':
            continue
        angles = {mode: {} for mode in labels}
        values = {mode: {} for mode in labels}
        jobs = [(dataset, files, source, need_angles) for dataset, files in datasets]
        print(f"Reading gamma {source}: {len(jobs)} samples, workers={min(ncpu, len(jobs))}", flush=True)
        if ncpu == 1:
            loaded = list(map(_read_dataset, jobs))
        else:
            with ProcessPoolExecutor(max_workers=min(ncpu, len(jobs))) as pool:
                loaded = list(pool.map(_read_dataset, jobs))
        for dataset, samples, scanned, diagnostics in loaded:
            print(diagnostics, end='', flush=True)
            for mode, value in samples.items():
                values[mode][dataset] = value
                angles[mode][dataset] = scanned[mode]
            if stage in (None, '00_data_validation'):
                representative = samples['Photon_PhotonConv' if source == 'nanoaod' else 'EmShower_Stic']
                if source == 'nanoaod':
                    native_gen[dataset] = {key: representative[key] for key in
                        ('event_keys', 'gen_signatures', 'events', 'metadata_excluded', 'reco_distributions', 'reco_label')}
                else:
                    validation.source_comparison(output_root, dataset, native_gen[dataset], representative)
        bins = dict(energy=np.r_[distributions.ENERGY_EDGES, 60],
                    cos_theta=distributions.COS_THETA_EDGES)
        for mode in labels:
            cut = angular.MATCHING_CUT_DEG if need_angles else None
            for dataset, value in values[mode].items():
                if stage in (None, '00_data_validation'):
                    validation.data_validation(output_root, dataset, value)
                if stage in (None, '01_sample_distribution'):
                    distributions.sample_distributions(output_root, dataset, value, bins)
                if stage in (None, '02_angular_matching'):
                    angular.plot(output_root, dataset, value, bins, angles[mode][dataset], cut)
                if source == 'nanoaod' and stage in (None, '03_truth_link_matching'):
                    truth.plot(output_root, dataset, value, bins)
                print(f'plots: {source}/{mode} {dataset}', flush=True)
            if need_angles:
                print(f'{source}/{mode}: opening-angle reference {cut:g} deg', flush=True)


def main(stage=None):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path,
                        default=DATASET_ROOT / '20261015_jongwon',
                        help='paired dataset campaign directory or one sample directory')
    parser.add_argument('-o', '--output', type=Path, default=None)
    parser.add_argument(
        "--stage", choices=("all", "00_data_validation", "01_sample_distribution",
                            "02_angular_matching", "03_truth_link_matching"),
        default=stage or "all", help="plotting stage (default: all)",
    )
    parser.add_argument("--split", choices=("train", "val", "test"), default="train")
    parser.add_argument("--ncpu", type=int, default=5, help="parallel sample readers (default: 5)")
    args = parser.parse_args()
    if args.ncpu < 1:
        parser.error("--ncpu must be positive")
    output = args.output or default_plot_root(PROJECT_ROOT, args.input, args.split)
    run(args.input, output, None if args.stage == "all" else args.stage, args.split, args.ncpu)


if __name__ == "__main__":
    main()
