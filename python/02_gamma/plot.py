"""Read both NanoAOD sources and execute the explicitly numbered photon stages."""
from importlib import import_module
from pathlib import Path

import numpy as np

from common.plot_utils import set_style
from .data import RECO_LABELS, SAMPLES, read_sample
from .raw_data import RECO_LABELS as RAW_RECO_LABELS

validation = import_module('.00_data_validation', __package__)
distributions = import_module('.01_sample_distribution', __package__)
angular = import_module('.02_angular_matching', __package__)
truth = import_module('.03_truth_link_matching', __package__)


def run(input_root: Path, output_root: Path, stage=None):
    """Reco definitions share fiducial Gen populations; only NanoAOD has truth links."""
    set_style()
    datasets = ([input_root] if (input_root / 'final_root').is_dir()
                else sorted(input_root.glob('*_photosFSR*')))
    datasets = [path for path in datasets if (path / 'final_root').is_dir()]
    if not datasets:
        raise FileNotFoundError(f'No PHOTOS-FSR datasets in {input_root}')
    need_angles = stage in (None, '02_angular_matching')
    native_reference = 3.0  # Previous five-sample truth-angle scan: 99% ISR retention.
    native_gen = {}
    for source, labels in (('nanoaod', RECO_LABELS), ('nanoaod_raw_sdst', RAW_RECO_LABELS)):
        if source == 'nanoaod_raw_sdst' and stage == '03_truth_link_matching':
            continue
        angles = {mode: {} for mode in labels}
        values = {mode: {} for mode in labels}
        for directory in datasets:
            dataset = directory.name
            sample = next(name for name in SAMPLES if f'_{name}_' in f'_{dataset}_')
            for mode in labels:
                angles[mode][dataset] = angular.new_accumulator()
            def scan(event):
                if need_angles:
                    for mode, objects in event['reco'].items():
                        angular.accumulate(angles[mode][dataset], event, objects)
            samples = read_sample(directory, sample, scan, source)
            for mode, value in samples.items():
                values[mode][dataset] = value
            if stage in (None, '00_data_validation'):
                representative = samples['Photon_PhotonConv' if source == 'nanoaod' else 'EmShower_Stic']
                if source == 'nanoaod':
                    native_gen[dataset] = {key: representative[key] for key in
                        ('event_keys', 'gen_signatures', 'events', 'metadata_excluded', 'reco_distributions', 'reco_label')}
                else:
                    validation.source_comparison(output_root, dataset, native_gen[dataset], representative)
        all_energy = np.concatenate([np.r_[value['gamma']['energy'], value['reco_distributions']['energy']]
                                     for samples in values.values() for value in samples.values()])
        bins = dict(energy=np.r_[distributions.ENERGY_EDGES, 60],
                    energy_log=distributions.energy_bins(all_energy, True),
                    cos_theta=distributions.COS_THETA_EDGES)
        for mode in labels:
            cut = (angular.choose_cut(angles[mode].values()) if source == 'nanoaod' else native_reference) if need_angles else None
            if mode == 'Photon_PhotonConv' and need_angles:
                native_reference = cut
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
                        default=Path(__file__).resolve().parents[2] / 'data/20260828_florian',
                        help='PHOTOS-FSR campaign directory or one dataset directory')
    parser.add_argument('-o', '--output', type=Path, default=Path(__file__).resolve().parents[2] / 'plots')
    args = parser.parse_args()
    run(args.input, args.output, stage)
