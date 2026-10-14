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
    """Reco definitions share Gen cuts; only standard NanoAOD has truth links."""
    set_style()
    need_angles = stage in (None, '02_angular_matching')
    native_reference = 3.0  # Previous five-sample truth-angle scan: 99% ISR retention.
    native_gen = {}
    for source, labels in (('nanoaod', RECO_LABELS), ('nanoaod_raw_sdst', RAW_RECO_LABELS)):
        if source == 'nanoaod_raw_sdst' and stage == '03_truth_link_matching':
            continue
        angles = {mode: {} for mode in labels}
        values = {mode: {} for mode in labels}
        for sample in SAMPLES:
            for mode in labels:
                angles[mode][sample] = angular.new_accumulator()
            def scan(event):
                if need_angles:
                    for mode, objects in event['reco'].items():
                        angular.accumulate(angles[mode][sample], event, objects)
            samples = read_sample(input_root, sample, scan, source)
            for mode, value in samples.items():
                values[mode][sample] = value
            if stage in (None, '00_data_validation'):
                representative = samples['Photon_PhotonConv' if source == 'nanoaod' else 'EmShower_Stic']
                if source == 'nanoaod':
                    native_gen[sample] = {key: representative[key] for key in
                        ('event_keys', 'gen_signatures', 'events', 'metadata_excluded', 'reco_distributions', 'reco_label')}
                else:
                    dataset = f'20260828_100kTest_{sample}_photosFSR'
                    validation.source_comparison(output_root, dataset, native_gen[sample], representative)
        all_energy = np.concatenate([np.r_[value['gamma']['energy'], value['reco_distributions']['energy']]
                                     for samples in values.values() for value in samples.values()])
        bins = dict(energy=np.r_[distributions.ENERGY_EDGES, 60],
                    energy_log=distributions.energy_bins(all_energy, True),
                    cos_theta=distributions.COS_THETA_EDGES)
        gen_energy = np.concatenate([value['gamma']['energy'] for value in values[next(iter(labels))].values()])
        gen_bins = dict(bins, energy=distributions.energy_bins(gen_energy))
        for mode in labels:
            cut = (angular.choose_cut(angles[mode].values()) if source == 'nanoaod' else native_reference) if need_angles else None
            if mode == 'Photon_PhotonConv' and need_angles:
                native_reference = cut
            for sample, value in values[mode].items():
                dataset = f'20260828_100kTest_{sample}_photosFSR'
                if stage in (None, '00_data_validation'):
                    validation.data_validation(output_root, dataset, value)
                if stage in (None, '01_sample_distribution'):
                    distributions.sample_distributions(output_root, dataset, value, bins)
                if stage in (None, '02_angular_matching'):
                    angular.plot(output_root, dataset, value, gen_bins, angles[mode][sample], cut)
                if source == 'nanoaod' and stage in (None, '03_truth_link_matching'):
                    truth.plot(output_root, dataset, value, gen_bins)
                print(f'plots: {source}/{mode} {dataset}', flush=True)
            if need_angles:
                print(f'{source}/{mode}: opening-angle reference {cut:g} deg', flush=True)


def main(stage=None):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path,
                        default=Path(__file__).resolve().parents[2] / 'data/20260828_florian')
    parser.add_argument('-o', '--output', type=Path, default=Path(__file__).resolve().parents[2] / 'plots')
    args = parser.parse_args()
    run(args.input, args.output, stage)
