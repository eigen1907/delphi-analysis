"""Opening-angle matching to all Parts carrying DELPHI gamma mass code 21."""

import numpy as np

from .angle_matching import accumulate, angular_stats, new_accumulator
from .data import SAMPLES, read_sample
from .plot_isr_study import detector_efficiency, energy_bins, matching_results, set_style


OPENING_ANGLE_DEG = 3


def plot_part_gamma_angular(input_root, output_root):
    set_style()
    study = output_root / 'isr_study' / 'part_gamma_angular'
    samples = {}
    for sample in SAMPLES:
        angles = new_accumulator()

        def scan(event):
            # All code-21 Parts participate, including candidates without truth links.
            accumulate(angles, event, event['reco']['gamma_pid'])

        values = read_sample(input_root, sample, scan)['gamma_pid']
        values['angular'] = angular_stats(angles, OPENING_ANGLE_DEG)
        samples[sample] = values

    bins = dict(energy=energy_bins(np.concatenate([values['gamma']['energy']
                                                  for values in samples.values()])),
                cos_theta=np.linspace(-1, 1, 21))
    for sample, values in samples.items():
        stats = values['angular']
        matching_results(study, sample, values, bins, stats, 'angular_matching')
        detector_efficiency(study, sample, values, stats, 'angular_matching')
    print(f'Part gamma angular study: {study}; opening angle <= {OPENING_ANGLE_DEG} deg', flush=True)
