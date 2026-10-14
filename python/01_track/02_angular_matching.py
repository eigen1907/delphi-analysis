"""Validate the same-charge angular cut and plot its one-to-one matches."""
import numpy as np
import matplotlib.pyplot as plt

from common.plot_utils import FIGURE_SIZE, clopper_pearson, finish, plot_path
from .data import ANGLE_CUT, SCAN_CUTS
from .results import matching_result


def plot(root, dataset, source, collection, values):
    output = plot_path(root, '01_track', '02_angular_matching', dataset, source, collection,
                       section='matching_validation', cut='cut_no')
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    # Greedy matching at the widest scan cut gives the same sorted assignment
    # at every narrower cut. Truth information does not enter this scan.
    success = np.array([(values['scan_angle'] < cut).sum() for cut in SCAN_CUTS])
    denominator = np.full(len(SCAN_CUTS), len(values['scan_angle']))
    lower, upper = clopper_pearson(success, denominator)
    ax.plot(SCAN_CUTS, success / denominator, color='C0', label='Same charge, one-to-one')
    ax.fill_between(SCAN_CUTS, lower, upper, color='C0', alpha=0.2)
    ax.axvline(ANGLE_CUT, color='black', linestyle='--', linewidth=2)
    ax.set_ylim(0, 1.35)
    ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
    finish(fig, ax, output / 'efficiency_vs_opening_angle.png', dataset, 'Opening-angle cut [rad]', 'Efficiency')
    matching_result(root, dataset, source, collection, values, 'angular')
