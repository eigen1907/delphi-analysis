"""Plot helpers and results shared by the two track matching stages."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from common.plot_utils import FIGURE_SIZE, finish, plot_efficiency, plot_path

STAGES = ("00_data_validation", "01_sample_distribution", "02_angular_matching", "03_truth_link_matching")
LABELS = {"p": r'$p$ [GeV]', "pt": r'$p_T$ [GeV]', "cos_theta": r'$\cos\theta$', "phi": r'$\phi$ [rad]'}


def momentum_bins(values, name):
    if name == "cos_theta":
        return np.linspace(-1, 1, 101)
    if name == "phi":
        return np.linspace(-np.pi, np.pi, 73)
    # Gen matching bins retain the complete observed momentum range.
    return np.r_[np.linspace(0, 60, 31), max(62, np.nextafter(max(values.max(initial=0), 60), np.inf))]


def histogram(ax, numbers, bins, label, color, denominator=1):
    counts = np.histogram(numbers, bins=bins)[0]
    assert counts.sum() == len(numbers)
    ax.stairs(counts / denominator, bins, label=label, color=color)
    shown = counts > 0
    centers = (bins[:-1] + bins[1:]) / 2
    ax.errorbar(centers[shown], counts[shown] / denominator,
                yerr=np.sqrt(counts[shown]) / denominator, fmt='none', color=color, capsize=2)
    return counts / denominator


def residual_bins(entries, scale):
    """Fine central bins and logarithmic tails retain outliers without hiding the peak."""
    low, high = entries.min(initial=0), np.nextafter(entries.max(initial=0), np.inf)
    left = -np.geomspace(scale, -low, 21)[1:][::-1] if low < -scale else []
    right = np.geomspace(scale, high, 21)[1:] if high > scale else []
    return np.r_[left, np.linspace(-scale, scale, 61), right]


def matching_result(root, dataset, source, collection, values, method):
    stage = STAGES[2] if method == 'angular' else STAGES[3]
    output = plot_path(root, '01_track', stage, dataset, source, collection,
                       section='matching_result', cut='cut_no')
    stats = values['methods'][method]
    for name in ('p', 'pt', 'cos_theta'):
        bins = momentum_bins(values['gen'][name], name)
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        plot_efficiency(ax, values['gen'][name], stats['matched'], bins, 'Gen → Reco', 'C0')
        ax.set_ylim(0, 1.35)
        ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
        finish(fig, ax, output / f'efficiency_vs_{name}.png', dataset, LABELS[name], 'Efficiency')
    for name, xlabel in (('pt_response', r'$(p_{T,\rm reco}-p_{T,\rm gen})/p_{T,\rm gen}$'),
                         ('theta_residual', r'$\theta_{\rm reco}-\theta_{\rm gen}$ [rad]'),
                         ('phi_residual', r'$\phi_{\rm reco}-\phi_{\rm gen}$ [rad]'),
                         ('angle', 'Opening angle [rad]')):
        entries = stats[name][stats['matched']]
        scale = 1 if 'response' in name else 0.05
        bins = residual_bins(entries, scale)
        if name == 'angle':
            bins = bins[bins >= 0]
        for log in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            histogram(ax, entries, bins, 'Matched track', 'C0')
            ax.set_xscale('symlog', linthresh=scale)
            if log:
                ax.set_yscale('log')
                ax.set_ylim(bottom=0.5)
            else:
                ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
            suffix = '_log' if log else ''
            finish(fig, ax, output / f'{name}{suffix}.png', dataset, xlabel, 'Gen tracks')
    ratios = stats['momentum_ratio']
    visible = ratios[ratios <= 3]
    unmatched = np.count_nonzero(~stats['matched']) / len(ratios)
    tail = np.count_nonzero(ratios > 3) / len(ratios)
    assert len(visible) + np.count_nonzero(ratios > 3) == len(ratios)
    label = f'Gen charged\nUnmatched: {unmatched:.3f}\nRatio >3: {tail:.3f}'
    for log in (False, True):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        histogram(ax, visible, np.linspace(0, 3, 61), label, 'C0')
        ax.set_xlim(0, 3)
        if log:
            ax.set_yscale('log')
            ax.set_ylim(bottom=0.5)
        else:
            ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
        suffix = '_log' if log else ''
        finish(fig, ax, output / f'momentum_ratio{suffix}.png', dataset,
               r'$p_{\rm reco}^{\rm leading}/p_{\rm gen}$', 'Gen tracks')
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        counts = stats['reco_count']
        histogram(ax, counts, np.arange(-0.5, counts.max(initial=0) + 1.5), 'Gen charged', 'C0')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if log:
            ax.set_yscale('log')
            ax.set_ylim(bottom=0.5)
        else:
            ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
        finish(fig, ax, output / f'matched_reco_multiplicity{suffix}.png', dataset,
               r'$N_{\rm reco}$ per Gen track', 'Gen tracks')
