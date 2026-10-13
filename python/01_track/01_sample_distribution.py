"""Compare Gen and Reco track distributions, with labelled momentum overflow."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from common.plot_utils import FIGURE_SIZE, finish, plot_path
from .results import STAGES, LABELS, histogram, momentum_bins


GEN_SPECIES = ((11, r'$e^\pm$'), (13, r'$\mu^\pm$'), (211, r'$\pi^\pm$'),
               (321, r'$K^\pm$'), (2212, r'$p/\bar p$'))


def gen_spectrum(ax, numbers, pdg, bins, events):
    """Disjoint charged Gen species form one color-only stack."""
    bottom = np.zeros(len(bins) - 1)
    for index, (code, label) in enumerate(GEN_SPECIES):
        counts = np.histogram(numbers[abs(pdg) == code], bins=bins)[0]
        ax.stairs(bottom + counts / events, bins, baseline=bottom, fill=True,
                  color=f'C{index}', label=f'Gen {label}')
        bottom += counts / events
    total = np.histogram(numbers, bins=bins)[0]
    assert total.sum() == len(numbers) and np.allclose(bottom * events, total)
    centers, shown = (bins[:-1] + bins[1:]) / 2, total > 0
    ax.errorbar(centers[shown], bottom[shown], yerr=np.sqrt(total[shown]) / events,
                fmt='none', color='grey', capsize=2)
    return bottom


def reco_spectrum(ax, numbers, bins, events, label='Reco'):
    counts = np.histogram(numbers, bins=bins)[0]
    assert counts.sum() == len(numbers)
    ax.stairs(counts / events, bins, fill=True, facecolor='none', edgecolor='black',
              hatch='///', linewidth=2.1, label=label)
    centers, shown = (bins[:-1] + bins[1:]) / 2, counts > 0
    ax.errorbar(centers[shown], counts[shown] / events, yerr=np.sqrt(counts[shown]) / events,
                fmt='none', color='black', capsize=2)
    return counts / events


def plot(root, dataset, source, collection, values):
    output = plot_path(root, '01_track', STAGES[1], dataset, source, collection)
    for name in LABELS:
        entries = np.r_[values['gen'][name], values['reco'][name]]
        bins = np.r_[np.linspace(0, 60, 31), 70] if name in ('p', 'pt') else momentum_bins(entries, name)
        for log in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            gen, reco = values['gen'][name], values['reco'][name]
            if name in ('p', 'pt'):
                # Cap display coordinates only; all >=60 GeV entries stay in
                # the labelled overflow bin. Matching and response use full p.
                gen, reco = np.minimum(gen, np.nextafter(70, -np.inf)), np.minimum(reco, np.nextafter(70, -np.inf))
            gen_height = gen_spectrum(ax, gen, values['gen']['pdg'], bins, values['counts']['events'])
            reco_height = reco_spectrum(ax, reco, bins, values['counts']['events'])
            height = max(gen_height.max(initial=0), reco_height.max(initial=0))
            if name in ('p', 'pt'):
                ax.set_xticks([0, 20, 40, 65], ['0', '20', '40', r'$\geq60$'])
            if log:
                ax.set_yscale('log')
                ax.set_ylim(0.5 / values['counts']['events'], max(height, 1 / values['counts']['events']) * 100)
            else:
                ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
            suffix = '_log' if log else ''
            finish(fig, ax, output / f'{name}{suffix}.png', dataset, LABELS[name], 'Tracks / bin per event')
    entries = np.r_[values['gen_counts'], values['reco_counts']]
    bins = np.arange(-0.5, entries.max(initial=0) + 1.5)
    for log in (False, True):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for index, level in enumerate(('gen', 'reco')):
            histogram(ax, values[f'{level}_counts'], bins, level.capitalize(), f'C{index}')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if log:
            ax.set_yscale('log')
            ax.set_ylim(bottom=0.5)
        else:
            ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
        suffix = '_log' if log else ''
        finish(fig, ax, output / f'multiplicity{suffix}.png', dataset, r'$N_{\rm track}$', 'Events')
