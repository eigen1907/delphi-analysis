"""Fiducial Gen efficiencies and event energy ratios in three energy bands."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from common.plot_utils import FIGURE_SIZE, finish, plot_efficiency, plot_path
from .data import ANALYSIS_SELECTIONS, ENERGY_LABELS, ENERGY_RANGES, GEN_LABELS
from .sample_distribution import ENERGY_EDGES

XLABELS = dict(energy=r'$E_\gamma^{\rm gen}$ [GeV]', cos_theta=r'$\cos\theta^{\rm gen}$')


def selection_note(ax, range_name):
    text = 'Gen fiducial'
    if ENERGY_LABELS[range_name]:
        text += '\n' + ENERGY_LABELS[range_name]
    ax.text(0.02, 0.89, text, transform=ax.transAxes, va='top', fontsize=20)


def event_energy_ratio(output, sample, values, stats, range_index):
    """Sum only Gen photons in this band; angular Reco sums use cone unions."""
    range_name = ENERGY_RANGES[range_index][0]
    energy = values['gamma']['energy']
    event_index = values['gen_event_index']
    bins = np.linspace(0, 3, 31)
    centers = (bins[:-1] + bins[1:]) / 2
    finite = np.isfinite(stats['energy_ratio'])
    histograms = []
    for index, (population, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
        selected = values['energy_selections'][range_name][population]
        gen_energy = np.bincount(event_index[selected], weights=energy[selected], minlength=values['events'])
        if 'event_energy' in stats:
            reco_energy = stats['event_energy'][:, range_index, index]
        else:
            usable = selected & finite
            reco_energy = np.bincount(event_index[usable],
                                      weights=energy[usable] * stats['energy_ratio'][usable],
                                      minlength=values['events']).astype(float, copy=False)
            ambiguous = np.bincount(event_index[selected & ~finite], minlength=values['events']) > 0
            reco_energy[ambiguous] = np.nan
        defined = (gen_energy > 0) & np.isfinite(reco_energy)
        ratio = reco_energy[defined] / gen_energy[defined]
        counts = np.histogram(ratio[(ratio >= 0) & (ratio <= 3)], bins=bins)[0]
        tail = np.count_nonzero(ratio > 3)
        assert counts.sum() + tail == len(ratio)
        mean_entries = ratio[(ratio > 0) & (ratio < 3)]
        histograms.append((counts, label, mean_entries.mean() if len(mean_entries) else np.nan,
                           tail / len(ratio) if len(ratio) else np.nan,
                           np.count_nonzero(ratio == 0) / len(ratio) if len(ratio) else np.nan))
        print(f'{sample} {values["collection"]} {range_name} {population}: '
              f'event ratios={defined.sum()}, no selected Gen={(gen_energy == 0).sum()}, '
              f'undefined energy={((gen_energy > 0) & ~np.isfinite(reco_energy)).sum()}, '
              f'zero ratio={np.count_nonzero(ratio == 0)}, >3={tail}', flush=True)
    maximum = max(1, max(counts.max() for counts, _, _, _, _ in histograms))
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for index, (counts, label, _, _, _) in enumerate(histograms):
        ax.stairs(counts, bins, label=label, color=f'C{index}')
        shown = counts > 0
        ax.errorbar(centers[shown], counts[shown], yerr=np.sqrt(counts[shown]),
                    fmt='none', color=f'C{index}', capsize=2)
    # Empty handles keep the statistics aligned with the three Gen labels.
    anchors = (0.51, 0.74, 0.87, 0.99)
    legend_style = dict(loc='upper right', frameon=False, fontsize=18,
                        title_fontsize=18, labelspacing=0.2, handleheight=1.6)
    groups = ax.legend(title='Gen', bbox_to_anchor=(anchors[0], 0.98), **legend_style)
    ax.add_artist(groups)
    empty = [Line2D([], [], color='none') for _ in histograms]
    means = ax.legend(empty, [f'{mean:.2f}' for _, _, mean, _, _ in histograms],
                      title=r'Mean $(0<R<3)$', bbox_to_anchor=(anchors[1], 0.98),
                      handlelength=0, handletextpad=0, **legend_style)
    ax.add_artist(means)
    zeros = ax.legend(empty, [f'{100 * fraction:.2f}' for _, _, _, _, fraction in histograms],
                      title='% (=0)', bbox_to_anchor=(anchors[2], 0.98),
                      handlelength=0, handletextpad=0, **legend_style)
    ax.add_artist(zeros)
    ax.legend(empty, [f'{100 * fraction:.2f}' for _, _, _, fraction, _ in histograms],
              title='% (>3)', bbox_to_anchor=(anchors[-1], 0.98),
              handlelength=0, handletextpad=0, **legend_style)
    selection_note(ax, range_name)
    ax.set_xlim(0, 3)
    ax.set_xticks(np.arange(0, 3.5, 0.5))
    ax.set_yscale('log')
    ax.set_ylim(0.5, 0.5 * (maximum / 0.5) ** 1.65)
    suffix = '' if range_name == 'all' else f'_{range_name}'
    finish(fig, ax, output / f'event_energy_ratio{suffix}_log.png', sample,
           r'$\sum E^{\rm reco}/\sum E_\gamma^{\rm gen}$', 'Events')


def efficiency(output, sample, values, stats, bins, coordinate, range_name):
    masks = values['energy_selections'][range_name]
    edges = ENERGY_EDGES if coordinate == 'energy' else bins['cos_theta']
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for index, (population, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
        selected = masks[population]
        plot_efficiency(ax, values['gamma'][coordinate][selected],
                        stats['reco_count'][selected] > 0, edges, label, f'C{index}')
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, 1.35)
    ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
    ax.set_yticks(np.linspace(0, 1, 6))
    selection_note(ax, range_name)
    suffix = '' if range_name == 'all' else f'_{range_name}'
    ylabel = 'Angular EM coverage' if values['source'] == 'nanoaod_raw_sdst' else 'Efficiency'
    finish(fig, ax, output / f'eff_vs_{coordinate}{suffix}.png', sample, XLABELS[coordinate], ylabel)


def matching_results(study, sample, values, bins, stats, stage):
    output = plot_path(study, '02_gamma', stage, sample, values['source'],
                       values['collection'], 'matching_result')
    efficiency(output, sample, values, stats, bins, 'energy', 'all')
    for range_index, (range_name, _, _) in enumerate(ENERGY_RANGES):
        efficiency(output, sample, values, stats, bins, 'cos_theta', range_name)
        event_energy_ratio(output, sample, values, stats, range_index)
