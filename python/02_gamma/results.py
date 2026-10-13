"""Shared Gen-level efficiencies and energy-ratio plots for both matchers."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from common.plot_utils import FIGURE_SIZE, finish, plot_efficiency, plot_path
from .data import ANALYSIS_SELECTIONS, GEN_CUTS, GEN_LABELS

XLABELS = dict(energy=r'$E_\gamma^{\rm gen}$ [GeV]', cos_theta=r'$\cos\theta^{\rm gen}$')


def energy_ratio_distributions(output, sample, values, stats, cut_index):
    """Gen-photon and event energy ratios, including unmatched photons at zero."""
    cut_name, min_gen_energy, _ = GEN_CUTS[cut_index]
    energy = values['gamma']['energy']
    event_index = values['gen_event_index']
    populations = dict(gen_energy_ratio=[], event_energy_ratio=[])
    for index, (population, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
        selected = values['cuts'][cut_name][population]
        matched = selected & (stats['reco_count'] > 0)
        finite = np.isfinite(stats['energy_ratio'])
        populations['gen_energy_ratio'].append((stats['energy_ratio'][selected & finite], label))
        gen_energy = np.bincount(event_index[selected], weights=energy[selected], minlength=values['events'])
        if 'event_energy' in stats:
            # Angular cones can share candidates: the event accumulator takes their union.
            reco_energy = stats['event_energy'][:, cut_index, index]
        else:
            usable = selected & finite
            reco_energy = np.bincount(event_index[usable],
                                      weights=energy[usable] * stats['energy_ratio'][usable],
                                      minlength=values['events'])
            ambiguous = np.bincount(event_index[selected & ~finite], minlength=values['events']) > 0
            reco_energy[ambiguous] = np.nan
        has_gen = gen_energy > 0
        defined = has_gen & np.isfinite(reco_energy)
        populations['event_energy_ratio'].append((reco_energy[defined] / gen_energy[defined], label))
        print(f'{sample} {output.parent.name}/{output.name} {population}: '
              f'Gen ratios={(selected & finite).sum()}, matched Gen={matched.sum()}, '
              f'zero Gen ratio={(selected & finite & (stats["energy_ratio"] == 0)).sum()}, '
              f'undefined Gen energy={(selected & ~finite).sum()}; '
              f'event ratios={defined.sum()}, no selected Gen={(~has_gen).sum()}, '
              f'undefined event energy={(has_gen & ~defined).sum()}, '
              f'zero recovered energy={(defined & (reco_energy == 0)).sum()}', flush=True)
    bins = np.linspace(0, 3, 31)
    centers = (bins[:-1] + bins[1:]) / 2
    for name, entries in populations.items():
        histograms = []
        for ratio, label in entries:
            in_range = ratio[(ratio >= 0) & (ratio <= 3)]
            counts = np.histogram(in_range, bins=bins)[0]
            tail = np.count_nonzero(ratio > 3)
            assert counts.sum() + tail == len(ratio)
            mean_entries = ratio[(ratio > 0) & (ratio < 3)]
            histograms.append((counts, label, mean_entries.mean() if len(mean_entries) else np.nan,
                               tail / len(ratio) if len(ratio) else np.nan,
                               np.count_nonzero(ratio == 0) / len(ratio) if len(ratio) else np.nan))
            print(f'{sample} {output.parent.name}/{output.name} {label} {name}: '
                  f'>3={tail}/{len(ratio)}', flush=True)
        maximum = max(1, max(counts.max() for counts, _, _, _, _ in histograms))
        for log in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for index, (counts, label, _, _, _) in enumerate(histograms):
                color = f'C{index}'
                ax.stairs(counts, bins, label=label, color=color)
                shown = counts > 0
                ax.errorbar(centers[shown], counts[shown], yerr=np.sqrt(counts[shown]),
                            fmt='none', color=color, capsize=2)
            # Statistics use separate legends with empty handles, aligned to the Gen labels.
            anchors = (0.27, 0.54, 0.77, 0.99)
            top = 0.86
            groups = ax.legend(title='Gen', loc='upper right', bbox_to_anchor=(anchors[0], top),
                               frameon=False, handleheight=1.6)
            ax.add_artist(groups)
            empty = [Line2D([], [], color='none') for _ in histograms]
            means = ax.legend(empty, [f'{mean:.2f}' for _, _, mean, _, _ in histograms],
                              title=r'Mean $(0<R<3)$', loc='upper right', bbox_to_anchor=(anchors[1], top),
                              frameon=False, handlelength=0, handletextpad=0, handleheight=1.6)
            ax.add_artist(means)
            zeros = ax.legend(empty, [f'{fraction:.4f}' for _, _, _, _, fraction in histograms],
                              title='Frac. (=0)', loc='upper right', bbox_to_anchor=(anchors[2], top),
                              frameon=False, handlelength=0, handletextpad=0, handleheight=1.6)
            ax.add_artist(zeros)
            ax.legend(empty, [f'{fraction:.4f}' for _, _, _, fraction, _ in histograms],
                      title='Frac. (>3)', loc='upper right', bbox_to_anchor=(anchors[-1], top),
                      frameon=False, handlelength=0, handletextpad=0, handleheight=1.6)
            if min_gen_energy:
                ax.text(0.02, 0.89, rf'$E_\gamma^{{\rm gen}}\geq{min_gen_energy}$ GeV',
                        transform=ax.transAxes, va='top')
            ax.set_xlim(0, 3)
            ax.set_xticks(np.arange(0, 3.5, 0.5))
            if log:
                ax.set_yscale('log')
                ax.set_ylim(0.5, 0.5 * (maximum / 0.5) ** 1.65)
            else:
                ax.set_ylim(0, maximum * 1.65)
            suffix = '_log' if log else ''
            xlabel = (r'$\sum E^{\rm reco}/E_\gamma^{\rm gen}$' if name == 'gen_energy_ratio'
                      else r'$\sum E^{\rm reco}/\sum E_\gamma^{\rm gen}$')
            ylabel = 'Gen photons' if name == 'gen_energy_ratio' else 'Events'
            finish(fig, ax, output / f'{name}_{sample}{suffix}.png', sample, xlabel, ylabel)


def matching_results(study, sample, values, bins, stats, stage):
    gamma = values['gamma']
    for cut_index, (cut_name, minimum, _) in enumerate(GEN_CUTS):
        output = plot_path(study, '02_gamma', stage, sample, values['source'], values['collection'],
                           'matching_result', cut_name)
        energy_ratio_distributions(output, sample, values, stats, cut_index)
        for coordinate in ('cos_theta', 'energy'):
            edges = (np.r_[minimum, bins[coordinate][bins[coordinate] > minimum]]
                     if coordinate == 'energy' else bins[coordinate])
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for index, (population, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
                selected = values['cuts'][cut_name][population]
                x = gamma[coordinate][selected]
                assert np.histogram(x, edges)[0].sum() == selected.sum()
                plot_efficiency(ax, x, stats['reco_count'][selected] > 0, edges, label, f'C{index}')
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0, 1.5)
            ax.set_yticks(np.linspace(0, 1, 6))
            if minimum:
                ax.text(0.02, 0.89, rf'$E_\gamma^{{\rm gen}}\geq{minimum}$ GeV',
                        transform=ax.transAxes, va='top')
            ylabel = 'Angular EM coverage' if values['source'] == 'nanoaod_raw_sdst' else 'Efficiency'
            finish(fig, ax, output / f'efficiency_vs_{coordinate}_{sample}.png', sample,
                   XLABELS[coordinate], ylabel)
