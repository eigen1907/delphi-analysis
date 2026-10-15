"""Gen-region efficiency and the saved Reco detector tags, kept as diagnostics."""
import matplotlib.pyplot as plt
import numpy as np

from common.plot_utils import FIGURE_SIZE, clopper_pearson, finish, plot_path
from .data import (
    ANALYSIS_SELECTIONS, DETECTORS, EM_REGIONS, ENERGY_LABELS, ENERGY_RANGES,
    GEN_LABELS, em_regions,
)


def detector_validation(root, dataset, values):
    """Inclusive Reco directions; overlapping flags are deliberately not stacked."""
    reco = values['reco_distributions']
    flags = np.asarray(reco['detector_flags'], dtype=bool).reshape(-1, len(DETECTORS))
    masks = [*flags.T, ~flags.any(axis=1)]
    bins = np.linspace(-1, 1, 101)
    centers = (bins[:-1] + bins[1:]) / 2
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    counts_by_tag = {}
    for index, (label, selected) in enumerate(zip((*DETECTORS, 'No tag'), masks, strict=True)):
        entries = reco['cos_theta'][selected]
        entries = entries[np.isfinite(entries)]
        counts = np.histogram(entries, bins)[0]
        counts_by_tag[label] = int(selected.sum())
        color = f'C{index}'
        ax.stairs(counts / values['events'], bins, label=label, color=color)
        shown = counts > 0
        ax.errorbar(centers[shown], counts[shown] / values['events'],
                    yerr=np.sqrt(counts[shown]) / values['events'],
                    fmt='none', color=color, capsize=2)
    ax.set_xlim(-1, 1)
    ax.set_yscale('log')
    ax.set_ylim(0.5 / values['events'], max(ax.get_ylim()[1] * 30, 1 / values['events']))
    output = plot_path(root, '02_gamma', '00_data_validation', dataset,
                       values['source'], values['collection']) / 'branch_checks'
    finish(fig, ax, output / 'detector_tag_cos_theta_log.png', dataset,
           r'$\cos\theta^{\rm reco}$', 'Reco objects / bin per event')
    return counts_by_tag


def detector_efficiencies(root, dataset, values, stats, stage):
    """Keep unmatched fiducial Gen photons in the denominator."""
    theta = np.rad2deg(np.arccos(np.clip(values['gamma']['cos_theta'], -1, 1)))
    regions = em_regions(theta)
    region_masks = [np.ones(len(theta), dtype=bool), *regions.values()]
    matched = stats['reco_count'] > 0
    tags = stats['detector_matched']
    tagged_matches = [tags.any(axis=1), *tags.T]
    ticklabels = ['All ranges\n' + r'$(0–180^\circ)$']
    for detector in DETECTORS:
        ranges = [rf'${low}–{high}^\circ$' for low, high in EM_REGIONS[detector]]
        ticklabels.append(detector + '\n' + '\n'.join(ranges))
    x = np.arange(len(region_masks))
    raw_counts = {}
    output = plot_path(root, '02_gamma', stage, dataset, values['source'],
                       values['collection'], 'matching_result')
    for range_name, _, _ in ENERGY_RANGES:
        raw_counts[range_name] = {}
        for tagged in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            kind = 'efficiency_vs_detector_flag' if tagged else 'efficiency_vs_detector_regional'
            raw_counts[range_name][kind] = {}
            for index, (population, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
                selected = values['energy_selections'][range_name][population]
                denominators = np.array([(selected & region).sum() for region in region_masks])
                numerators = np.array([(selected & region & (tag if tagged else matched)).sum()
                                       for region, tag in zip(region_masks, tagged_matches, strict=True)])
                shown = denominators > 0
                ratio = numerators[shown] / denominators[shown]
                lower, upper = clopper_pearson(numerators[shown], denominators[shown])
                ax.errorbar(x[shown], ratio,
                            yerr=np.maximum([ratio - lower, upper - ratio], 0),
                            fmt='o', color=f'C{index}', label=label,
                            markersize=9, alpha=0.7, capsize=3)
                raw_counts[range_name][kind][population] = dict(numerator=numerators, denominator=denominators)
            ax.set_xticks(x, ticklabels, fontsize=19)
            ax.set_xlabel('Region\n' + r'$(\theta)$', fontsize=22)
            ax.set_xlim(-0.5, len(x) - 0.5)
            ax.set_ylim(0, 1.35)
            ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
            ax.set_yticks(np.linspace(0, 1, 6))
            note = 'Gen fiducial'
            if ENERGY_LABELS[range_name]:
                note += '\n' + ENERGY_LABELS[range_name]
            ax.text(0.02, 0.89, note, transform=ax.transAxes, va='top', fontsize=20)
            if tagged and values['collection'] == 'Photon_PhotonConv':
                ax.text(0.02, 0.70 if ENERGY_LABELS[range_name] else 0.82, 'Conversions: no detector tag',
                        transform=ax.transAxes, va='top', fontsize=18)
            raw = values['source'] == 'nanoaod_raw_sdst'
            ylabel = ('Tagged coverage' if tagged else 'Angular EM coverage') if raw else (
                'Tagged efficiency' if tagged else 'Efficiency')
            filename = 'eff_det_flag' if tagged else 'eff_det_regional'
            suffix = '' if range_name == 'all' else f'_{range_name}'
            finish(fig, ax, output / f'{filename}{suffix}.png', dataset,
                   'Region\n' + r'$(\theta)$', ylabel)
    return raw_counts
