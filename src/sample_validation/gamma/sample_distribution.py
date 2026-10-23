"""Overlay stable Gen populations and the inclusive Reco collection."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from common.plot_utils import FIGURE_SIZE, finish, plot_path
from .data import GEN_SELECTIONS

ENERGY_EDGES = np.array([0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50])

COS_THETA_EDGES = np.linspace(-1, 1, 101)

GEN_COMPONENTS = (('beam_isr', r'$\gamma_{\mathrm{BeamISR}}$'),
                  ('nonbeam_isr', r'$\gamma_{\mathrm{NonBeamISR}}$'),
                  ('fsr', r'$\gamma_{\mathrm{FSR}}$'), ('decayed', r'$\gamma_{\mathrm{Decayed}}$'))

GEN_STACKS = dict(gen_gamma=GEN_COMPONENTS, gen_isr=GEN_COMPONENTS[:2],
                  gen_no_isr=GEN_COMPONENTS[2:], gen_isr_non_beam=GEN_COMPONENTS[1:2])


def log_energy_bins(entries):
    positive = entries[np.isfinite(entries) & (entries > 0)]
    low = positive.min() / 2 if len(positive) else 0.1
    high = np.nextafter(positive.max(), np.inf) if len(positive) else 60
    count = max(1, int(np.ceil(6 * np.log10(high / low))))
    return np.r_[0, np.geomspace(low, high, count + 1)]


def bin_centers(ax, bins, logarithmic=False):
    ax.set_xlim(bins[0], bins[-1])
    if logarithmic:
        # Preserve zero-energy candidates in a linear first bin, as in the old log-x plot.
        ax.set_xscale('symlog', linthresh=bins[1])
        ax.set_xticks(np.r_[0, ax.get_xticks()[2:]])
        return np.r_[0, np.sqrt(bins[1:-1] * bins[2:])]
    return (bins[:-1] + bins[1:]) / 2


def gen_spectrum(ax, entries, category, components, bins, events, logarithmic=False):
    centers = bin_centers(ax, bins, logarithmic)
    bottom = np.zeros(len(bins) - 1)
    for key, label in components:
        index = GEN_COMPONENTS.index((key, label))
        selected = entries[(category == key) & np.isfinite(entries)]
        counts = np.histogram(selected, bins=bins)[0]
        assert counts.sum() == len(selected), 'Spectrum must retain all finite entries'
        ax.stairs(bottom + counts / events, bins, baseline=bottom, fill=True,
                  label=f'Gen {label}', color=f'C{index}', linewidth=0.9)
        bottom = bottom + counts / events
    total = np.histogram(entries[np.isfinite(entries)], bins=bins)[0]
    assert np.allclose(bottom * events, total)
    shown = total > 0
    ax.errorbar(centers[shown], bottom[shown], yerr=np.sqrt(total[shown]) / events,
                fmt='none', ecolor='grey', capsize=2, elinewidth=0.8)


def reco_spectrum(ax, entries, bins, events, label, is_conversion=None, logarithmic=False):
    """Stack native photons and conversions without changing the inclusive total."""
    centers = bin_centers(ax, bins, logarithmic)
    finite = np.isfinite(entries)
    entries = entries[finite]
    counts = np.histogram(entries, bins=bins)[0]
    assert counts.sum() == len(entries)
    if is_conversion is None:
        components = [(counts, f'Reco {label}', '///')]
    else:
        is_conversion = np.asarray(is_conversion, dtype=bool)[finite]
        photons = np.histogram(entries[~is_conversion], bins=bins)[0]
        conversions = np.histogram(entries[is_conversion], bins=bins)[0]
        assert np.array_equal(photons + conversions, counts), 'Reco stack must preserve the inclusive spectrum'
        components = [(conversions, r'Reco $\gamma_{\mathrm{Conv}}$', '\\' * 3),
                      (photons, r'Reco $\gamma$', '///')]
    bottom = np.zeros(len(counts))
    for component, component_label, hatch in components:
        ax.stairs(bottom + component / events, bins, baseline=bottom, fill=True,
                  facecolor='none', edgecolor='black', hatch=hatch, linewidth=2.1,
                  label=component_label)
        bottom = bottom + component / events
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / events, yerr=np.sqrt(counts[shown]) / events,
                fmt='none', ecolor='black', capsize=2, elinewidth=0.8)
    return counts / events


def sample_distributions(study, sample, values, bins):
    """Change only the Gen selection; the inclusive Reco population stays the same."""
    gamma = values['gamma']
    candidates = values['reco_distributions']
    is_conversion = (candidates['is_conversion']
                     if values['source'] == 'nanoaod' and values['collection'] == 'Photon_PhotonConv'
                     else None)
    energy_log_bins = log_energy_bins(np.r_[gamma['energy'], candidates['energy']])
    for gen_name in GEN_SELECTIONS:
        output = plot_path(study, '02_gamma', '01_sample_distribution', sample, values['source'],
                           values['collection']) / gen_name
        output.mkdir(parents=True, exist_ok=True)
        gen_selected = values['selections'][gen_name]
        for coordinate in ('cos_theta', 'energy'):
            scales = ((False, False), (True, False), (True, True)) if coordinate == 'energy' else (
                (False, False), (True, False))
            for log, log_x in scales:
                edges = energy_log_bins if log_x else bins[coordinate]
                gen_entries = gamma[coordinate]
                reco_entries = candidates[coordinate]
                if coordinate == 'energy' and not log_x:
                    # E >= 50 GeV goes into the labelled overflow bin, never out of the study.
                    last = np.nextafter(edges[-1], -np.inf)
                    gen_entries = np.minimum(gen_entries, last)
                    reco_entries = np.minimum(reco_entries, last)
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                gen_spectrum(ax, gen_entries[gen_selected], gamma['category'][gen_selected],
                             GEN_STACKS[gen_name], edges, values['events'], log_x)
                reco_height = reco_spectrum(ax, reco_entries, edges, values['events'],
                                            values['reco_label'], is_conversion, log_x)
                all_gen_entries = gen_entries[values['selections']['gen_gamma']]
                all_gen_counts = np.histogram(all_gen_entries[np.isfinite(all_gen_entries)], bins=edges)[0]
                # Use one y range for all four Gen selections.
                height = max(all_gen_counts.max() / values['events'], reco_height.max(), 1 / values['events'])
                if log:
                    ax.set_yscale('log')
                    ax.set_ylim(0.5 / values['events'], height * 1000)
                else:
                    ax.set_ylim(0, height * 1.8)
                if coordinate == 'energy' and not log_x:
                    ax.set_xticks([0, 10, 20, 30, 40, 55], ['0', '10', '20', '30', '40', r'$\geq 50$'])
                suffix = '_loglog' if log_x else '_log' if log else ''
                xlabel = r'$E_\gamma$ [GeV]' if coordinate == 'energy' else r'$\cos\theta$'
                finish(fig, ax, output / f'{coordinate}{suffix}.png', sample, xlabel,
                       r'$N_\gamma$ / bin per event' if values['source'] == 'nanoaod' else 'Objects / bin per event')
        event_multiplicity(output, sample, values, GEN_STACKS[gen_name])


def event_multiplicity(output, sample, values, components):
    """Photon counts per event, including N_gamma = 0; overlapping curves are not stacked."""
    entries = [(values['gen_event_counts'][name], f'Gen {label}', f'C{GEN_COMPONENTS.index((name, label))}', '-')
               for name, label in components]
    entries.append((values['reco_event_counts'], f"Reco {values['reco_label']}", 'C4', '--'))
    bins = np.arange(-0.5, max(item[0].max() for item in entries) + 1.5)
    centers = (bins[:-1] + bins[1:]) / 2
    for log in (False, True):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        height = 0
        for numbers, label, color, line in entries:
            counts = np.histogram(numbers, bins=bins)[0]
            assert counts.sum() == len(numbers) == values['events']
            shown = counts > 0
            ax.stairs(counts, bins, label=label, color=color, linestyle=line)
            ax.errorbar(centers[shown], counts[shown], yerr=np.sqrt(counts[shown]),
                        fmt='none', color=color, capsize=2)
            height = max(height, (counts + np.sqrt(counts)).max())
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlim(bins[0], bins[-1])
        if log:
            ax.set_yscale('log')
            ax.set_ylim(0.5, height * 100)
        else:
            ax.set_ylim(0, height * 1.8)
        suffix = '_log' if log else ''
        finish(fig, ax, output / f'multiplicity{suffix}.png',
               sample, r'$N_\gamma$' if values['source'] == 'nanoaod' else r'$N$', 'Events')
