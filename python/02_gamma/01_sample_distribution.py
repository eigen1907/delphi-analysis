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


def energy_bins(entries, logarithmic=False):
    maximum = np.nextafter(entries.max(), np.inf)
    if logarithmic:
        positive = entries[entries > 0]
        low = positive.min() / 2
        count = max(1, int(np.ceil(6 * np.log10(maximum / low))))
        return np.r_[0, np.geomspace(low, maximum, count + 1)]
    return np.r_[ENERGY_EDGES[ENERGY_EDGES < maximum], maximum]


def bin_centers(ax, bins, logarithmic=False):
    ax.set_xlim(bins[0], bins[-1])
    if logarithmic:
        # The first interval contains zero only, rather than dropping zero-energy candidates.
        ax.set_xscale('symlog', linthresh=bins[1])
        centers = np.r_[0, np.sqrt(bins[1:-1] * bins[2:])]
        # Space the labels around zero; the bins and all entries are unchanged.
        ax.set_xticks(np.r_[0, ax.get_xticks()[2:]])
    else:
        centers = (bins[:-1] + bins[1:]) / 2
    return centers


def gen_spectrum(ax, entries, category, components, bins, events, logarithmic):
    centers = bin_centers(ax, bins, logarithmic)
    bottom = np.zeros(len(bins) - 1)
    for key, label in components:
        index = GEN_COMPONENTS.index((key, label))
        selected = entries[(category == key) & np.isfinite(entries)]
        counts = np.histogram(selected, bins=bins)[0]
        assert counts.sum() == len(selected), 'Spectrum must retain all finite entries'
        ax.stairs(bottom + counts / events, bins, baseline=bottom, fill=True,
                  label=f'Gen {label}', color=f'C{index}', linewidth=0.9)
        bottom += counts / events
    total = np.histogram(entries[np.isfinite(entries)], bins=bins)[0]
    assert np.allclose(bottom * events, total)
    shown = total > 0
    ax.errorbar(centers[shown], bottom[shown], yerr=np.sqrt(total[shown]) / events,
                fmt='none', ecolor='grey', capsize=2, elinewidth=0.8)


def reco_spectrum(ax, entries, bins, events, logarithmic, label):
    """One inclusive Reco curve, with no truth-origin selection."""
    centers = bin_centers(ax, bins, logarithmic)
    entries = entries[np.isfinite(entries)]
    counts = np.histogram(entries, bins=bins)[0]
    assert counts.sum() == len(entries)
    ax.stairs(counts / events, bins, fill=True, facecolor='none', edgecolor='black',
              hatch='///', linewidth=2.1, label=f'Reco {label}')
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / events, yerr=np.sqrt(counts[shown]) / events,
                fmt='none', ecolor='black', capsize=2, elinewidth=0.8)
    return counts / events


def sample_distributions(study, sample, values, bins):
    """Change only the Gen selection; the inclusive Reco curve stays the same."""
    gamma = values['gamma']
    candidates = values['reco_distributions']
    for gen_name in GEN_SELECTIONS:
        output = plot_path(study, '02_gamma', '01_sample_distribution', sample, values['source'],
                           values['collection']) / gen_name
        output.mkdir(parents=True, exist_ok=True)
        gen_selected = values['selections'][gen_name]
        for coordinate in ('cos_theta', 'energy'):
            for log in (False, True):
                log_energy = log and coordinate == 'energy'
                edges = bins['energy_log'] if log_energy else bins[coordinate]
                gen_entries = gamma[coordinate]
                reco_entries = candidates[coordinate]
                if coordinate == 'energy' and not log:
                    # E >= 50 GeV goes into the labelled overflow bin, never out of the study.
                    last = np.nextafter(edges[-1], -np.inf)
                    gen_entries = np.minimum(gen_entries, last)
                    reco_entries = np.minimum(reco_entries, last)
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                gen_spectrum(ax, gen_entries[gen_selected], gamma['category'][gen_selected],
                             GEN_STACKS[gen_name], edges, values['events'], log_energy)
                reco_height = reco_spectrum(ax, reco_entries, edges, values['events'], log_energy, values['reco_label'])
                all_gen_entries = gen_entries[values['selections']['gen_gamma']]
                all_gen_counts = np.histogram(all_gen_entries[np.isfinite(all_gen_entries)], bins=edges)[0]
                # Use one y range for all four Gen selections.
                height = max(all_gen_counts.max() / values['events'], reco_height.max(), 1 / values['events'])
                if log:
                    ax.set_yscale('log')
                    ax.set_ylim(0.5 / values['events'], height * 1000)
                else:
                    ax.set_ylim(0, height * 1.8)
                if coordinate == 'energy' and not log:
                    ax.set_xticks([0, 10, 20, 30, 40, 55], ['0', '10', '20', '30', '40', r'$\geq 50$'])
                suffix = '_log' if log else ''
                xlabel = r'$E_\gamma$ [GeV]' if coordinate == 'energy' else r'$\cos\theta$'
                finish(fig, ax, output / f'{coordinate}_{sample}{suffix}.png', sample, xlabel,
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
        finish(fig, ax, output / f'multiplicity_{sample}{suffix}.png',
               sample, r'$N_\gamma$' if values['source'] == 'nanoaod' else r'$N$', 'Events')
