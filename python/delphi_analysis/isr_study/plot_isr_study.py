"""Read the five samples, then draw distributions, link checks, and efficiencies."""
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from matplotlib.ticker import MaxNLocator
from scipy.stats import beta

from .angle_matching import accumulate, angular_stats, choose_cut, new_accumulator, plot_angle_matching
from .data import GEN_SELECTIONS, LINK_STATUSES, RECO_DEFINITIONS, SAMPLES, read_sample

FIGURE_SIZE = (11, 9)
COVERAGE = 0.6826894921370859
ENERGY_EDGES = np.array([0, 0.1, 0.2, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50])
GEN_COMPONENTS = (('beam_isr', 'BeamISR'), ('nonbeam_isr', 'NonBeamISR'),
                  ('fsr', 'FSR'), ('decayed', 'Decayed'))
GEN_STACKS = dict(gen_gamma=GEN_COMPONENTS, gen_isr=GEN_COMPONENTS[:2],
                  gen_no_isr=GEN_COMPONENTS[2:], gen_isr_non_beam=GEN_COMPONENTS[1:2])
RECO_LABELS = ('Photon', 'Photon+conv')
XLABELS = dict(energy=r'$E_\gamma^{\rm gen}$ [GeV]', cos_theta=r'$\cos\theta^{\rm gen}$')


def clopper_pearson(numerator, denominator):
    alpha = 1 - COVERAGE
    lower, upper = np.zeros_like(denominator, dtype=float), np.ones_like(denominator, dtype=float)
    nonzero, not_all = numerator > 0, numerator < denominator
    lower[nonzero] = beta.ppf(alpha / 2, numerator[nonzero], denominator[nonzero] - numerator[nonzero] + 1)
    upper[not_all] = beta.ppf(1 - alpha / 2, numerator[not_all] + 1, denominator[not_all] - numerator[not_all])
    return lower, upper


def finish(fig, ax, path, sample, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    legend = ax.get_legend()
    if legend is None:
        legend = ax.legend(title=sample, loc='best', frameon=True, framealpha=1, edgecolor='none')
    legend.get_title().set_fontweight('bold')
    legend.get_frame().set_facecolor('white')
    ax.grid(alpha=0.2)
    mh.label.exp_label(exp='DELPHI', llabel='Simulation', rlabel='LEP 1 (91.2 GeV)', loc=0, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def energy_bins(entries, logarithmic=False):
    maximum = np.nextafter(entries.max(), np.inf)
    if logarithmic:
        positive = entries[entries > 0]
        low = positive.min() / 2
        count = max(1, int(np.ceil(6 * np.log10(maximum / low))))
        return np.r_[0, np.geomspace(low, maximum, count + 1)]
    return np.r_[ENERGY_EDGES[ENERGY_EDGES < maximum], maximum]


def bin_centers(ax, bins, logarithmic=False):
    if logarithmic:
        # The first interval contains zero only, rather than dropping zero-energy candidates.
        ax.set_xscale('symlog', linthresh=bins[1])
        centers = np.r_[0, np.sqrt(bins[1:-1] * bins[2:])]
    else:
        centers = (bins[:-1] + bins[1:]) / 2
    ax.set_xlim(bins[0], bins[-1])
    if logarithmic:
        # Space the labels around zero; the bins and all entries are unchanged.
        ax.set_xticks(np.r_[0, ax.get_xticks()[2:]])
    return centers


def stacked_spectrum(ax, entries, category, components, bins, events, logarithmic, reco):
    centers = bin_centers(ax, bins, logarithmic)
    bottom = np.zeros(len(bins) - 1)
    for key, label in components:
        index = {'beam_isr': 0, 'nonbeam_isr': 1, 'fsr': 2, 'decayed': 3,
                 'gamma': 0, 'gamma_conv': 1}[key]
        selected = entries[(category == key) & np.isfinite(entries)]
        counts = np.histogram(selected, bins=bins)[0]
        assert counts.sum() == len(selected), 'Spectrum must retain all finite entries'
        ax.stairs(bottom + counts / events, bins, baseline=bottom, fill=True,
                  label=f'Reco {label}' if reco else f'Gen {label}',
                  facecolor='none' if reco else f'C{index}',
                  edgecolor='black' if reco else f'C{index}',
                  hatch=('///', '\\\\\\')[index] if reco else None,
                  linewidth=1.8 if reco else 0.7)
        bottom += counts / events
    total = np.histogram(entries[np.isfinite(entries)], bins=bins)[0]
    assert np.allclose(bottom * events, total)
    shown = total > 0
    ax.errorbar(centers[shown], bottom[shown], yerr=np.sqrt(total[shown]) / events,
                fmt='none', ecolor='black' if reco else 'grey', capsize=2, elinewidth=0.8)
    return bottom


def sample_distributions(study, sample, values, bins):
    """Change only the Gen selection; the inclusive Reco stack stays the same."""
    gamma = values['gamma']
    candidates = values['reco_distributions']['reco_gamma_plus_conv']
    for gen_name in GEN_SELECTIONS:
        output = study / '01_sample_distribution' / gen_name
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
                stacked_spectrum(ax, gen_entries[gen_selected], gamma['category'][gen_selected],
                                 GEN_STACKS[gen_name], edges, values['events'], log_energy, False)
                reco_height = stacked_spectrum(ax, reco_entries, candidates['kind'],
                                               (('gamma', r'$\gamma$'), ('gamma_conv', r'$\gamma_{\rm conv}$')),
                                               edges, values['events'], log_energy, True)
                all_gen_counts = np.histogram(gen_entries[np.isfinite(gen_entries)], bins=edges)[0]
                # Use one y range for all four Gen selections.
                height = max(all_gen_counts.max() / values['events'], reco_height.max(), 1 / values['events'])
                if log:
                    ax.set_yscale('log')
                    ax.set_ylim(0.5 / values['events'], height * 1000)
                else:
                    ax.set_ylim(0, height * 1.8)
                ax.legend(title=sample, loc='upper center', ncol=2,
                          frameon=True, framealpha=1, edgecolor='none')
                if coordinate == 'energy' and not log:
                    ax.set_xticks([0, 10, 20, 30, 40, 55], ['0', '10', '20', '30', '40', r'$\geq 50$'])
                suffix = '_log' if log else ''
                xlabel = r'$E_\gamma$ [GeV]' if coordinate == 'energy' else r'$\cos\theta$'
                finish(fig, ax, output / f'{coordinate}_{sample}{suffix}.png', sample, xlabel,
                       r'$N_\gamma$ / bin per event')
    event_multiplicity(study, sample, values)


def event_multiplicity(study, sample, values):
    """Photon counts per event, including N_gamma = 0; overlapping curves are not stacked."""
    populations = (('all_gamma', r'Gen all $\gamma$'), ('all_isr', 'Gen ISR'),
                   ('beam_isr', 'Gen BeamISR'), ('nonbeam_isr', 'Gen NonBeamISR'),
                   ('fsr', 'Gen FSR'), ('decayed', 'Gen Decayed'))
    entries = [(values['gen_event_counts'][name], label, f'C{index}', '-')
               for index, (name, label) in enumerate(populations)]
    entries += [(values['reco_event_counts'][name], label, f'C{index + 6}', '--')
                for index, (name, label) in enumerate((('reco_gamma', r'Reco $\gamma$'),
                                                     ('reco_gamma_plus_conv', r'Reco $\gamma$+conv')))]
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
        ax.legend(title=sample, loc='upper center', ncol=2, frameon=True, framealpha=1, edgecolor='none')
        suffix = '_log' if log else ''
        finish(fig, ax, study / '01_sample_distribution' / f'multiplicity_{sample}{suffix}.png',
               sample, r'$N_\gamma$', 'Events')


def truth_link_validation(study, sample, values):
    labels = ('Agree', 'Forward\nonly', 'Reverse\nonly', 'Conflict', 'No origin')
    for gen_name in GEN_SELECTIONS:
        output = study / '02_truth_link_matching_validation' / gen_name
        output.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        x = np.arange(len(labels))
        maximum = 0
        for index, reco_name in enumerate(RECO_DEFINITIONS):
            counts = np.array([values['matching_validation'][gen_name][reco_name][key] for key in LINK_STATUSES])
            positions = x + (index - 0.5) * 0.3
            ax.bar(positions, counts, width=0.3, color=f'C{index}', label=RECO_LABELS[index])
            ax.errorbar(positions, counts, yerr=np.sqrt(counts), fmt='none', color=f'C{index}', capsize=2)
            maximum = max(maximum, counts.max())
        ax.set_xticks(x, labels)
        ax.set_yscale('log')
        ax.set_ylim(top=max(maximum * 10, 1))
        finish(fig, ax, output / f'truth_link_status_{sample}.png', sample, 'Link status', 'Reco objects')


def profile_points(coordinate, response, bins):
    defined = np.isfinite(response)
    count = np.histogram(coordinate[defined], bins=bins)[0]
    sums = np.histogram(coordinate[defined], bins=bins, weights=response[defined])[0]
    squares = np.histogram(coordinate[defined], bins=bins, weights=response[defined] ** 2)[0]
    assert count.sum() == defined.sum()
    mean = np.divide(sums, count, out=np.full(len(count), np.nan), where=count > 0)
    error = np.full(len(count), np.nan)
    measured = count > 1
    error[measured] = np.sqrt(np.maximum(squares[measured] - sums[measured] ** 2 / count[measured], 0)
                              / (count[measured] * (count[measured] - 1)))
    return mean, error


def matching_results(study, sample, values, bins, cut_degrees):
    gamma = values['gamma']
    series = []
    for source, method, line in (('reco', 'Truth', '-'), ('angular', f'Angle {cut_degrees:g}°', '--')):
        for index, (name, label) in enumerate(zip(RECO_DEFINITIONS, RECO_LABELS, strict=True)):
            series.append((values[source][name], f'{label} · {method}', f'C{index}', line))
    for gen_name in GEN_SELECTIONS:
        output = study / '04_matching_result' / gen_name
        output.mkdir(parents=True, exist_ok=True)
        selected = values['selections'][gen_name]
        for coordinate in ('cos_theta', 'energy'):
            edges = bins[coordinate]
            centers = (edges[:-1] + edges[1:]) / 2
            x = gamma[coordinate][selected]
            denominator = np.histogram(x, bins=edges)[0]
            assert denominator.sum() == selected.sum()
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for stats, label, color, line in series:
                matched = stats['reco_count'][selected] > 0
                numerator = np.histogram(x[matched], bins=edges)[0]
                assert numerator.sum() == matched.sum() and np.all(numerator <= denominator)
                shown = denominator > 0
                ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=shown)
                lower, upper = clopper_pearson(numerator, denominator)
                ax.stairs(ratio, edges, baseline=None, color=color, linestyle=line, label=label)
                ax.errorbar(centers[shown], ratio[shown], yerr=np.maximum([ratio[shown] - lower[shown], upper[shown] - ratio[shown]], 0),
                            fmt='none', color=color, capsize=2)
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0, 1.5)
            ax.set_yticks(np.linspace(0, 1, 6))
            ax.legend(title=sample, loc='upper center', ncol=2, frameon=True, framealpha=1, edgecolor='none')
            finish(fig, ax, output / f'efficiency_vs_{coordinate}_{sample}.png', sample, XLABELS[coordinate], 'Efficiency')
            for log_y in (False, True):
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                means, upper_errors = [], []
                for stats, label, color, line in series:
                    mean, error = profile_points(x, stats['energy_ratio'][selected], edges)
                    means.append(mean)
                    upper_errors.append(mean + np.nan_to_num(error, nan=0))
                    ax.stairs(mean, edges, baseline=None, color=color, linestyle=line, label=label)
                    ax.errorbar(centers, mean, yerr=error, fmt='none', color=color, capsize=2)
                ax.set_xlim(edges[0], edges[-1])
                if log_y:
                    ax.set_yscale('log')
                    positive = np.concatenate(means)
                    positive = positive[np.isfinite(positive) & (positive > 0)]
                    ax.set_ylim(positive.min() / 3, np.nanmax(upper_errors) * 30)
                else:
                    ax.set_ylim(0, np.nanmax(upper_errors) * 1.65)
                ax.legend(title=sample, loc='upper center', ncol=2, frameon=True, framealpha=1, edgecolor='none')
                suffix = '_log' if log_y else ''
                finish(fig, ax, output / f'energy_response_vs_{coordinate}_{sample}{suffix}.png', sample,
                       XLABELS[coordinate], r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$')


def plot_isr_study(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    study = output_root / 'isr_study'
    study.mkdir(parents=True, exist_ok=True)
    # Read each event once: Gen selection, truth association, and angular scan.
    angles, values = {}, {}
    for sample in SAMPLES:
        angles[sample] = new_accumulator()
        values[sample] = read_sample(input_root, sample, lambda event: accumulate(angles[sample], event))
    all_energy = np.concatenate([np.r_[value['gamma']['energy'], value['reco_distributions']['reco_gamma_plus_conv']['energy']]
                                 for value in values.values()])
    bins = dict(energy=np.r_[ENERGY_EDGES, 60], energy_log=energy_bins(all_energy, True),
                cos_theta=np.linspace(-1, 1, 21))
    gen_bins = dict(bins, energy=energy_bins(np.concatenate([value['gamma']['energy'] for value in values.values()])))
    # Set the angular reference from pooled non-beam ISR truth retention.
    cut = choose_cut(angles.values())
    for sample in SAMPLES:
        values[sample]['angular'] = angular_stats(angles[sample], cut)
    # Draw the same four sections for every sample.
    for sample, value in values.items():
        sample_distributions(study, sample, value, bins)
        truth_link_validation(study, sample, value)
        plot_angle_matching(study / '03_angular_matching_validation', sample, angles[sample], cut)
        matching_results(study, sample, value, gen_bins, cut)
        print(f'plots: {sample}', flush=True)
    print(f'ISR study: {study}; opening-angle reference {cut:g} deg', flush=True)
