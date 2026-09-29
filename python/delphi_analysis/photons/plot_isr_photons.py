"""Inclusive photon reconstruction and ISR energy study, without gen kinematic cuts."""
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from scipy.stats import beta

from delphi_analysis.plot_utils import add_delphi_label
from .data import MAX_ANGLE, SAMPLES, read_sample
from .plot_2d import POPULATIONS, energy_bins, plot_maps

COS_BINS = np.linspace(-1, 1, 41)
PHI_BINS = np.linspace(0, 2 * np.pi, 41)
COVERAGE = 0.6826894921370859
FIGURE_SIZE = (11, 9)


def clopper_pearson(numerator, denominator):
    alpha = 1 - COVERAGE
    lower, upper = np.zeros_like(denominator, dtype=float), np.ones_like(denominator, dtype=float)
    nonzero, not_all = numerator > 0, numerator < denominator
    lower[nonzero] = beta.ppf(alpha / 2, numerator[nonzero], denominator[nonzero] - numerator[nonzero] + 1)
    upper[not_all] = beta.ppf(1 - alpha / 2, numerator[not_all] + 1, denominator[not_all] - numerator[not_all])
    return lower, upper


def count_curve(ax, values, bins, n_events, color, label):
    counts = np.histogram(values, bins=bins)[0]
    assert counts.sum() == len(values), "Histogram range must include every entry"
    ax.stairs(counts / n_events, bins, color=color, linewidth=1.8, label=label)
    centers = (bins[:-1] + bins[1:]) / 2
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / n_events,
                yerr=np.sqrt(counts[shown]) / n_events,
                fmt='.', color=color, markersize=3, capsize=1.5, linewidth=0.8)
    return counts


def finish(ax, sample, xlabel, ylabel, legend_loc='upper right'):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=sample, loc=legend_loc, framealpha=1, edgecolor='none')
    ax.grid(alpha=0.2)
    add_delphi_label(ax)


def plot_population(output, kind, name, values_by_sample, sample):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    if name == 'energy':
        energies = np.concatenate([v[kind][name] for v in values_by_sample.values()])
        bins = np.geomspace(10 ** np.floor(np.log10(energies.min())),
                            10 ** np.ceil(np.log10(energies.max())), 61)
        xlabel = r'$E_\gamma$ [GeV]'
    elif name == 'cos_theta':
        bins, xlabel = COS_BINS, r'$\cos\theta_\gamma$'
    elif name == 'phi':
        bins, xlabel = PHI_BINS, r'$\phi_\gamma$ [rad]'
    else:
        maximum = max(v[kind][name].max() for v in values_by_sample.values())
        bins, xlabel = np.arange(-0.5, maximum + 1.5), r'$N_\gamma$'
    population = values_by_sample[sample][kind]
    values = population[name]
    label = POPULATIONS[kind]
    if name == 'phi' and np.any(np.isnan(values)):
        label += '\n(φ undef. → 0)'
        values = np.where(np.isnan(values), 0, values)
    n_events = len(population['multiplicity'])
    counts = count_curve(ax, values, bins, n_events, 'C0', label)
    peak = counts.max() / n_events
    if name == 'energy':
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(top=peak * 5)
    elif name == 'cos_theta' or name == 'phi':
        ax.set_yscale('log')
        ax.set_ylim(top=peak * 5)
    else:
        ax.set_ylim(0, peak * 1.5)
    ylabel = r'$P(N_\gamma)$' if name == 'multiplicity' else 'Photons / event / bin'
    finish(ax, sample, xlabel, ylabel)
    fig.tight_layout()
    fig.savefig(output / f'{kind}_{name}.png', dpi=150)
    plt.close(fig)


COMPONENTS = {
    'isr': (('noncollinear_isr', 'C0', '..'), ('collinear_isr', 'C1', '///')),
    'stable_gen': (
        ('noncollinear_isr', 'C0', '..'), ('collinear_isr', 'C1', '///'),
        ('fsr', 'C2', '\\\\'), ('others', 'C3', 'xx'),
    ),
}


def plot_stack(output, total_kind, name, values_by_sample, sample):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    if name == 'energy':
        energies = np.concatenate([v[total_kind][name] for v in values_by_sample.values()])
        bins = np.geomspace(10 ** np.floor(np.log10(energies.min())),
                            10 ** np.ceil(np.log10(energies.max())), 61)
        xlabel = r'$E_\gamma$ [GeV]'
    elif name == 'cos_theta':
        bins, xlabel = COS_BINS, r'$\cos\theta_\gamma$'
    elif name == 'phi':
        bins, xlabel = PHI_BINS, r'$\phi_\gamma$ [rad]'
    else:
        maximum = max(v[total_kind][name].max() for v in values_by_sample.values())
        bins, xlabel = np.arange(-0.5, maximum + 1.5), r'$N_\gamma$'
    n_events = len(values_by_sample[sample][total_kind]['multiplicity'])
    bottom = np.zeros(len(bins) - 1)
    for kind, color, hatch in COMPONENTS[total_kind]:
        values = values_by_sample[sample][kind][name]
        undefined = name == 'phi' and np.any(np.isnan(values))
        if name == 'phi':
            values = np.where(np.isnan(values), 0, values)
        label = POPULATIONS[kind]
        if undefined:
            label += '\n(φ undef. → 0)'
        counts = np.histogram(values, bins=bins)[0]
        assert counts.sum() == len(values)
        if name == 'multiplicity':
            count_curve(ax, values, bins, n_events, color, label)
        else:
            heights = counts / n_events
            ax.stairs(bottom + heights, bins, baseline=bottom, fill=True,
                      facecolor=color, edgecolor=color, hatch=hatch,
                      alpha=0.75, label=label)
            shown = counts > 0
            centers = ((bins[:-1] + bins[1:]) / 2)[shown]
            ax.errorbar(centers, (bottom + heights)[shown],
                        yerr=np.sqrt((bottom * n_events + counts)[shown]) / n_events,
                        fmt='none', color=color, linewidth=0.8, capsize=1)
            bottom += heights
    if name == 'multiplicity':
        count_curve(ax, values_by_sample[sample][total_kind][name], bins,
                    n_events, 'black', f'Total {POPULATIONS[total_kind].lower()}')
    else:
        total_values = values_by_sample[sample][total_kind][name]
        if name == 'phi':
            total_values = np.where(np.isnan(total_values), 0, total_values)
        assert np.allclose(bottom * n_events, np.histogram(total_values, bins=bins)[0])
    if name == 'energy':
        ax.set_xscale('log')
    if name in ('energy', 'cos_theta', 'phi'):
        ax.set_yscale('log')
        ax.set_ylim(top=max(0.1, bottom.max() * (100 if name == 'energy' else 5)))
    ylabel = r'$P(N_\gamma)$' if name == 'multiplicity' else 'Photons / event / bin'
    finish(ax, sample, xlabel, ylabel, 'upper center' if name == 'cos_theta' else 'upper right')
    fig.tight_layout()
    fig.savefig(output / f'{total_kind}_{name}.png', dpi=150)
    plt.close(fig)


def plot_efficiency(output, name, gen_kind, matched_kind, values, sample, e_bins):
    bins = e_bins if name == 'energy' else np.linspace(-1, 1, 21)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    denominator = np.histogram(values[gen_kind][name], bins=bins)[0]
    numerator = np.histogram(values[matched_kind][name], bins=bins)[0]
    assert denominator.sum() == len(values[gen_kind][name])
    assert numerator.sum() == len(values[matched_kind][name])
    assert np.all(numerator <= denominator)
    valid = denominator > 0
    ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=valid)
    lower, upper = clopper_pearson(numerator, denominator)
    centers = np.sqrt(bins[:-1] * bins[1:]) if name == 'energy' else (bins[:-1] + bins[1:]) / 2
    label = 'Truth-linked match' if 'truth' in matched_kind else 'Angular match'
    ax.stairs(ratio, bins, color='C0', label=label)
    ax.errorbar(centers[valid], ratio[valid],
                yerr=[ratio[valid] - lower[valid], upper[valid] - ratio[valid]],
                fmt='.', color='C0', capsize=2, linewidth=0.8)
    if name == 'energy':
        ax.set_xscale('log')
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.05)
    xlabel = r'$E_\gamma^{\mathrm{gen}}$ [GeV]' if name == 'energy' else r'$\cos\theta_\gamma^{\mathrm{gen}}$'
    finish(ax, sample, xlabel, 'Efficiency')
    fig.tight_layout()
    prefix = 'photon' if gen_kind == 'stable_gen' else 'isr'
    fig.savefig(output / f'{prefix}_efficiency_vs_{name}.png', dpi=150)
    plt.close(fig)


def plot_energy(output, values_by_sample, sample, reco_key):
    match_label = 'Truth-linked reco' if reco_key == 'truth_reco_energy' else 'Angular-associated reco'
    e_max = max(v['events']['isr_energy'].max() for v in values_by_sample.values())
    bins = np.r_[0, np.geomspace(1e-8, max(50, e_max), 61)]
    event = values_by_sample[sample]['events']
    gen, reco, cm = event['isr_energy'], event[reco_key], event['cm_energy']
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    count_curve(ax, gen, bins, len(gen), 'C0', 'Generated ISR energy')
    ax.set_xscale('symlog', linthresh=1e-8)
    ax.set_yscale('log')
    ax.set_ylim(top=3)
    finish(ax, sample, r'$\sum E_\mathrm{ISR}^{\mathrm{gen}}$ [GeV]', 'Event fraction')
    fig.tight_layout()
    fig.savefig(output / 'isr_total_energy_per_event.png', dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    denominator = np.histogram(gen, bins=bins, weights=gen)[0]
    numerator = np.histogram(gen, bins=bins, weights=reco)[0]
    ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=denominator > 0)
    errors = np.full(len(ratio), np.nan)
    for index in range(len(ratio)):
        selected = (gen >= bins[index]) & (gen < bins[index + 1])
        if index == len(ratio) - 1:
            selected = (gen >= bins[index]) & (gen <= bins[index + 1])
        n = selected.sum()
        if n > 1 and denominator[index] > 0:
            errors[index] = np.std(reco[selected] - ratio[index] * gen[selected], ddof=1) * np.sqrt(n) / denominator[index]
    centers = (bins[:-1] + bins[1:]) / 2
    valid = np.isfinite(errors)
    ax.stairs(ratio, bins, color='C0', label=match_label)
    ax.errorbar(centers[valid], ratio[valid], yerr=errors[valid], color='C0', fmt='.', capsize=2)
    ax.set_xscale('symlog', linthresh=1e-8)
    ax.set_yscale('symlog', linthresh=0.1)
    finish(ax, sample, r'$\sum E_\mathrm{ISR}^{\mathrm{gen}}$ [GeV]',
           r'$\sum E_\mathrm{reco} / \sum E_\mathrm{ISR}^\mathrm{gen}$')
    fig.tight_layout()
    fig.savefig(output / 'isr_energy_recovery_vs_total_gen_energy.png', dpi=150)
    plt.close(fig)

    for name, energy, label in (
        ('radiated', gen, 'Generated ISR'),
        ('recovered', reco, match_label),
        ('difference', gen - reco, 'Generated ISR − matched reco'),
    ):
        fractions = {
            item: ({'radiated': values['events']['isr_energy'],
                    'recovered': values['events'][reco_key],
                    'difference': values['events']['isr_energy'] - values['events'][reco_key]}[name]
                   / values['events']['cm_energy'])
            for item, values in values_by_sample.items()
        }
        limits = np.concatenate(list(fractions.values()))
        fraction_bins = np.linspace(min(0, limits.min()), np.nextafter(limits.max(), np.inf), 61)
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        values = energy / cm
        count_curve(ax, values, fraction_bins, len(values), 'C0', label)
        ax.set_yscale('log')
        ax.set_ylim(top=3)
        finish(ax, sample, 'Energy / √s', 'Event fraction')
        fig.tight_layout()
        fig.savefig(output / f'isr_{name}_fraction_of_sqrts.png', dpi=150)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for index, (energy, hatch, label) in enumerate(((gen, '', 'Generated ISR'),
                                                    (reco, '//', match_label),
                                                    (gen - reco, 'xx', 'ISR − matched reco'))):
        fraction = 100 * energy / cm
        ax.bar(index, fraction.mean(), yerr=fraction.std(ddof=1) / np.sqrt(len(fraction)),
               color='C0', edgecolor='black', hatch=hatch, capsize=3, label=label)
    ax.set_xticks(range(3), ['Generated', 'Matched reco', 'Difference'])
    ax.set_ylim(top=ax.get_ylim()[1] * 1.5)
    finish(ax, sample, 'Energy component', 'Mean energy / √s [%]')
    fig.tight_layout()
    fig.savefig(output / 'isr_mean_energy_fractions_of_sqrts.png', dpi=150)
    plt.close(fig)


def plot_matching(output, values_by_sample, sample):
    for name in ('nearest_angle', 'energy_response'):
        all_values = np.concatenate([v['matching'][name] for v in values_by_sample.values()])
        positive = all_values[all_values > 0]
        edge = 10 ** np.floor(np.log10(positive.min()))
        bins = np.r_[0, np.geomspace(edge, np.nextafter(all_values.max(), np.inf), 61)]
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        values = values_by_sample[sample]['matching'][name]
        label = 'Gen photons' if name == 'nearest_angle' else 'Matched pairs'
        count_curve(ax, values, bins, len(values_by_sample[sample]['stable_gen']['multiplicity']), 'C0', label)
        ax.set_xscale('symlog', linthresh=edge)
        ax.set_yscale('log')
        reference = MAX_ANGLE if name == 'nearest_angle' else 1
        ax.axvline(reference, color='black', linestyle='--')
        xlabel = r'$\Delta\theta_{\min}$ [rad]' if name == 'nearest_angle' else r'$E_\gamma^{\mathrm{reco}} / E_\gamma^{\mathrm{gen}}$'
        finish(ax, sample, xlabel, 'Photons / event / bin' if name == 'nearest_angle' else 'Pairs / event / bin')
        fig.tight_layout()
        fig.savefig(output / f'photon_{name}.png', dpi=150)
        plt.close(fig)


def plot_isr_photons(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    study = output_root / 'isr_photons'
    directories = {name: study / folder for name, folder in (
        ('gen', '01_gen'), ('reco', '02_reco'), ('matching', '03_matching'),
        ('efficiency', '04_efficiency'),
    )}
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    values = {sample: read_sample(input_root, sample) for sample in SAMPLES}
    e_bins = energy_bins(values)
    for sample in SAMPLES:
        sample_dirs = {name: directory / sample for name, directory in directories.items()}
        truth_dir = directories['efficiency'] / 'truth_linked' / sample
        for directory in (*sample_dirs.values(), truth_dir):
            directory.mkdir(parents=True, exist_ok=True)
        for kind, directory in (
            ('reco', 'reco'), ('matched_gen', 'matching'),
            ('matched_reco', 'matching'), ('isr_matched_gen', 'matching'), ('isr_matched_reco', 'matching'),
        ):
            for name in ('multiplicity', 'energy', 'cos_theta', 'phi'):
                plot_population(sample_dirs[directory], kind, name, values, sample)
        for name in ('multiplicity', 'energy', 'cos_theta', 'phi'):
            for kind in ('isr', 'stable_gen'):
                plot_stack(sample_dirs['gen'], kind, name, values, sample)
        for gen_kind, matched_kind in (('stable_gen', 'matched_gen'), ('isr', 'isr_matched_gen')):
            for name in ('energy', 'cos_theta'):
                plot_efficiency(sample_dirs['efficiency'], name, gen_kind, matched_kind, values[sample], sample, e_bins)
        plot_energy(sample_dirs['efficiency'], values, sample, 'matched_reco_energy')
        for gen_kind, matched_kind in (('stable_gen', 'truth_matched_gen'), ('isr', 'isr_truth_matched_gen')):
            for name in ('energy', 'cos_theta'):
                plot_efficiency(truth_dir, name, gen_kind, matched_kind, values[sample], sample, e_bins)
        plot_energy(truth_dir, values, sample, 'truth_reco_energy')
        plot_matching(sample_dirs['matching'], values, sample)
    plot_maps(directories, values)
    print(f'plots: {study}')
