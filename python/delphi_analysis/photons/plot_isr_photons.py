"""Inclusive photon reconstruction and ISR energy study, without gen kinematic cuts."""
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from scipy.stats import beta

from delphi_analysis.plot_utils import add_delphi_label, sample_styles
from .data import MAX_ANGLE, SAMPLES, read_sample
from .plot_2d import POPULATIONS, energy_bins, plot_maps

COS_BINS = np.linspace(-1, 1, 41)
PHI_BINS = np.linspace(0, 2 * np.pi, 41)
COVERAGE = 0.6826894921370859
FIGURE_SIZE = (12, 8)


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


def finish(ax, title, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=title, loc='upper right', fontsize=12, title_fontsize=12,
              framealpha=0.9, edgecolor='none')
    ax.grid(alpha=0.2)
    add_delphi_label(ax)


def plot_population(output, kind, name, values_by_sample, styles):
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
        bins, xlabel = np.arange(-0.5, maximum + 1.5), 'Photons per event'
    peaks = []
    undefined = False
    for sample in SAMPLES:
        population = values_by_sample[sample][kind]
        values = population[name]
        if name == 'phi':
            undefined |= bool(np.any(np.isnan(values)))
            values = np.where(np.isnan(values), 0, values)
        n_events = len(population['multiplicity'])
        counts = count_curve(ax, values, bins, n_events, styles[sample][0], sample)
        peaks.append(counts.max() / n_events)
    if name == 'energy':
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(top=max(peaks) * 5)
    elif name == 'cos_theta' or name == 'phi':
        ax.set_yscale('log')
        ax.set_ylim(top=max(peaks) * 5)
    else:
        ax.set_ylim(0, max(peaks) * 1.5)
    ylabel = 'Fraction of events / bin' if name == 'multiplicity' else 'Photons / event / bin'
    finish(ax, POPULATIONS[kind], xlabel, ylabel)
    if undefined:
        ax.text(0.02, 0.05, 'Undefined φ displayed at 0', transform=ax.transAxes, fontsize=11)
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


def plot_stack(output, total_kind, name, values_by_sample):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), layout='constrained')
    axes[1, 2].set_visible(False)
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
        bins, xlabel = np.arange(-0.5, maximum + 1.5), 'Photons per event'
    for panel, sample in enumerate(SAMPLES):
        ax = axes.flat[panel]
        n_events = len(values_by_sample[sample][total_kind]['multiplicity'])
        bottom = np.zeros(len(bins) - 1)
        for kind, color, hatch in COMPONENTS[total_kind]:
            values = values_by_sample[sample][kind][name]
            if name == 'phi':
                values = np.where(np.isnan(values), 0, values)
            label = POPULATIONS[kind]
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
                        n_events, 'black', 'Total')
        else:
            total_values = values_by_sample[sample][total_kind][name]
            if name == 'phi':
                total_values = np.where(np.isnan(total_values), 0, total_values)
            assert np.allclose(bottom * n_events, np.histogram(total_values, bins=bins)[0])
        if name == 'energy':
            ax.set_xscale('log')
        if name in ('energy', 'cos_theta', 'phi'):
            ax.set_yscale('log')
            ax.set_ylim(top=max(0.1, bottom.max() * 5))
        ax.set_title(sample, fontsize=14)
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel('Fraction of events / bin' if name == 'multiplicity' else 'Photons / event / bin', fontsize=12)
        ax.tick_params(labelsize=10)
        ax.legend(loc='upper center' if name == 'cos_theta' else 'upper right',
                  fontsize=10, framealpha=0.9, edgecolor='none')
        ax.grid(alpha=0.2)
    title = 'DELPHI Simulation · ISR components' if total_kind == 'isr' else 'DELPHI Simulation · stable gen photons'
    if name == 'phi':
        title += ' · undefined φ displayed at 0'
    if name == 'multiplicity':
        title += ' · component distributions overlaid'
    fig.suptitle(title, fontsize=17)
    fig.savefig(output / f'{total_kind}_{name}.png', dpi=150)
    plt.close(fig)


def plot_efficiency(output, name, gen_kind, matched_kind, values_by_sample, styles, e_bins):
    bins = e_bins if name == 'energy' else np.linspace(-1, 1, 21)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for sample in SAMPLES:
        values = values_by_sample[sample]
        denominator = np.histogram(values[gen_kind][name], bins=bins)[0]
        numerator = np.histogram(values[matched_kind][name], bins=bins)[0]
        assert denominator.sum() == len(values[gen_kind][name])
        assert numerator.sum() == len(values[matched_kind][name])
        assert np.all(numerator <= denominator)
        valid = denominator > 0
        ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=valid)
        lower, upper = clopper_pearson(numerator, denominator)
        color = styles[sample][0]
        centers = np.sqrt(bins[:-1] * bins[1:]) if name == 'energy' else (bins[:-1] + bins[1:]) / 2
        ax.stairs(ratio, bins, color=color, label=sample)
        ax.errorbar(centers[valid], ratio[valid],
                    yerr=[ratio[valid] - lower[valid], upper[valid] - ratio[valid]],
                    fmt='.', color=color, capsize=2, linewidth=0.8)
    if name == 'energy':
        ax.set_xscale('log')
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.05)
    xlabel = r'Gen $E_\gamma$ [GeV]' if name == 'energy' else r'Gen $\cos\theta_\gamma$'
    title = 'All photons' if gen_kind == 'stable_gen' else 'All ISR photons'
    if 'truth' in matched_kind:
        title += ' · truth-linked'
    finish(ax, title, xlabel, 'Matching efficiency')
    fig.tight_layout()
    prefix = 'photon' if gen_kind == 'stable_gen' else 'isr'
    fig.savefig(output / f'{prefix}_efficiency_vs_{name}.png', dpi=150)
    plt.close(fig)


def plot_energy(output, values_by_sample, styles, reco_key):
    match_label = 'Truth-linked reco' if reco_key == 'truth_reco_energy' else 'Angular-associated reco'
    e_max = max(v['events']['isr_energy'].max() for v in values_by_sample.values())
    bins = np.r_[0, np.geomspace(1e-8, max(50, e_max), 61)]
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for sample in SAMPLES:
        total = values_by_sample[sample]['events']['isr_energy']
        count_curve(ax, total, bins, len(total), styles[sample][0], sample)
    ax.set_xscale('symlog', linthresh=1e-8)
    ax.set_yscale('log')
    ax.set_ylim(top=3)
    finish(ax, 'All ISR photons', 'Total ISR energy per event [GeV]', 'Fraction of events / bin')
    fig.tight_layout()
    fig.savefig(output / 'isr_total_energy_per_event.png', dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for sample in SAMPLES:
        event = values_by_sample[sample]['events']
        gen, reco = event['isr_energy'], event[reco_key]
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
        color = styles[sample][0]
        ax.stairs(ratio, bins, color=color, label=sample)
        ax.errorbar(centers[valid], ratio[valid], yerr=errors[valid], color=color, fmt='.', capsize=2)
    ax.set_xscale('symlog', linthresh=1e-8)
    ax.set_yscale('symlog', linthresh=0.1)
    finish(ax, match_label, 'Total ISR energy per event [GeV]', r'$\sum E_{\mathrm{matched\ reco}} / \sum E_{\mathrm{ISR}}$')
    fig.tight_layout()
    fig.savefig(output / 'isr_energy_recovery_vs_total_gen_energy.png', dpi=150)
    plt.close(fig)

    for name in ('radiated', 'recovered', 'difference'):
        fractions = {}
        for sample in SAMPLES:
            event = values_by_sample[sample]['events']
            energies = {'radiated': event['isr_energy'], 'recovered': event[reco_key],
                        'difference': event['isr_energy'] - event[reco_key]}
            fractions[sample] = energies[name] / event['cm_energy']
        limits = np.concatenate(list(fractions.values()))
        fraction_bins = np.linspace(min(0, limits.min()), np.nextafter(limits.max(), np.inf), 61)
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for sample, values in fractions.items():
            count_curve(ax, values, fraction_bins, len(values), styles[sample][0], sample)
        ax.set_yscale('log')
        ax.set_ylim(top=3)
        finish(ax, f'{match_label}: {name}', f'{name.capitalize()} energy / √s', 'Fraction of events / bin')
        fig.tight_layout()
        fig.savefig(output / f'isr_{name}_fraction_of_sqrts.png', dpi=150)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for index, sample in enumerate(SAMPLES):
        event = values_by_sample[sample]['events']
        gen, reco, cm = event['isr_energy'], event[reco_key], event['cm_energy']
        for offset, energy, hatch, label in ((-0.24, gen, '', 'All ISR'), (0, reco, '//', 'Matched reco'),
                                             (0.24, gen - reco, 'xx', 'ISR − matched reco')):
            fraction = 100 * energy / cm
            ax.bar(index + offset, fraction.mean(), 0.24, yerr=fraction.std(ddof=1) / np.sqrt(len(fraction)),
                   color=styles[sample][0], edgecolor='black', hatch=hatch, capsize=3,
                   label=label if index == 0 else None)
    ax.set_xticks(range(len(SAMPLES)), SAMPLES)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.5)
    finish(ax, match_label, '', 'Mean energy / √s per event [%]')
    fig.tight_layout()
    fig.savefig(output / 'isr_mean_energy_fractions_of_sqrts.png', dpi=150)
    plt.close(fig)


def plot_matching(output, values_by_sample, styles):
    for name in ('nearest_angle', 'energy_response'):
        all_values = np.concatenate([v['matching'][name] for v in values_by_sample.values()])
        positive = all_values[all_values > 0]
        edge = 10 ** np.floor(np.log10(positive.min()))
        bins = np.r_[0, np.geomspace(edge, np.nextafter(all_values.max(), np.inf), 61)]
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for sample in SAMPLES:
            values = values_by_sample[sample]['matching'][name]
            count_curve(ax, values, bins, len(values_by_sample[sample]['stable_gen']['multiplicity']), styles[sample][0], sample)
        ax.set_xscale('symlog', linthresh=edge)
        ax.set_yscale('log')
        reference = MAX_ANGLE if name == 'nearest_angle' else 1
        ax.axvline(reference, color='black', linestyle='--')
        xlabel = 'Nearest gen–reco opening angle [rad]' if name == 'nearest_angle' else r'$E_\gamma^{reco} / E_\gamma^{gen}$'
        finish(ax, 'Photon matching', xlabel, 'Photons / event / bin')
        fig.tight_layout()
        fig.savefig(output / f'photon_{name}.png', dpi=150)
        plt.close(fig)


def plot_isr_photons(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    plt.rcParams.update({'font.size': 17, 'axes.labelsize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
    study = output_root / 'isr_photons'
    directories = {name: study / folder for name, folder in (
        ('gen', '01_gen'), ('reco', '02_reco'), ('matching', '03_matching'),
        ('efficiency', '04_efficiency'),
    )}
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    values = {sample: read_sample(input_root, sample) for sample in SAMPLES}
    styles = sample_styles(SAMPLES)
    for kind, directory in (
        ('reco', 'reco'), ('matched_gen', 'matching'),
        ('matched_reco', 'matching'), ('isr_matched_gen', 'matching'), ('isr_matched_reco', 'matching'),
    ):
        for name in ('multiplicity', 'energy', 'cos_theta', 'phi'):
            plot_population(directories[directory], kind, name, values, styles)
    for name in ('multiplicity', 'energy', 'cos_theta', 'phi'):
        for kind in ('isr', 'stable_gen'):
            plot_stack(directories['gen'], kind, name, values)
    e_bins = energy_bins(values)
    for gen_kind, matched_kind in (('stable_gen', 'matched_gen'), ('isr', 'isr_matched_gen')):
        for name in ('energy', 'cos_theta'):
            plot_efficiency(directories['efficiency'], name, gen_kind, matched_kind, values, styles, e_bins)
    plot_energy(directories['efficiency'], values, styles, 'matched_reco_energy')
    truth_dir = directories['efficiency'] / 'truth_linked'
    truth_dir.mkdir(exist_ok=True)
    for gen_kind, matched_kind in (('stable_gen', 'truth_matched_gen'), ('isr', 'isr_truth_matched_gen')):
        for name in ('energy', 'cos_theta'):
            plot_efficiency(truth_dir, name, gen_kind, matched_kind, values, styles, e_bins)
    plot_energy(truth_dir, values, styles, 'truth_reco_energy')
    plot_matching(directories['matching'], values, styles)
    plot_maps(directories, values)
    print(f'plots: {study}')
