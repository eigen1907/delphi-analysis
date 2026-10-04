"""Gen-photon reconstruction, recovered energy, and reco-linked Sim diagnostics."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from matplotlib.ticker import MaxNLocator
from scipy.stats import beta

from .data import POPULATIONS, SAMPLES, read_sample

FIGURE_SIZE = (11, 9)
COVERAGE = 0.6826894921370859


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
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles, labels, title=sample, loc='upper right', framealpha=1, edgecolor='none')
    legend.get_title().set_fontweight('bold')
    legend.set_frame_on(True)
    legend.get_frame().set_facecolor('white')
    ax.grid(alpha=0.2)
    mh.label.exp_label(exp='DELPHI', llabel='Simulation', rlabel='LEP 1 (91.2 GeV)', loc=0, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def log_bins(entries, bins_per_decade):
    low = 10 ** np.floor(np.log10(entries.min()))
    high = np.nextafter(entries.max(), np.inf)
    return np.geomspace(low, high, int(np.ceil(bins_per_decade * np.log10(high / low))) + 1)


def count_plot(path, sample, entries, bins, normalization, xlabel, ylabel, label, linear, xscale='linear'):
    counts = np.histogram(entries, bins=bins)[0]
    assert counts.sum() == len(entries), 'Histogram must retain every selected entry'
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.stairs(counts / normalization, bins, color='C0', label=label)
    centers = np.sqrt(bins[:-1] * bins[1:]) if xscale == 'log' else (bins[:-1] + bins[1:]) / 2
    if xscale == 'symlog' and not linear:
        centers[0] = 0  # The first bin contains only unmatched photons with ratio zero.
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / normalization,
                yerr=np.sqrt(counts[shown]) / normalization,
                fmt='o', color='C0', markersize=4, capsize=2, linewidth=0.8)
    if not linear:
        ax.set_yscale('log')
        ax.set_ylim(top=counts.max() / normalization * 10)
        ax.set_xscale(xscale)
        if xscale == 'symlog':
            ax.set_xscale('symlog', linthresh=bins[1])  # Keep zero visible; display scale only.
    else:
        ax.set_ylim(0, 1.4 * (counts + np.sqrt(counts)).max() / normalization)
    ax.set_xlim(bins[0], bins[-1])
    if np.issubdtype(entries.dtype, np.integer):
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    finish(fig, ax, path, sample, xlabel, ylabel)


def efficiency_plot(path, sample, coordinates, bins, numerator, denominator, linear):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    centers = np.sqrt(bins[:-1] * bins[1:]) if coordinates == 'energy' and not linear else (bins[:-1] + bins[1:]) / 2
    assert np.all(numerator <= denominator)
    shown = denominator > 0
    ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=shown)
    lower, upper = clopper_pearson(numerator, denominator)
    ax.stairs(ratio, bins, baseline=None, color='C0', label='Gen → Reco')
    ax.errorbar(centers[shown], ratio[shown],
                yerr=[ratio[shown] - lower[shown], upper[shown] - ratio[shown]],
                fmt='o', color='C0', markersize=4, capsize=2, linewidth=0.8)
    if coordinates == 'energy' and not linear:
        ax.set_xscale('log')
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.25)
    ax.set_yticks(np.linspace(0, 1, 6))
    xlabel = r'$E_\gamma^{\rm gen}$ [GeV]' if coordinates == 'energy' else r'$\cos\theta_\gamma^{\rm gen}$'
    finish(fig, ax, path, sample, xlabel, 'Efficiency')


def recovery_profile(path, sample, coordinates, ratios, bins, xlabel, linear):
    """Unweighted mean ratio per Gen photon, including unmatched zeros, with SEM."""
    counts = np.histogram(coordinates, bins=bins)[0]
    sums = np.histogram(coordinates, bins=bins, weights=ratios)[0]
    squares = np.histogram(coordinates, bins=bins, weights=ratios ** 2)[0]
    assert counts.sum() == len(ratios)
    mean = np.divide(sums, counts, out=np.full(len(counts), np.nan), where=counts > 0)
    measured = counts > 1
    sem = np.sqrt(np.maximum(squares[measured] - sums[measured] ** 2 / counts[measured], 0)
                  / (counts[measured] * (counts[measured] - 1)))
    centers = (bins[:-1] + bins[1:]) / 2 if linear else np.sqrt(bins[:-1] * bins[1:])
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.stairs(mean, bins, baseline=None, color='C0', label='All gen (unmatched = 0)')
    ax.plot(centers[counts > 0], mean[counts > 0], 'o', color='C0', markersize=4)
    ax.errorbar(centers[measured], mean[measured], yerr=sem, fmt='none', color='C0', capsize=2)
    if not linear:
        ax.set_xscale('log')
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(bottom=0)
    finish(fig, ax, path, sample, xlabel, r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$')


def species_plot(path, sample, links, linear):
    codes = links['mass_code']
    groups = (codes == 21, np.abs(codes) == 2, (codes != 21) & (np.abs(codes) != 2))
    bottom = np.zeros(3)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for index, (selected, label) in enumerate(((links['depth'] == 0, 'Direct'), (links['depth'] > 0, 'Descendant'))):
        counts = np.array([np.count_nonzero(group & selected) for group in groups])
        ax.bar(np.arange(3), counts / len(codes), bottom=bottom, label=label,
               color=f'C{index}', hatch='' if index == 0 else '//')
        bottom += counts / len(codes)
    counts = np.array([np.count_nonzero(group) for group in groups])
    ax.errorbar(np.arange(3), bottom, yerr=np.sqrt(counts) / len(codes), fmt='o', color='black', capsize=2)
    ax.set_xticks(np.arange(3), (r'$\gamma$', r'$e^\pm$', 'Other'))
    if not linear:
        ax.set_yscale('log')
        ax.set_ylim(top=10)
    else:
        ax.set_ylim(0, 1.4)
    finish(fig, ax, path, sample, 'Linked Sim species', 'Reco fraction')


def plot_population(output, sample, pop, links, all_populations, all_links):
    """Nine observables; logarithmic figures also have fully linear siblings."""
    energies = np.concatenate([other['energy'] for other in all_populations])
    efficiency_energy = log_bins(energies, 3)
    linear_energy = np.arange(0, 10 * np.ceil(energies.max() / 10) + 1)
    with_reco = pop['reco_count'] > 0
    n_gen = len(pop['energy'])
    for linear in (False, True):
        suffix = '_linear' if linear else ''
        for coordinates, bins in (('energy', linear_energy if linear else efficiency_energy),
                                  ('cos_theta', np.linspace(-1, 1, 21))):
            if coordinates == 'cos_theta' and linear:
                continue  # Both cos(theta) figures already have only linear axes.
            gen_counts = np.histogram(pop[coordinates], bins=bins)[0]
            reco_counts = np.histogram(pop[coordinates][with_reco], bins=bins)[0]
            assert gen_counts.sum() == n_gen and reco_counts.sum() == with_reco.sum()
            efficiency_plot(output / f'gen_to_reco_efficiency_vs_{coordinates}{suffix}.png', sample,
                            coordinates, bins, reco_counts, gen_counts, linear)
            xlabel = r'$E_\gamma^{\rm gen}$ [GeV]' if coordinates == 'energy' else r'$\cos\theta_\gamma^{\rm gen}$'
            recovery_profile(output / f'energy_recovery_vs_{coordinates}{suffix}.png', sample,
                             pop[coordinates], pop['summed_ratio'], bins, xlabel, linear or coordinates == 'cos_theta')
        bins = np.arange(-0.5, max(other['reco_count'].max() for other in all_populations) + 1.5)
        count_plot(output / f'reco_gamma_multiplicity{suffix}.png', sample, pop['reco_count'], bins, n_gen,
                   r'$N_\gamma^{\rm reco}$ / gen', 'Gen fraction', 'All gen', linear)
        for field, filename, xlabel in (
            ('leading_ratio', 'leading_energy_response', r'$E_{\rm leading}^{\rm reco}/E_\gamma^{\rm gen}$'),
            ('summed_ratio', 'summed_energy_response', r'$\sum E_\gamma^{\rm reco}/E_\gamma^{\rm gen}$'),
        ):
            entries = pop[field] if field == 'summed_ratio' else pop[field][with_reco]
            all_entries = np.concatenate([other[field][other['reco_count'] > 0] for other in all_populations])
            assert np.all(np.isfinite(entries))
            bins = np.linspace(0, np.nextafter(all_entries.max(), np.inf), 81) if linear else log_bins(all_entries, 6)
            summed = field == 'summed_ratio'
            if summed and not linear:
                bins = np.r_[0, bins]
            count_plot(output / f'{filename}{suffix}.png', sample, entries, bins,
                       n_gen if summed else with_reco.sum(), xlabel, 'Gen fraction / bin',
                       'All gen (unmatched = 0)' if summed else 'Matched gen', linear, 'symlog' if summed else 'log')
        depth_bins = np.arange(-0.5, max(other['depth'].max() for other in all_links) + 1.5)
        count_plot(output / f'linked_sim_depth{suffix}.png', sample, links['depth'], depth_bins, len(links['depth']),
                   'Sim ancestry depth', 'Reco fraction', 'Truth-associated reco', linear)
        species_plot(output / f'linked_sim_species{suffix}.png', sample, links, linear)


def write_summary(output, values):
    summary = {}
    for sample, values_sample in values.items():
        populations = {}
        for name, pop in values_sample['populations'].items():
            n_gen = len(pop['energy'])
            matched = pop['reco_count'] > 0
            n_reco = int(matched.sum())
            lower, upper = clopper_pearson(np.array([n_reco]), np.array([n_gen]))
            links = values_sample['linked_sim'][name]
            codes, counts = np.unique(links['mass_code'], return_counts=True)
            ratio = pop['summed_ratio']
            populations[name] = dict(
                gen_photons=n_gen,
                gen_to_reco_efficiency_68_percent_cp=dict(numerator=n_reco, denominator=n_gen,
                    efficiency=n_reco / n_gen, lower=float(lower[0]), upper=float(upper[0])),
                associated_reco_candidates=int(pop['reco_count'].sum()),
                gen_with_multiple_reco=int(np.count_nonzero(pop['reco_count'] > 1)),
                energy_recovery=dict(mean_all_gen=float(ratio.mean()), mean_matched_gen=float(ratio[matched].mean()),
                    matched_quantiles_16_50_84=np.quantile(ratio[matched], [0.16, 0.5, 0.84]).tolist(),
                    gen_with_ratio_above_one=int(np.count_nonzero(ratio > 1)), max_ratio=float(ratio.max()),
                    total_gen_energy_GeV=float(pop['energy'].sum()),
                    total_associated_reco_energy_GeV=float(np.sum(pop['energy'] * ratio)),
                    ratio_of_total_energies=float(np.sum(pop['energy'] * ratio) / pop['energy'].sum())),
                linked_sim=dict(direct=int(np.count_nonzero(links['depth'] == 0)),
                    descendant=int(np.count_nonzero(links['depth'] > 0)),
                    nonterminal=int(np.count_nonzero(links['has_children'])),
                    mass_codes=dict(zip(map(str, codes), map(int, counts), strict=True)),
                    depth_counts=np.bincount(links['depth']).tolist()),
                saved_lineage_diagnostic=dict(gen_with_lineage=int(np.count_nonzero(pop['sim_count'])),
                    gen_with_secondary_sim=int(np.count_nonzero(pop['descendant_count']))),
            )
        summary[sample] = dict(events=values_sample['events'], association=values_sample['association'], populations=populations)
    (output / 'study_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')


def plot_isr_photons(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    study = output_root / 'photon_study'
    values = {sample: read_sample(input_root, sample) for sample in SAMPLES}
    for population in POPULATIONS:
        all_populations = [values[sample]['populations'][population] for sample in SAMPLES]
        all_links = [values[sample]['linked_sim'][population] for sample in SAMPLES]
        for sample in SAMPLES:
            output = study / population / sample
            output.mkdir(parents=True, exist_ok=True)
            plot_population(output, sample, values[sample]['populations'][population],
                            values[sample]['linked_sim'][population], all_populations, all_links)
    write_summary(study, values)
    print(f'plots: {study}')
