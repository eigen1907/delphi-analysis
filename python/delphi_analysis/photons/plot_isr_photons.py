"""The same Gen-unit truth study for photons, ISR, and ISR without beam photons."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import MaxNLocator
from scipy.stats import beta

from .data import POPULATIONS, SAMPLES, read_sample

FIGURE_SIZE = (11, 9)
COVERAGE = 0.6826894921370859
COUNT_LABEL = r'$N_\gamma$ / bin per event'
COS_BINS = np.linspace(-1, 1, 41)
TOPOLOGIES = ('no_sim', 'direct', 'shower')


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
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / normalization,
                yerr=np.sqrt(counts[shown]) / normalization,
                fmt='o', color='C0', markersize=4, capsize=2, linewidth=0.8)
    if not linear:
        ax.set_yscale('log')
        ax.set_ylim(top=counts.max() / normalization * 10)
        ax.set_xscale(xscale)
        if xscale == 'symlog':
            ax.set_xscale('symlog', linthresh=0.001)  # Display scale only, never an angle cut.
    else:
        ax.set_ylim(0, 1.4 * (counts + np.sqrt(counts)).max() / normalization)
    ax.set_xlim(bins[0], bins[-1])
    if np.issubdtype(entries.dtype, np.integer):
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    finish(fig, ax, path, sample, xlabel, ylabel)


def fraction_plot(path, sample, coordinates, bins, numerators, denominators, labels, linear):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    centers = np.sqrt(bins[:-1] * bins[1:]) if coordinates == 'energy' and not linear else (bins[:-1] + bins[1:]) / 2
    for index, (numerator, denominator, label) in enumerate(zip(numerators, denominators, labels, strict=True)):
        assert np.all(numerator <= denominator)
        shown = denominator > 0
        ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=shown)
        lower, upper = clopper_pearson(numerator, denominator)
        ax.stairs(ratio, bins, baseline=None, color=f'C{index}', label=label)
        ax.errorbar(centers[shown], ratio[shown],
                    yerr=[ratio[shown] - lower[shown], upper[shown] - ratio[shown]],
                    fmt='o', color=f'C{index}', markersize=4, capsize=2, linewidth=0.8)
    if coordinates == 'energy' and not linear:
        ax.set_xscale('log')
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.45 if len(labels) > 1 else 1.25)
    ax.set_yticks(np.linspace(0, 1, 6))
    xlabel = r'$E_\gamma^{\rm gen}$ [GeV]' if coordinates == 'energy' else r'$\cos\theta_\gamma^{\rm gen}$'
    finish(fig, ax, path, sample, xlabel, 'Fraction' if len(labels) > 1 else 'Efficiency')


def map_plot(path, sample, entries, bins, normalization, xlabel, ylabel, colorlabel, linear, energy_axis=False):
    counts = np.histogram2d(*entries, bins=bins)[0]
    assert counts.sum() == len(entries[0]), 'Map must retain every selected Gen photon'
    rates = counts / normalization
    norm = Normalize(vmin=0, vmax=rates.max()) if linear else LogNorm(vmin=rates[rates > 0].min(), vmax=rates.max())
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    mesh = ax.pcolormesh(*bins, rates.T if linear else np.ma.masked_equal(rates.T, 0),
                         cmap='viridis', norm=norm, shading='flat', rasterized=True)
    fig.colorbar(mesh, ax=ax, pad=0.025).set_label(colorlabel)
    if energy_axis and not linear:
        ax.set_xscale('log')
    if not energy_axis:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlim(bins[0][0], bins[0][-1])
    ax.set_ylim(bins[1][0], bins[1][-1])
    finish(fig, ax, path, sample, xlabel, ylabel)


def plot_population(output, sample, pop, all_populations, n_events):
    """Seventeen requested plots, plus fully linear versions of logarithmic ones."""
    energies = np.concatenate([other['energy'] for other in all_populations])
    log_energy = log_bins(energies, 6)
    efficiency_energy = log_bins(energies, 3)
    linear_energy = np.arange(0, 10 * np.ceil(energies.max() / 10) + 1)
    with_sim, with_reco = pop['sim_count'] > 0, pop['reco_count'] > 0
    n_gen = len(pop['energy'])
    assert np.all(~with_reco | with_sim)
    assert np.isin(pop['topology'], [0, 1, 2]).all()
    for linear in (False, True):
        suffix = '_linear' if linear else ''
        e_bins = linear_energy if linear else log_energy
        count_plot(output / f'gen_energy{suffix}.png', sample, pop['energy'], e_bins, n_events,
                   r'$E_\gamma^{\rm gen}$ [GeV]', COUNT_LABEL, 'Gen', linear, 'log')
        count_plot(output / f'gen_cos_theta{suffix}.png', sample, pop['cos_theta'], COS_BINS, n_events,
                   r'$\cos\theta_\gamma^{\rm gen}$', COUNT_LABEL, 'Gen', linear)
        map_plot(output / f'gen_energy_cos_theta{suffix}.png', sample,
                 (pop['energy'], pop['cos_theta']), (e_bins, COS_BINS), n_events,
                 r'$E_\gamma^{\rm gen}$ [GeV]', r'$\cos\theta_\gamma^{\rm gen}$', COUNT_LABEL, linear, True)
        for coordinates, bins in (('energy', linear_energy if linear else efficiency_energy),
                                  ('cos_theta', np.linspace(-1, 1, 21))):
            if coordinates == 'cos_theta' and linear:
                continue  # These fraction/efficiency figures already have only linear axes.
            gen_counts = np.histogram(pop[coordinates], bins=bins)[0]
            sim_counts = np.histogram(pop[coordinates][with_sim], bins=bins)[0]
            reco_counts = np.histogram(pop[coordinates][with_reco], bins=bins)[0]
            assert gen_counts.sum() == n_gen and sim_counts.sum() == with_sim.sum() and reco_counts.sum() == with_reco.sum()
            valid = sim_counts > 0
            assert np.allclose(sim_counts[valid] / gen_counts[valid] * reco_counts[valid] / sim_counts[valid],
                               reco_counts[valid] / gen_counts[valid])
            topology = [np.histogram(pop[coordinates][pop['topology'] == code], bins=bins)[0] for code in range(3)]
            assert np.array_equal(np.sum(topology, axis=0), gen_counts)
            fraction_plot(output / f'topology_fraction_vs_{coordinates}{suffix}.png', sample,
                          coordinates, bins, topology, [gen_counts] * 3, TOPOLOGIES, linear)
            for name, numerator, denominator, label in (
                ('gen_to_sim', sim_counts, gen_counts, 'Gen → Sim'),
                ('sim_to_reco', reco_counts, sim_counts, 'Sim → Reco'),
                ('gen_to_reco', reco_counts, gen_counts, 'Gen → Reco'),
            ):
                fraction_plot(output / f'{name}_efficiency_vs_{coordinates}{suffix}.png', sample,
                              coordinates, bins, [numerator], [denominator], [label], linear)
        multiplicity_bins = []
        for field, filename, xlabel in (
            ('sim_gamma_count', 'sim_gamma_multiplicity', r'$N_\gamma^{\rm Sim}$'),
            ('reco_count', 'reco_gamma_multiplicity', r'$N_\gamma^{\rm reco}$'),
        ):
            bins = np.arange(-0.5, max(other[field].max() for other in all_populations) + 1.5)
            multiplicity_bins.append(bins)
            count_plot(output / f'{filename}{suffix}.png', sample, pop[field], bins, n_gen,
                       xlabel, 'Gen fraction', 'Gen', linear)
        map_plot(output / f'sim_vs_reco_multiplicity{suffix}.png', sample,
                 (pop['sim_gamma_count'], pop['reco_count']), multiplicity_bins, n_gen,
                 r'$N_\gamma^{\rm Sim}$', r'$N_\gamma^{\rm reco}$', 'Gen fraction', linear)
        for field, filename, xlabel in (
            ('leading_ratio', 'leading_energy_response', r'$E_{\rm leading}^{\rm reco}/E_\gamma^{\rm gen}$'),
            ('summed_ratio', 'summed_energy_response', r'$\sum E_\gamma^{\rm reco}/E_\gamma^{\rm gen}$'),
            ('opening_angle', 'angular_response', r'$\Delta\theta$ [rad]'),
        ):
            entries = pop[field][with_reco]
            all_entries = np.concatenate([other[field][other['reco_count'] > 0] for other in all_populations])
            assert np.all(np.isfinite(entries))
            if field == 'opening_angle':
                bins = np.linspace(0, np.pi, 81) if linear else np.r_[0, np.geomspace(1e-8, np.pi, 81)]
                xscale = 'symlog'
            else:
                assert np.all(entries > 0)
                bins = np.linspace(0, np.nextafter(all_entries.max(), np.inf), 81) if linear else log_bins(all_entries, 6)
                xscale = 'log'
            count_plot(output / f'{filename}{suffix}.png', sample, entries, bins, with_reco.sum(),
                       xlabel, 'Fraction / bin', 'Matched gen', linear, xscale)


def write_summary(output, values):
    summary = {}
    for sample, values_sample in values.items():
        populations = {}
        for name, pop in values_sample['populations'].items():
            n_gen = len(pop['energy'])
            n_sim, n_reco = int(np.count_nonzero(pop['sim_count'])), int(np.count_nonzero(pop['reco_count']))
            efficiencies = {}
            for label, numerator, denominator in (
                ('gen_to_sim', n_sim, n_gen), ('sim_to_reco', n_reco, n_sim), ('gen_to_reco', n_reco, n_gen),
            ):
                lower, upper = clopper_pearson(np.array([numerator]), np.array([denominator]))
                efficiencies[label] = dict(numerator=numerator, denominator=denominator,
                                          efficiency=numerator / denominator if denominator else None,
                                          lower=float(lower[0]) if denominator else None,
                                          upper=float(upper[0]) if denominator else None)
            populations[name] = dict(
                gen_photons=n_gen,
                topology=dict(zip(TOPOLOGIES, map(int, np.bincount(pop['topology'], minlength=3)), strict=True)),
                efficiency_68_percent_cp=efficiencies,
                sim_gamma_descendants=int(pop['sim_gamma_count'].sum()),
                associated_reco_candidates=int(pop['reco_count'].sum()),
                gen_with_multiple_reco=int(np.count_nonzero(pop['reco_count'] > 1)),
            )
        summary[sample] = dict(events=values_sample['events'], association=values_sample['association'], populations=populations)
    (output / 'study_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')


def plot_isr_photons(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    study = output_root / 'photon_study'
    values = {sample: read_sample(input_root, sample) for sample in SAMPLES}
    for population in POPULATIONS:
        all_populations = [values[sample]['populations'][population] for sample in SAMPLES]
        for sample in SAMPLES:
            output = study / population / sample
            output.mkdir(parents=True, exist_ok=True)
            plot_population(output, sample, values[sample]['populations'][population], all_populations, values[sample]['events'])
    write_summary(study, values)
    print(f'plots: {study}')
