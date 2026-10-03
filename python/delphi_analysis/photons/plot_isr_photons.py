"""Generator, reco, matching, and ISR-efficiency plots for the Florian samples."""
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from scipy.stats import beta

from .data import MATCH_METHODS, MAX_ANGLE, SCAN_MAX_ANGLE, SAMPLES, read_sample
from .plot_2d import POPULATIONS, energy_bins, plot_maps

COS_BINS = np.linspace(-1, 1, 41)
PHI_BINS = np.linspace(0, 2 * np.pi, 41)
COVERAGE = 0.6826894921370859
FIGURE_SIZE = (11, 9)
COUNT_LABEL = r'$N_\gamma$ / bin per event'
METHOD_LABELS = {
    'truth': 'Truth only',
    'angle': f'Angle only (< {MAX_ANGLE:g} rad)',
    'hybrid': f'Truth + angle (< {MAX_ANGLE:g} rad)',
}


def clopper_pearson(numerator, denominator):
    alpha = 1 - COVERAGE
    lower, upper = np.zeros_like(denominator, dtype=float), np.ones_like(denominator, dtype=float)
    nonzero, not_all = numerator > 0, numerator < denominator
    lower[nonzero] = beta.ppf(alpha / 2, numerator[nonzero], denominator[nonzero] - numerator[nonzero] + 1)
    upper[not_all] = beta.ppf(1 - alpha / 2, numerator[not_all] + 1, denominator[not_all] - numerator[not_all])
    return lower, upper


def count_curve(ax, values, bins, n_events, color, label, linewidth=1.8, markersize=3,
                linestyle='-'):
    counts = np.histogram(values, bins=bins)[0]
    assert counts.sum() == len(values), "Histogram range must include every entry"
    ax.stairs(counts / n_events, bins, color=color, linewidth=linewidth,
              linestyle=linestyle, label=label)
    centers = (bins[:-1] + bins[1:]) / 2
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / n_events,
                yerr=np.sqrt(counts[shown]) / n_events,
                fmt='.', color=color, markersize=markersize, capsize=1.5, linewidth=0.8)
    return counts


def finish(ax, sample, xlabel, ylabel, legend_loc='upper right', legend_columns=1):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    legend = ax.legend(title=sample, loc=legend_loc, ncol=legend_columns,
                       framealpha=1, edgecolor='none')
    legend.get_title().set_fontweight('bold')
    ax.grid(alpha=0.2)
    mh.label.exp_label(exp='DELPHI', llabel='Simulation',
                       rlabel='LEP 1 (91.2 GeV)', loc=0, ax=ax)


def plot_bins(name, arrays, linear=False):
    if name == 'energy':
        values = np.concatenate(arrays)
        if linear:
            return np.arange(0, 10 * np.ceil(values.max() / 10) + 1), r'$E_\gamma$ [GeV]'
        return (np.geomspace(10 ** np.floor(np.log10(values.min())),
                            10 ** np.ceil(np.log10(values.max())), 61), r'$E_\gamma$ [GeV]')
    if name == 'cos_theta':
        return COS_BINS, r'$\cos\theta_\gamma$'
    if name == 'phi':
        return PHI_BINS, r'$\phi_\gamma$ [rad]'
    maximum = max(array.max() for array in arrays)
    return np.arange(-0.5, maximum + 1.5), r'$N_\gamma$'


def set_population_scale(ax, name, peak, log_y, linear=False):
    if name == 'energy' and not linear:
        ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')
        ax.set_ylim(top=peak * (100 if name == 'energy' else 5))
    else:
        ax.set_ylim(0, peak * (1.8 if name == 'phi' else 1.3))


COMPONENTS = (
    ('noncollinear_isr', 'C0', '..'), ('collinear_isr', 'C1', '///'),
    ('fsr', 'C2', '\\\\'), ('others', 'C3', 'xx'),
)


def plot_gen_reco(output, name, values_by_sample, sample, gen_kind, linear=False):
    bins, xlabel = plot_bins(name, [v[kind][name] for v in values_by_sample.values()
                                     for kind in ('stable_gen', 'reco')], linear)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    n_events = len(values_by_sample[sample]['stable_gen']['multiplicity'])
    bottom = np.zeros(len(bins) - 1)
    for kind, color, hatch in COMPONENTS:
        if gen_kind == 'stable_gen_wo_beam' and kind == 'collinear_isr':
            continue
        values = values_by_sample[sample][kind][name]
        if name == 'phi':
            values = np.nan_to_num(values, nan=0)
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
        total_label = 'Stable gen total' if gen_kind == 'stable_gen' else 'Gen total (no Beam ISR)'
        total = count_curve(ax, values_by_sample[sample][gen_kind][name], bins,
                            n_events, 'black', total_label, linewidth=3, markersize=5,
                            linestyle='--')
        peak = total.max() / n_events
    else:
        total_values = values_by_sample[sample][gen_kind][name]
        if name == 'phi':
            total_values = np.nan_to_num(total_values, nan=0)
        assert np.allclose(bottom * n_events, np.histogram(total_values, bins=bins)[0])
        peak = bottom.max()
    reco = values_by_sample[sample]['reco'][name]
    if name == 'phi':
        reco = np.nan_to_num(reco, nan=0)
    reco_counts = np.histogram(reco, bins=bins)[0]
    assert reco_counts.sum() == len(reco)
    centers = np.sqrt(bins[:-1] * bins[1:]) if name == 'energy' and not linear else (bins[:-1] + bins[1:]) / 2
    shown = reco_counts > 0
    ax.errorbar(centers[shown], reco_counts[shown] / n_events,
                yerr=np.sqrt(reco_counts[shown]) / n_events,
                fmt='o', color='black', markersize=6, capsize=2, linewidth=1,
                zorder=5, label=POPULATIONS['reco'])
    peak = max(peak, reco_counts.max() / n_events)
    log_y = not linear and name != 'multiplicity'
    set_population_scale(ax, name, peak, log_y, linear)
    if name == 'energy' and log_y:
        ax.set_ylim(top=peak * 50)
    if name == 'multiplicity':
        ax.set_ylim(0, 1.05)
    finish(ax, sample, xlabel, 'Event fraction' if name == 'multiplicity' else COUNT_LABEL,
           'upper right' if name == 'energy' and linear else 'upper left' if name == 'energy'
           else 'upper center' if name == 'cos_theta' else 'upper right',
           2 if name == 'energy' else 1)
    fig.tight_layout()
    fig.savefig(output / f'gen_reco_{name}{"_linear" if linear else ""}.png', dpi=150)
    plt.close(fig)


def plot_matched_isr(output, name, values_by_sample, sample, linear=False):
    keys = [f'isr_matched_gen_{method}' for method in MATCH_METHODS]
    bins, xlabel = plot_bins(name, [v[key][name] for v in values_by_sample.values() for key in keys], linear)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    n_events = len(values_by_sample[sample]['stable_gen']['multiplicity'])
    peak = 0
    for index, method in enumerate(MATCH_METHODS):
        values = values_by_sample[sample][f'isr_matched_gen_{method}'][name]
        if name == 'phi':
            values = np.nan_to_num(values, nan=0)
        counts = count_curve(ax, values, bins, n_events, f'C{index}', METHOD_LABELS[method])
        peak = max(peak, counts.max() / n_events)
    set_population_scale(ax, name, peak, name == 'energy' and not linear, linear)
    finish(ax, sample, xlabel, 'Event fraction' if name == 'multiplicity' else COUNT_LABEL)
    fig.tight_layout()
    fig.savefig(output / f'isr_matched_gen_{name}{"_linear" if linear else ""}.png', dpi=150)
    plt.close(fig)


def plot_efficiency(output, name, gen_kind, matched_prefix, filename, values, sample, e_bins, linear=False):
    bins = e_bins if name == 'energy' else np.linspace(-1, 1, 21)
    denominator = np.histogram(values[gen_kind][name], bins=bins)[0]
    assert denominator.sum() == len(values[gen_kind][name])
    centers = np.sqrt(bins[:-1] * bins[1:]) if name == 'energy' and not linear else (bins[:-1] + bins[1:]) / 2
    valid = denominator > 0
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for index, method in enumerate(MATCH_METHODS):
        selected = values[f'{matched_prefix}_{method}'][name]
        numerator = np.histogram(selected, bins=bins)[0]
        assert numerator.sum() == len(selected)
        assert np.all(numerator <= denominator)
        ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=valid)
        lower, upper = clopper_pearson(numerator, denominator)
        color = f'C{index}'
        ax.stairs(ratio, bins, color=color, linewidth=1.8, label=METHOD_LABELS[method])
        ax.errorbar(centers[valid], ratio[valid],
                    yerr=[ratio[valid] - lower[valid], upper[valid] - ratio[valid]],
                    fmt='.', color=color, markersize=4, capsize=2, linewidth=0.8)
    if name == 'energy' and not linear:
        ax.set_xscale('symlog', linthresh=0.1)
    ax.set_xlim(bins[0], bins[-1])
    if linear:
        ax.set_ylim(0, 1.05)
    else:
        ax.set_ylim(0, 0.30 if name == 'energy' else 0.08)
    xlabel = r'$E_\gamma^{\mathrm{gen}}$ [GeV]' if name == 'energy' else r'$\cos\theta_\gamma^{\mathrm{gen}}$'
    finish(ax, sample, xlabel, 'Efficiency', 'upper left')
    fig.tight_layout()
    fig.savefig(output / f'{filename}_efficiency_vs_{name}{"_linear" if linear else ""}.png', dpi=150)
    plt.close(fig)


def plot_angle_scan(output, values_by_sample, linear=False):
    cuts = np.linspace(0.001, SCAN_MAX_ANGLE, 100)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for name, color, label in (
        ('truth', 'C0', 'Truth'),
        ('no_link', 'C1', 'No gen link'),
        ('wrong_link', 'C2', 'Other gen'),
        ('beam_no_link', 'C3', 'Beam subset'),
    ):
        angles = np.sort(np.concatenate([v['angle_scan'][name] for v in values_by_sample.values()]))
        ax.plot(cuts, np.searchsorted(angles, cuts), color=color, label=label)
    ax.axvline(MAX_ANGLE, color='black', linestyle='--', label=f'{MAX_ANGLE:g} rad cut')
    if not linear:
        ax.set_yscale('log')
    else:
        ax.set_ylim(bottom=0)
    ax.set_xlim(0, SCAN_MAX_ANGLE)
    finish(ax, '5 samples', 'Opening-angle cut [rad]', 'ISR pairs',
           'center' if linear else 'lower left')
    fig.tight_layout()
    fig.savefig(output / f'isr_opening_angle_cut_scan{"_linear" if linear else ""}.png', dpi=150)
    plt.close(fig)


def plot_isr_photons(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    study = output_root / 'isr_photons'
    directories = {name: study / folder for name, folder in (
        ('gen_reco', '01_gen_reco'), ('wo_beam', '02_gen_reco_wo_beamISR'),
        ('matching', '03_matching'), ('efficiency', '04_efficiency'),
    )}
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    values = {sample: read_sample(input_root, sample) for sample in SAMPLES}
    e_bins = energy_bins(values)
    linear_e_bins = plot_bins('energy', [v['stable_gen']['energy'] for v in values.values()], True)[0]
    for sample in SAMPLES:
        sample_dirs = {name: directory / sample for name, directory in directories.items()}
        for directory in sample_dirs.values():
            directory.mkdir(parents=True, exist_ok=True)
        for name in ('multiplicity', 'energy', 'cos_theta', 'phi'):
            plot_gen_reco(sample_dirs['gen_reco'], name, values, sample, 'stable_gen')
            plot_gen_reco(sample_dirs['wo_beam'], name, values, sample, 'stable_gen_wo_beam')
            plot_matched_isr(sample_dirs['matching'], name, values, sample)
            if name != 'multiplicity':
                plot_gen_reco(sample_dirs['gen_reco'], name, values, sample, 'stable_gen', True)
                plot_gen_reco(sample_dirs['wo_beam'], name, values, sample, 'stable_gen_wo_beam', True)
            if name == 'energy':
                plot_matched_isr(sample_dirs['matching'], name, values, sample, True)
        for name in ('energy', 'cos_theta'):
            plot_efficiency(sample_dirs['efficiency'], name, 'isr', 'isr_matched_gen',
                            'all_isr', values[sample], sample, e_bins)
            plot_efficiency(sample_dirs['efficiency'], name, 'noncollinear_isr',
                            'nonbeam_isr_matched_gen', 'nonbeam_isr', values[sample], sample, e_bins)
            if name == 'energy':
                plot_efficiency(sample_dirs['efficiency'], name, 'isr', 'isr_matched_gen',
                                'all_isr', values[sample], sample, linear_e_bins, True)
                plot_efficiency(sample_dirs['efficiency'], name, 'noncollinear_isr',
                                'nonbeam_isr_matched_gen', 'nonbeam_isr', values[sample], sample, linear_e_bins, True)
    plot_angle_scan(directories['matching'], values)
    plot_angle_scan(directories['matching'], values, True)
    plot_maps(directories, values)
    print(f'plots: {study}')
