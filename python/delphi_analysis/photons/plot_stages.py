"""Gen–Sim–Reco stages, with unique gen photons as the efficiency unit."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from .plot_isr_photons import (
    COUNT_LABEL, FIGURE_SIZE, clopper_pearson, count_curve, finish, plot_bins,
)


def spectrum(output, name, arrays, bins, sample, n_events, objects, linear):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    peak = 0
    labels = ('Gen', 'Saved Sim roots', 'Associated reco') if objects else ('Gen', 'With Sim lineage', 'With reco')
    if objects and len(arrays) == 4:
        labels += ('All reco',)
    for index, (entries, label) in enumerate(zip(arrays, labels, strict=True)):
        # Retain beam-axis photons at phi=0, as in the existing photon study.
        if name == 'phi':
            entries = np.nan_to_num(entries, nan=0)
        if objects and index == len(arrays) - 1:
            counts = np.histogram(entries, bins=bins)[0]
            assert counts.sum() == len(entries)
            centers = np.sqrt(bins[:-1] * bins[1:]) if name == 'energy' and not linear else (bins[:-1] + bins[1:]) / 2
            shown = counts > 0
            ax.errorbar(centers[shown], counts[shown] / n_events,
                        yerr=np.sqrt(counts[shown]) / n_events, fmt='o', color='black',
                        markersize=6, capsize=2, linewidth=1, label=label, zorder=5)
        else:
            counts = count_curve(ax, entries, bins, n_events,
                                 ('0.4', 'C0', 'C3')[index], label,
                                 linewidth=3 if index == 0 else 2,
                                 linestyle=('-', '--', ':')[index])
        peak = max(peak, (counts + np.sqrt(counts)).max() / n_events)
    if not linear and name != 'multiplicity':
        ax.set_yscale('log')
        ax.set_ylim(top=peak * (200 if len(arrays) == 4 else 20))
        if name == 'energy':
            ax.set_xscale('log')
    else:
        ax.set_ylim(0, peak * 1.8)
    if name == 'multiplicity':
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlim(bins[0], bins[-1])
    xlabel = {'energy': r'$E_\gamma$ [GeV]', 'cos_theta': r'$\cos\theta_\gamma$',
              'phi': r'$\phi_\gamma$ [rad]', 'multiplicity': r'$N_\gamma$'}[name]
    ylabel = 'Event fraction' if name == 'multiplicity' else 'N / bin per event' if objects else COUNT_LABEL
    finish(ax, sample, xlabel, ylabel)
    fig.tight_layout()
    fig.savefig(output / f'{"gen_sim_reco" if objects else "gen_lineage"}_{name}{"_linear" if linear else ""}.png', dpi=150)
    plt.close(fig)


def efficiency(output, name, stage, bins, sample, linear):
    gen = stage[f'gen_{name}']
    with_sim, with_reco = stage['sim_count'] > 0, stage['reco_count'] > 0
    counts = [np.histogram(entries, bins=bins)[0] for entries in (gen, gen[with_sim], gen[with_reco])]
    n_gen, n_sim, n_reco = counts
    assert n_gen.sum() == len(gen)
    assert n_sim.sum() == np.count_nonzero(with_sim)
    assert n_reco.sum() == np.count_nonzero(with_reco)
    assert np.all(n_reco <= n_sim) and np.all(n_sim <= n_gen)
    valid = n_sim > 0
    assert np.allclose((n_sim[valid] / n_gen[valid]) * (n_reco[valid] / n_sim[valid]),
                       n_reco[valid] / n_gen[valid])
    centers = np.sqrt(bins[:-1] * bins[1:]) if name == 'energy' and not linear else (bins[:-1] + bins[1:]) / 2
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    for numerator, denominator, label, color, linestyle, marker, size in (
        (n_sim, n_gen, 'Gen → Sim', 'C0', '--', '.', 3),
        (n_reco, n_sim, 'Reco | Sim', 'C1', '-', 's', 5),
        (n_reco, n_gen, 'Gen → Reco', 'C2', ':', 'o', 3),
    ):
        shown = denominator > 0
        ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=shown)
        lower, upper = clopper_pearson(numerator, denominator)
        ax.stairs(ratio, bins, baseline=None, color=color, linestyle=linestyle, linewidth=2, label=label)
        ax.errorbar(centers[shown], ratio[shown],
                    yerr=[ratio[shown] - lower[shown], upper[shown] - ratio[shown]],
                    fmt=marker, color=color, markersize=size, capsize=2, linewidth=0.8)
    if name == 'energy' and not linear:
        ax.set_xscale('log')
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.45)
    ax.set_yticks(np.linspace(0, 1, 6))
    xlabel = r'$E_\gamma^{\rm gen}$ [GeV]' if name == 'energy' else r'$\cos\theta_\gamma^{\rm gen}$'
    finish(ax, sample, xlabel, 'Efficiency')
    fig.tight_layout()
    fig.savefig(output / f'stage_efficiency_vs_{name}{"_linear" if linear else ""}.png', dpi=150)
    plt.close(fig)


def lineage_counts(output, stage, sample):
    maximum = max(stage['sim_count'].max(), stage['reco_count'].max())
    bins = np.arange(-0.5, maximum + 1.5)
    for linear in (False, True):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        peak = 0
        for entries, color, label in ((stage['sim_count'], 'C0', 'Sim lineage'),
                                      (stage['reco_count'], 'C3', 'Reco')):
            counts = count_curve(ax, entries, bins, len(stage['gen_energy']), color, label)
            peak = max(peak, (counts + np.sqrt(counts)).max() / len(stage['gen_energy']))
        if linear:
            ax.set_ylim(0, peak * 1.8)
        else:
            ax.set_xscale('symlog', linthresh=1)
            ax.set_yscale('log')
            ax.set_ylim(top=peak * 20)
        ax.set_xlim(bins[0], bins[-1])
        if linear:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        finish(ax, sample, r'$N_{\rm objects}$ per gen photon', 'Gen fraction')
        fig.tight_layout()
        fig.savefig(output / f'lineage_multiplicity{"_linear" if linear else ""}.png', dpi=150)
        plt.close(fig)


def plot_stages(output: Path, values_by_sample: dict, group: str):
    """Compare all entries in a gen group, without energy or angular cuts."""
    stages = {sample: values['stages'][group] for sample, values in values_by_sample.items()}
    for sample, stage in stages.items():
        n_events = len(values_by_sample[sample]['stable_gen']['multiplicity'])
        destination = output / sample
        destination.mkdir(parents=True, exist_ok=True)
        assert len(stage['gen_energy']) == len(stage['sim_count']) == len(stage['reco_count'])
        assert np.all((stage['reco_count'] == 0) | (stage['sim_count'] > 0))
        for name in ('energy', 'cos_theta', 'phi', 'multiplicity'):
            with_sim, with_reco = stage['sim_count'] > 0, stage['reco_count'] > 0
            gen_arrays = ([stage['gen_multiplicity'], stage['sim_multiplicity'], stage['reco_multiplicity']]
                          if name == 'multiplicity' else
                          [stage[f'gen_{name}'], stage[f'gen_{name}'][with_sim], stage[f'gen_{name}'][with_reco]])
            for linear in ((False,) if name == 'multiplicity' else (False, True)):
                arrays_for_bins = [other[f'gen_{name}'] for other in stages.values()]
                bins = plot_bins(name, arrays_for_bins, linear)[0]
                spectrum(destination, name, gen_arrays, bins, sample, n_events, False, linear)
                if name != 'multiplicity':
                    arrays = [stage[f'{prefix}_{name}'] for prefix in ('gen', 'root_sim', 'reco')]
                    object_arrays = [other[f'{prefix}_{name}'] for other in stages.values()
                                     for prefix in ('gen', 'root_sim', 'reco')]
                    if group == 'stable':
                        arrays.append(values_by_sample[sample]['reco'][name])
                        object_arrays.extend(values['reco'][name] for values in values_by_sample.values())
                    object_bins = plot_bins(name, object_arrays, linear)[0]
                    spectrum(destination, name, arrays, object_bins, sample, n_events, True, linear)
            if name in ('energy', 'cos_theta'):
                for linear in ((False, True) if name == 'energy' else (False,)):
                    bins = plot_bins(name, [other[f'gen_{name}'] for other in stages.values()], linear)[0]
                    if name == 'energy' and not linear:
                        high = np.nextafter(max(other['gen_energy'].max() for other in stages.values()), np.inf)
                        bins = np.geomspace(bins[0], high, int(np.ceil(3 * np.log10(high / bins[0]))) + 1)
                    efficiency(destination, name, stage, bins, sample, linear)
        lineage_counts(destination, stage, sample)
