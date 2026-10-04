"""ISR response from reco photons linked through the simulation ancestry."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import MaxNLocator

from .plot_isr_photons import FIGURE_SIZE, count_curve, finish


def plot_count(output, entries, bins, linear_bins, n_events, sample, xlabel, ylabel,
               label, symlog=False):
    for edges, suffix, linear in ((bins, '', False), (linear_bins, '_linear', True)):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        counts = count_curve(ax, entries, edges, n_events, 'C0', label)
        ax.set_xlim(edges[0], edges[-1])
        if np.issubdtype(entries.dtype, np.integer):
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if linear:
            ax.set_ylim(0, 1.3 * (counts + np.sqrt(counts)).max() / n_events)
        else:
            ax.set_yscale('log')
            ax.set_ylim(top=5 * counts.max() / n_events)
            if symlog:
                ax.set_xscale('symlog', linthresh=0.01)
        finish(ax, sample, xlabel, ylabel)
        fig.tight_layout()
        fig.savefig(output.parent / f'{output.stem}{suffix}.png', dpi=150)
        plt.close(fig)


def residual_bins(residual):
    """Keep the complete residual range, including large positive tails."""
    negative = np.geomspace(0.01, max(1, -residual.min()), 41)
    positive = np.geomspace(0.01, max(1, residual.max()), 61)
    return np.concatenate((-negative[::-1], [0], positive))


def plot_response(output: Path, values: dict, sample: str):
    """Use every truth-associated pair, without energy or opening-angle cuts."""
    response = values['response']
    n_events = len(values['stable_gen']['multiplicity'])
    gen, reco = response['gen_energy'], response['reco_energy']
    assert len(gen) == len(reco) == len(response['opening_angle']) == len(response['depth'])
    assert np.all(np.isfinite(gen) & (gen > 0))
    assert np.all(np.isfinite(reco) & (reco > 0))

    for linear in (False, True):
        if linear:
            gen_bins = np.arange(0, 10 * np.ceil(gen.max() / 10) + 1)
            reco_bins = np.arange(0, 10 * np.ceil(reco.max() / 10) + 1)
        else:
            gen_bins = np.geomspace(10 ** np.floor(np.log10(gen.min())),
                                   10 ** np.ceil(np.log10(gen.max())), 61)
            reco_bins = np.geomspace(10 ** np.floor(np.log10(reco.min())),
                                    10 ** np.ceil(np.log10(reco.max())), 61)
        counts = np.histogram2d(gen, reco, bins=(gen_bins, reco_bins))[0]
        assert counts.sum() == len(gen), 'Energy map must include every associated reco photon'
        rates = counts / n_events
        norm = Normalize(vmin=0, vmax=rates.max()) if linear else LogNorm(
            vmin=rates[rates > 0].min(), vmax=rates.max())
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        mesh = ax.pcolormesh(gen_bins, reco_bins,
                             rates.T if linear else np.ma.masked_equal(rates.T, 0),
                             norm=norm, cmap='viridis', shading='flat', rasterized=True)
        upper = max(gen_bins[-1], reco_bins[-1])
        ax.plot([0, upper], [0, upper], color='black', linestyle='--',
                label=r'$E_\gamma^{\rm reco}=E_\gamma^{\rm gen}$')
        if not linear:
            ax.set_xscale('log')
            ax.set_yscale('log')
        ax.set_xlim(gen_bins[0], gen_bins[-1])
        ax.set_ylim(reco_bins[0], reco_bins[-1])
        fig.colorbar(mesh, ax=ax, pad=0.025).set_label(r'$N_{\rm reco}$ / bin per event')
        finish(ax, sample, r'$E_\gamma^{\rm gen}$ [GeV]', r'$E_\gamma^{\rm reco}$ [GeV]', 'upper left')
        ax.get_legend().set_frame_on(True)
        ax.get_legend().get_frame().set_facecolor('white')
        ax.grid(False)
        fig.tight_layout()
        fig.savefig(output / f'isr_energy_response{"_linear" if linear else ""}.png', dpi=150)
        plt.close(fig)

    for name, gen_energy, reco_energy, label, ylabel, xlabel in (
        ('isr_energy_residual', gen, reco, 'Truth ancestry', r'$N_{\rm reco}$ / bin per event',
         r'$(E_\gamma^{\rm reco}-E_\gamma^{\rm gen})/E_\gamma^{\rm gen}$'),
        ('isr_summed_energy_residual', response['group_gen_energy'], response['group_reco_energy'],
         'Sum of linked reco photons', r'$N_{\rm ISR}$ / bin per event',
         r'$(\sum E_\gamma^{\rm reco}-E_\gamma^{\rm gen})/E_\gamma^{\rm gen}$'),
    ):
        residual = (reco_energy - gen_energy) / gen_energy
        assert np.all(np.isfinite(residual))
        plot_count(output / f'{name}.png', residual, residual_bins(residual),
                   np.linspace(min(0, residual.min()), max(0, residual.max()), 81),
                   n_events, sample, xlabel, ylabel, label, symlog=True)

    angle_bins = np.linspace(0, np.nextafter(response['opening_angle'].max(), np.inf), 81)
    plot_count(output / 'isr_opening_angle.png', response['opening_angle'], angle_bins, angle_bins,
               n_events, sample, r'$\Delta\theta$ [rad]', r'$N_{\rm reco}$ / bin per event', 'Truth ancestry')
    for name, entries, xlabel, ylabel, label in (
        ('isr_reco_descendants', response['reco_count'], r'$N_{\rm reco}$ per gen ISR',
         r'$N_{\rm ISR}$ / bin per event', 'All gen ISR'),
        ('isr_sim_ancestry_depth', response['depth'], 'Sim parent steps',
         r'$N_{\rm reco}$ / bin per event', 'Truth ancestry'),
    ):
        bins = np.arange(-0.5, entries.max() + 1.5)
        plot_count(output / f'{name}.png', entries, bins, bins,
                   n_events, sample, xlabel, ylabel, label)
