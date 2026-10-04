"""Qualitative soft/collinear ISR checks, without a precision QED prediction."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .plot_2d import POPULATIONS


def plot_generator(output_root: Path, values_by_sample: dict):
    from .plot_isr_photons import COMPONENTS, FIGURE_SIZE, finish

    fractions = np.concatenate([values['isr']['x_gamma'] for values in values_by_sample.values()])
    assert np.all(np.isfinite(fractions) & (fractions > 0))
    log_bins = np.geomspace(10 ** np.floor(np.log10(fractions.min())),
                           10 ** np.ceil(np.log10(fractions.max())), 61)
    linear_bins = np.linspace(0, np.ceil(10 * fractions.max()) / 10, 51)
    cosine_bins = np.linspace(0, 1, 41)
    for name, logarithmic_bins, uniform_bins, xlabel, ylabel in (
        ('x_gamma', log_bins, linear_bins, r'$x_\gamma$', r'$dN_\gamma/dx_\gamma$ per event'),
        ('abs_cos_theta', cosine_bins, cosine_bins, r'$|\cos\theta_\gamma|$',
         r'$dN_\gamma/d|\cos\theta_\gamma|$ per event'),
    ):
        for linear in (False, True):
            bins = uniform_bins if linear else logarithmic_bins
            widths = np.diff(bins)
            centers = np.sqrt(bins[:-1] * bins[1:]) if name == 'x_gamma' and not linear else (bins[:-1] + bins[1:]) / 2
            for sample, values in values_by_sample.items():
                destination = output_root / '01_gen_reco' / sample
                destination.mkdir(parents=True, exist_ok=True)
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                n_events = len(values['stable_gen']['multiplicity'])
                bottom = np.zeros(len(widths))
                cumulative = np.zeros(len(widths))
                for kind, color, hatch in COMPONENTS[:2]:
                    data = values[kind]['x_gamma'] if name == 'x_gamma' else np.abs(values[kind]['cos_theta'])
                    counts = np.histogram(data, bins=bins)[0]
                    assert counts.sum() == len(data), 'Generator bins must include every ISR photon'
                    heights = counts / (n_events * widths)
                    ax.stairs(bottom + heights, bins, baseline=bottom, fill=True,
                              facecolor=color, edgecolor=color, hatch=hatch,
                              alpha=0.75, label=POPULATIONS[kind])
                    cumulative += counts
                    shown = counts > 0
                    ax.errorbar(centers[shown], (bottom + heights)[shown],
                                yerr=np.sqrt(cumulative[shown]) / (n_events * widths[shown]),
                                fmt='none', color=color, capsize=1, linewidth=0.8)
                    bottom += heights
                ax.set_xlim(bins[0], bins[-1])
                if linear:
                    ax.set_ylim(0, bottom.max() * 1.35)
                else:
                    ax.set_yscale('log')
                    ax.set_ylim(top=bottom.max() * 10)
                    if name == 'x_gamma':
                        ax.set_xscale('log')
                finish(ax, sample, xlabel, ylabel,
                       'upper right' if name == 'x_gamma' else 'upper left')
                fig.tight_layout()
                fig.savefig(destination / f'gen_isr_{name}{"_linear" if linear else ""}.png', dpi=150)
                plt.close(fig)
