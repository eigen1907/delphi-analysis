"""Energy–angle maps for gen, reco, and matched ISR photons."""

from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from matplotlib.colors import LogNorm

from .data import SAMPLES


POPULATIONS = {
    "stable_gen": "Stable gen photons",
    "reco": "Reco photons",
    "isr": "All gen ISR photons",
    "collinear_isr": "Beam ISR",
    "noncollinear_isr": "Non-beam ISR",
    "fsr": "FSR",
    "others": "Others",
    "isr_matched_gen_truth": "Truth match",
    "isr_matched_gen_angle": "Opening angle match",
    "isr_matched_gen_hybrid": "Truth + angle fallback",
}


def energy_bins(values_by_sample: dict) -> np.ndarray:
    """Cover all stable gen and reco photons, including the soft beam tail."""
    energies = np.concatenate([
        values[kind]["energy"]
        for values in values_by_sample.values()
        for kind in ("stable_gen", "reco")
    ])
    assert np.all(np.isfinite(energies) & (energies > 0)), "Log energy requires finite positive energies"
    low = 10.0 ** np.floor(np.log10(energies.min()))
    high = np.nextafter(energies.max(), np.inf)
    end = min(5.0, high)
    edges = np.geomspace(low, end, int(np.ceil(3 * np.log10(end / low))) + 1)
    return np.append(edges, high) if high > end else edges


def photon_counts(values: dict, e_bins: np.ndarray, cos_bins: np.ndarray) -> np.ndarray:
    energy = np.asarray(values["energy"])
    cosine = np.asarray(values["cos_theta"])
    counts = np.histogram2d(energy, cosine, bins=(e_bins, cos_bins))[0]
    assert counts.sum() == len(energy), "A photon falls outside the plotted bins"
    return counts


def draw_maps(output: Path, maps: dict, e_bins: np.ndarray, cos_bins: np.ndarray, norm: LogNorm) -> None:
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e4e4e4")
    for sample in SAMPLES:
        fig, ax = plt.subplots(figsize=(11, 9), layout="constrained")
        mesh = ax.pcolormesh(e_bins, cos_bins, np.ma.masked_invalid(maps[sample].T),
                             cmap=cmap, norm=norm, shading="flat", rasterized=True)
        ax.set_xscale("log")
        ax.set_xlim(e_bins[0], e_bins[-1])
        ax.set_ylim(-1, 1)
        ax.set_yticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xlabel(r"$E_\gamma$ [GeV]")
        ax.set_ylabel(r"$\cos\theta_\gamma$")
        ax.grid(False)
        ax.text(0.96, 0.94, sample, transform=ax.transAxes, ha="right", va="top",
                color="white", bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.6})
        fig.colorbar(mesh, ax=ax, pad=0.025).set_label(r"$N_\gamma$ / bin per event")
        mh.label.exp_label(exp="DELPHI", llabel="Simulation",
                           rlabel="LEP 1 (91.2 GeV)", loc=0, ax=ax)
        destination = output.parent / sample
        destination.mkdir(parents=True, exist_ok=True)
        fig.savefig(destination / output.name, dpi=150)
        plt.close(fig)


def draw_gen_reco_maps(output: Path, values_by_sample: dict, e_bins: np.ndarray,
                       cos_bins: np.ndarray) -> None:
    rates = {
        (sample, kind): photon_counts(values[kind], e_bins, cos_bins)
        / len(values["stable_gen"]["multiplicity"])
        for sample, values in values_by_sample.items()
        for kind in ("stable_gen", "reco")
    }
    positive = np.concatenate([rate[rate > 0] for rate in rates.values()])
    norm = LogNorm(vmin=positive.min(), vmax=positive.max())
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e4e4e4")
    for sample in SAMPLES:
        fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey=True,
                                 layout="constrained")
        for ax, kind in zip(axes, ("stable_gen", "reco"), strict=True):
            rate = rates[sample, kind]
            mesh = ax.pcolormesh(e_bins, cos_bins, np.ma.masked_equal(rate.T, 0),
                                 cmap=cmap, norm=norm, shading="flat", rasterized=True)
            ax.set_xscale("log")
            ax.set_xlim(e_bins[0], e_bins[-1])
            ax.set_ylim(-1, 1)
            ax.set_yticks([-1, -0.5, 0, 0.5, 1])
            ax.set_xlabel(r"$E_\gamma$ [GeV]")
            ax.grid(False)
            label = "Gen" if kind == "stable_gen" else f"Reco · {sample}"
            ax.text(0.01, 0.52, label, transform=ax.transAxes, va="center")
        axes[0].set_ylabel(r"$\cos\theta_\gamma$")
        fig.colorbar(mesh, ax=axes, pad=0.02).set_label(r"$N_\gamma$ / bin per event")
        mh.label.exp_label(exp="DELPHI", llabel="Simulation", rlabel="", loc=0, ax=axes[0])
        axes[1].text(1, 1.02, "LEP 1 (91.2 GeV)", transform=axes[1].transAxes,
                     ha="right", va="bottom")
        destination = output.parent / sample
        destination.mkdir(parents=True, exist_ok=True)
        fig.savefig(destination / output.name, dpi=150)
        plt.close(fig)


def plot_maps(directories: dict[str, Path], values_by_sample: dict) -> None:
    e_bins = energy_bins(values_by_sample)
    cos_bins = np.linspace(-1, 1, 21)
    draw_gen_reco_maps(directories["gen_reco"] / "gen_reco_energy_cos_theta.png",
                       values_by_sample, e_bins, cos_bins)
    for kind, stage in (
        ("isr_matched_gen_truth", "matching"),
        ("isr_matched_gen_angle", "matching"),
        ("isr_matched_gen_hybrid", "matching"),
    ):
        rates = {
            sample: photon_counts(values[kind], e_bins, cos_bins) / len(values["stable_gen"]["multiplicity"])
            for sample, values in values_by_sample.items()
        }
        positive = np.concatenate([rate[rate > 0] for rate in rates.values()])
        norm = LogNorm(vmin=positive.min(), vmax=positive.max())
        maps = {sample: np.where(rate > 0, rate, np.nan) for sample, rate in rates.items()}
        draw_maps(directories[stage] / f"{kind}_energy_cos_theta.png", maps, e_bins, cos_bins, norm)
