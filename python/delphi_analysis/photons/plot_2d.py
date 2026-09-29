"""Energy–angle distributions and efficiencies for the photon study."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from scipy.stats import beta

from delphi_analysis.plot_utils import add_heatmap_label
from .data import SAMPLES


POPULATIONS = {
    "stable_gen": "Stable gen photons",
    "reco": "Reco photons",
    "isr": "All gen ISR photons",
    "collinear_isr": "Beam-collinear ISR",
    "noncollinear_isr": "Non-collinear ISR",
    "fsr": "FSR from hard pair",
    "others": "Other gen photons",
    "matched_gen": "Matched gen photons",
    "matched_reco": "Matched reco photons",
    "isr_matched_gen": "Matched gen ISR photons",
    "isr_matched_reco": "Reco photons matched to ISR",
}


def energy_bins(values_by_sample: dict) -> np.ndarray:
    """Cover every population, including soft beam photons and the energy tail."""
    energies = np.concatenate([
        values[kind]["energy"]
        for values in values_by_sample.values()
        for kind in POPULATIONS
    ])
    assert np.all(np.isfinite(energies) & (energies > 0)), "Log energy requires finite positive energies"
    low = 10.0 ** np.floor(np.log10(energies.min()))
    high = np.nextafter(energies.max(), np.inf)
    # Three bins per decade below 5 GeV, then one bin for the sparse energy tail.
    end = min(5.0, high)
    edges = np.geomspace(low, end, int(np.ceil(3 * np.log10(end / low))) + 1)
    return np.append(edges, high) if high > end else edges


def photon_counts(values: dict, e_bins: np.ndarray, cos_bins: np.ndarray) -> np.ndarray:
    energy = np.asarray(values["energy"])
    cosine = np.asarray(values["cos_theta"])
    assert len(energy) == len(cosine), "Photon coordinates must describe the same photons"
    counts = np.histogram2d(energy, cosine, bins=(e_bins, cos_bins))[0]
    assert counts.sum() == len(energy), "A photon falls outside the plotted bins"
    return counts


def draw_maps(output: Path, maps: dict, e_bins: np.ndarray, cos_bins: np.ndarray,
              title: str, color_label: str, norm, empty_label: str) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True,
                             layout="constrained")
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e4e4e4")
    used_axes = []
    for ax, sample in zip(axes.flat, SAMPLES):
        values = maps[sample]
        mesh = ax.pcolormesh(e_bins, cos_bins, np.ma.masked_invalid(values.T),
                             cmap=cmap, norm=norm, shading="flat", rasterized=True)
        ax.set_xscale("log")
        ax.set_xlim(e_bins[0], e_bins[-1])
        ax.set_ylim(-1, 1)
        ax.set_yticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xlabel(r"$E_\gamma$ [GeV]", fontsize=13)
        ax.set_ylabel(r"$\cos\theta_\gamma$", fontsize=13)
        ax.tick_params(which="both", direction="out", top=False, right=False)
        ax.tick_params(which="major", length=4, labelsize=10, labelbottom=True)
        ax.tick_params(which="minor", length=2)
        ax.grid(False)
        add_heatmap_label(ax, sample, fontsize=12)
        used_axes.append(ax)
    for ax in list(axes.flat)[len(SAMPLES):]:
        ax.set_visible(False)
    colorbar = fig.colorbar(mesh, ax=used_axes, shrink=0.86, pad=0.025)
    colorbar.set_label(color_label, fontsize=12)
    colorbar.ax.tick_params(labelsize=10)
    fig.suptitle(f"DELPHI Simulation · {title}\n{empty_label}", fontsize=14)
    fig.savefig(output, dpi=150)
    plt.close(fig)


def plot_maps(directories: dict[str, Path], values_by_sample: dict) -> None:
    e_bins = energy_bins(values_by_sample)
    cos_bins = np.linspace(-1, 1, 21)
    destinations = {
        "stable_gen": "gen", "isr": "gen", "collinear_isr": "gen",
        "noncollinear_isr": "gen", "fsr": "gen", "others": "gen",
        "reco": "reco",
        "matched_gen": "matching", "matched_reco": "matching",
        "isr_matched_gen": "matching", "isr_matched_reco": "matching",
    }
    counts = {
        kind: {
            sample: photon_counts(values[kind], e_bins, cos_bins)
            for sample, values in values_by_sample.items()
        }
        for kind in POPULATIONS
    }
    for kind, title in POPULATIONS.items():
        rates = {
            sample: hist / len(values_by_sample[sample]["stable_gen"]["multiplicity"])
            for sample, hist in counts[kind].items()
        }
        positive = np.concatenate([rate[rate > 0] for rate in rates.values()])
        norm = LogNorm(vmin=positive.min(), vmax=positive.max())
        maps = {sample: np.where(rate > 0, rate, np.nan) for sample, rate in rates.items()}
        draw_maps(directories[destinations[kind]] / f"{kind}_energy_cos_theta.png", maps, e_bins, cos_bins,
                  title, "Photons / event / bin", norm, "Gray: no photons")

    for name, gen_kind, matched_kind, output_dir in (
        ("photon", "stable_gen", "matched_gen", directories["efficiency"]),
        ("isr", "isr", "isr_matched_gen", directories["efficiency"]),
        ("photon_truth_linked", "stable_gen", "truth_matched_gen", directories["efficiency"] / "truth_linked"),
        ("isr_truth_linked", "isr", "isr_truth_matched_gen", directories["efficiency"] / "truth_linked"),
    ):
        efficiencies, widths, denominators = {}, {}, {}
        for sample, denominator in counts[gen_kind].items():
            numerator = photon_counts(values_by_sample[sample][matched_kind], e_bins, cos_bins)
            assert np.all(numerator <= denominator), "Matched gen photons must be a subset of gen photons"
            valid = denominator > 0
            efficiency = np.full(denominator.shape, np.nan)
            efficiency[valid] = numerator[valid] / denominator[valid]
            lower, upper = np.zeros(denominator.shape), np.ones(denominator.shape)
            alpha = 1 - 0.6826894921370859
            nonzero = valid & (numerator > 0)
            not_all = valid & (numerator < denominator)
            lower[nonzero] = beta.ppf(alpha / 2, numerator[nonzero], denominator[nonzero] - numerator[nonzero] + 1)
            upper[not_all] = beta.ppf(1 - alpha / 2, numerator[not_all] + 1, denominator[not_all] - numerator[not_all])
            efficiencies[sample] = efficiency
            widths[sample] = np.where(valid, upper - lower, np.nan)
            denominators[sample] = np.where(valid, denominator, np.nan)

        title = "Photon matching" if name.startswith("photon") else "ISR matching"
        if "truth_linked" in name:
            title += " · truth-linked"
        draw_maps(output_dir / f"{name}_efficiency_energy_cos_theta.png",
                  efficiencies, e_bins, cos_bins, title, "Matching efficiency",
                  Normalize(0, 1), "Gen coordinates · Gray: no gen photons")
        draw_maps(output_dir / f"{name}_efficiency_cp_width_energy_cos_theta.png",
                  widths, e_bins, cos_bins, f"{title} uncertainty",
                  "68.27% Clopper–Pearson interval width", Normalize(0, 1),
                  "Gen coordinates · Gray: no gen photons")
        maximum = max(hist.max() for hist in counts[gen_kind].values())
        draw_maps(output_dir / f"{name}_efficiency_denominator_energy_cos_theta.png",
                  denominators, e_bins, cos_bins, f"{title} denominator",
                  "Gen photons / bin", LogNorm(1, maximum),
                  "Gen coordinates · Gray: no gen photons")
