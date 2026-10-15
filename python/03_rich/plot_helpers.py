"""Shared plot operations for the two RICH association methods."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from common.plot_utils import FIGURE_SIZE, finish, plot_efficiency

EDGES = {"p": np.array([0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 60, 70]),
         "cos_theta": np.linspace(-1, 1, 101)}
LABELS = {"p": r"$p$ [GeV]", "cos_theta": r"$\cos\theta$"}
CAPTION = r"Saved $\theta_C > 0$; not PID efficiency"


def histogram(ax, entries, edges, label, color, denominator):
    counts = np.histogram(entries, bins=edges)[0]
    assert counts.sum() == len(entries), "Every entry must appear, including the labelled overflow bin"
    centers = (edges[:-1] + edges[1:]) / 2
    ax.stairs(counts / denominator, edges, label=label, color=color)
    ax.errorbar(centers, counts / denominator, yerr=np.sqrt(counts) / denominator,
                fmt="none", color=color, capsize=2)

def momentum_axis(ax):
    ax.set_xticks([0, 10, 20, 30, 40, 50, 65], ["0", "10", "20", "30", "40", "50", r"$\geq60$"])
    ax.set_xlim(0, 70)

def opening_angles(values, method, output, sample):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    angles = np.rad2deg(values['methods'][method]['angle'])
    edges = np.append(np.linspace(0, 10, 101), 12)
    histogram(ax, np.minimum(angles, edges[-1]), edges, 'Associated track', 'C0', 1)
    ax.set_xticks([0, 2, 4, 6, 8, 11], ['0', '2', '4', '6', '8', '$\\geq10$'])
    ax.set_xlim(0, 12)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
    finish(fig, ax, output / 'opening_angle.png', sample, 'Opening angle [deg]', 'Track pairs / bin')

def matching_results(values, method, output, sample):
    result = values["methods"][method]
    for name, edges in EDGES.items():
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        entries = values["gen"][name]
        if name == "p":
            entries = np.minimum(entries, edges[-1])
        for index, (key, label) in enumerate((("matched", "Track"), ("gas", "Track + saved Gas"),
                                             ("liquid", "Track + saved Liquid"), ("either", "Track + saved RICH"))):
            plot_efficiency(ax, entries, result[key], edges, label, f"C{index}")
        if name == "p":
            momentum_axis(ax)
        ax.set_ylim(0, 1.35)
        ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
        ax.text(0.02, 0.86, CAPTION, transform=ax.transAxes, fontsize=18)
        finish(fig, ax, output / f"track_plus_saved_angle_vs_gen_{name}.png", sample,
               "Gen " + LABELS[name], "Gen coverage")
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        entries = result[f"track_{name}"]
        if name == "p":
            entries = np.minimum(entries, edges[-1])
        for index, (key, label) in enumerate((("gas", "Gas"), ("liquid", "Liquid"), ("either", "Gas or Liquid"))):
            plot_efficiency(ax, entries, result[f"track_{key}"], edges, label, f"C{index + 1}")
        if name == "p":
            momentum_axis(ax)
        ax.set_ylim(0, 1.35)
        ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
        ax.text(0.02, 0.86, CAPTION, transform=ax.transAxes, fontsize=18)
        finish(fig, ax, output / f"saved_angle_given_track_vs_reco_{name}.png", sample,
               "Reco " + LABELS[name], "Saved angle / matched track")
    bins = np.arange(-0.5, result["reco_count"].max(initial=0) + 1.5)
    for logarithmic in (False, True):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for index, (key, label) in enumerate((("reco_count", "Associated tracks"), ("saved_count", "With saved RICH"))):
            histogram(ax, result[key], bins, label, f"C{index}", 1)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if logarithmic:
            ax.set_yscale("log")
            ax.set_ylim(0.5, len(result["reco_count"]) * 100)
        else:
            ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
        ax.text(0.02, 0.86, CAPTION, transform=ax.transAxes, fontsize=18)
        suffix = "_log" if logarithmic else ""
        finish(fig, ax, output / f"associated_track_multiplicity{suffix}.png", sample,
               r"$N_{\mathrm{track}}$ / Gen", "Gen particles / bin")
    print(f"{sample} {method}: Gen={len(result['matched'])}, with_track={result['matched'].sum()}, "
          f"with_saved_gas={result['gas'].sum()}, with_saved_liquid={result['liquid'].sum()}, "
          f"with_saved_either={result['either'].sum()}, split_Gen={(result['reco_count'] > 1).sum()}, "
          f"matched_tracks={len(result['track_p'])}, saved_gas_tracks={result['track_gas'].sum()}, "
          f"saved_liquid_tracks={result['track_liquid'].sum()}, saved_either_tracks={result['track_either'].sum()}")
