"""Charged Gen and reconstructed-track distributions."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from common.plot_utils import FIGURE_SIZE, finish, plot_path
from ..track import sample_distribution as track_plots
from .data import SOURCES
from .plot_helpers import EDGES, LABELS, CAPTION, histogram, momentum_axis

def distributions(values, output, sample):
    for name, edges in EDGES.items():
        for logarithmic in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            gen, reco = values["gen"][name], values["reco"][name]
            if name == "p":
                gen, reco = np.minimum(gen, edges[-1]), np.minimum(reco, edges[-1])
            gen_height = track_plots.gen_spectrum(ax, gen, values["gen"]["pdg"], edges, values["counts"]["events"])
            reco_height = track_plots.reco_spectrum(ax, reco, edges, values["counts"]["events"], "Reco track")
            if name == "p":
                momentum_axis(ax)
            if logarithmic:
                ax.set_yscale("log")
                height = max(gen_height.max(initial=0), reco_height.max(initial=0))
                ax.set_ylim(0.5 / values["counts"]["events"], max(height, 1 / values["counts"]["events"]) * 100)
            else:
                ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
            suffix = "_log" if logarithmic else ""
            finish(fig, ax, output / f"{name}{suffix}.png", sample, LABELS[name], "Tracks / bin / event")
    bins = np.arange(-0.5, max(values["multiplicity"]["gen"].max(), values["multiplicity"]["reco"].max()) + 1.5)
    for logarithmic in (False, True):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for index, (key, label) in enumerate((("gen", "Gen charged"), ("reco", "Reco track"),
                                             ("gas", "Gas: saved angle"), ("liquid", "Liquid: saved angle"),
                                             ("either", "Gas or Liquid: saved angle"))):
            histogram(ax, values["multiplicity"][key], bins, label, f"C{index}", 1)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if logarithmic:
            ax.set_yscale("log")
            ax.set_ylim(0.5, values["counts"]["events"] * 100)
        else:
            ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
        ax.text(0.02, 0.86, CAPTION, transform=ax.transAxes, fontsize=18)
        suffix = "_log" if logarithmic else ""
        finish(fig, ax, output / f"multiplicity{suffix}.png", sample, r"$N_{\mathrm{track}}$", "Events")


def run(values, files, source, base, sample):
    output = plot_path(base, "03_rich", "01_sample_distribution", sample, source, SOURCES[source])
    distributions(values, output, sample)


if __name__ == "__main__":
    from .plot import main
    main("01_sample_distribution")
