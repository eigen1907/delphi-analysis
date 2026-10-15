"""Saved unique-origin, same-charge associations from Gen to track to RICH."""
import matplotlib.pyplot as plt
import numpy as np
from common.plot_utils import FIGURE_SIZE, finish, plot_path
from .data import LINK_LABELS, SOURCES
from .plot_helpers import opening_angles, matching_results

def matching_validation(values, output, sample):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    x = np.arange(len(LINK_LABELS))
    counts = np.asarray([values['links'][label] for label in LINK_LABELS])
    ax.bar(x, counts, yerr=np.sqrt(counts), capsize=2, color='C0')
    ax.set_xticks(x, [label.replace(' only', '\nonly') for label in LINK_LABELS])
    ax.set_ylim(0, max(counts.max(), 1) * 1.3)
    finish(fig, ax, output / 'truth_link_status.png', sample, 'Link status', 'Reco tracks')
    opening_angles(values, "truth", output, sample)


def run(values, files, source, base, sample):
    if "truth" not in values["methods"]:
        return  # Raw SDST has no Part/Sim truth links.
    validation = plot_path(base, "03_rich", "03_truth_link_matching", sample, source,
                           SOURCES[source], section="matching_validation", cut="cut_no")
    result = plot_path(base, "03_rich", "03_truth_link_matching", sample, source,
                       SOURCES[source], section="matching_result", cut="cut_no")
    matching_validation(values, validation, sample)
    matching_results(values, "truth", result, sample)


if __name__ == "__main__":
    from .plot import main
    main("03_truth_link_matching")
