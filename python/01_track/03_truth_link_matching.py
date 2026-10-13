"""Validate saved Gen origins and plot same-charge truth-associated tracks."""
import numpy as np
import matplotlib.pyplot as plt

from common.plot_utils import FIGURE_SIZE, finish, plot_path
from .data import LINK_LABELS
from .results import matching_result


def plot(root, dataset, source, collection, values):
    output = plot_path(root, '01_track', '03_truth_link_matching', dataset, source, collection,
                       section='matching_validation', cut='cut_no')
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    counts = np.array([values['links'][label] for label in LINK_LABELS])
    ax.bar(np.arange(len(counts)), counts, color='C0')
    ax.errorbar(np.arange(len(counts)), counts, yerr=np.sqrt(counts), fmt='none', color='black', capsize=3)
    ax.set_xticks(np.arange(len(counts)), [label.replace(' ', '\n') for label in LINK_LABELS])
    ax.set_ylim(0, counts.max(initial=0) * 1.7)
    finish(fig, ax, output / 'truth_link_status.png', dataset, 'Link status', 'Usable tracks')
    matching_result(root, dataset, source, collection, values, 'truth')
