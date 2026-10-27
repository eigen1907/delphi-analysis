"""Validate both truth-link directions and draw the nominal truth result."""
import matplotlib.pyplot as plt
import numpy as np

from common.plot_utils import FIGURE_SIZE, finish, plot_path
from .data import ANALYSIS_SELECTIONS, GEN_LABELS, LINK_STATUSES
from .detector_plots import detector_efficiencies
from .results import matching_results, selection_views

def truth_link_validation(study, sample, values):
    """Count candidates with known selected origins; No origin is shown in stage 00."""
    labels = ('Agree', 'Forward\nonly', 'Reverse\nonly', 'Conflict', 'No origin')
    scope = 'fiducial' if values['fiducial'] else 'all'
    output = plot_path(study, '02_gamma', '03_truth_link_matching', sample, values['source'],
                       values['collection'], f'matching_validation/{scope}')
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    x = np.arange(len(labels))
    maximum = 0
    for index, (population, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
        counts = np.array([values['matching_validation'][population][key] for key in LINK_STATUSES])
        positions = x + (index - 1) * 0.25
        ax.bar(positions, counts, width=0.25, color=f'C{index}', label=label)
        ax.errorbar(positions, counts, yerr=np.sqrt(counts), fmt='none', color=f'C{index}', capsize=2)
        maximum = max(maximum, counts.max())
    ax.set_xticks(x, labels)
    ax.set_yscale('log')
    ax.set_ylim(0.5, max(maximum * 100, 1))
    if values['fiducial']:
        ax.text(0.02, 0.89, 'Gen fiducial', transform=ax.transAxes, va='top', fontsize=20)
    finish(fig, ax, output / 'truth_link_status.png', sample,
           'Link status', f"# of linked {values['reco_label']}")


def plot(study, sample, values, bins):
    for selected in selection_views(values):
        truth_link_validation(study, sample, selected)
        detector_efficiencies(study, sample, selected, values['reco'], '03_truth_link_matching')
        matching_results(study, sample, selected, bins, values['reco'], '03_truth_link_matching')
