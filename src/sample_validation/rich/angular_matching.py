"""Same-charge track matching; Cherenkov angle is never a direction match."""
import matplotlib.pyplot as plt
import numpy as np
from common.plot_utils import FIGURE_SIZE, clopper_pearson, finish, plot_path
from .data import ANGLE_CUT, SCAN_CUTS, SOURCES
from .plot_helpers import CAPTION, opening_angles, matching_results

def matching_validation(values, output, sample):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    cuts = np.rad2deg(SCAN_CUTS)
    for index, (name, label) in enumerate((('track', 'Track'), ('gas', 'Track + saved Gas'),
                                         ('liquid', 'Track + saved Liquid'), ('either', 'Track + saved RICH'))):
        entries = values['scan'][name]
        numerator = np.asarray([(entries < cut).sum() for cut in SCAN_CUTS])
        denominator = np.full(len(cuts), len(entries))
        lower, upper = clopper_pearson(numerator, denominator)
        efficiency = numerator / max(len(entries), 1)
        ax.plot(cuts, efficiency, color=f'C{index}', label=label)
        ax.fill_between(cuts, lower, upper, color=f'C{index}', alpha=0.15)
    ax.axvline(np.rad2deg(ANGLE_CUT), color='black', linestyle='--', linewidth=2)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1.35)
    ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
    ax.text(0.02, 0.86, CAPTION, transform=ax.transAxes, fontsize=18)
    finish(fig, ax, output / 'opening_angle_scan.png', sample, 'Opening-angle cut [deg]', 'Gen coverage')
    opening_angles(values, "angular", output, sample)


def run(values, files, source, base, sample):
    validation = plot_path(base, "03_rich", "02_angular_matching", sample, source,
                           SOURCES[source], section="matching_validation", cut="cut_no")
    result = plot_path(base, "03_rich", "02_angular_matching", sample, source,
                       SOURCES[source], section="matching_result", cut="cut_no")
    matching_validation(values, validation, sample)
    matching_results(values, "angular", result, sample)


if __name__ == "__main__":
    from .plot import main
    main("02_angular_matching")
