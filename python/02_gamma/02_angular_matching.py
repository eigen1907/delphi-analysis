"""Scan independent opening-angle cones and retain truth angles as a cut reference."""

import numpy as np

from .data import ANALYSIS_SELECTIONS, ENERGY_RANGES, GEN_LABELS
from common.matching import angle_matrix
from common.plot_utils import FIGURE_SIZE, clopper_pearson, finish, plot_path

CUTS_DEG = np.array([0.25, 0.5, 1, 2, 3, 5, 10, 20])


def new_accumulator():
    """Store truth retention, per-Gen cones, and event energy from their union."""
    return dict(populations={
        population: dict(
            gen_count=0,
            angular_success=np.zeros(len(CUTS_DEG), dtype=np.int64),
            angular_energy_count=np.zeros(len(CUTS_DEG), dtype=np.int64),
            angular_energy_sum=np.zeros(len(CUTS_DEG)),
            angular_energy_square_sum=np.zeros(len(CUTS_DEG)),
        ) for population in ANALYSIS_SELECTIONS
    }, truth_reference=dict(success=0, recovered=np.zeros(len(CUTS_DEG), dtype=np.int64)),
                angular={name: dict(pending=[], blocks=[])
                for name in ("reco_count", "energy_ratio", "event_energy", "detector_matched")})


def accumulate(accumulator, event, objects):
    """Every Gen photon gets its own cone; candidates can enter several cones."""
    gen_idx, gen_p4 = event["all_gamma_idx"], event["all_gamma_p4"]
    momentum = np.linalg.norm(objects["p4"][:, :3], axis=1)
    directional = momentum > 0
    p4 = objects["p4"][directional]
    angles = angle_matrix(gen_p4[:, :3], p4[:, :3])
    cones = angles[:, :, None] <= np.deg2rad(CUTS_DEG)[None, None, :]
    multiplicity = cones.sum(axis=1)
    reco_energy = p4[None, :, 3, None]
    energy_ratio = (cones * reco_energy).sum(axis=1) / gen_p4[:, 3, None]
    # Geometry ignores truth conflicts, but saved-object energy overlaps
    # still invalidate a cone's energy. Empty cones have zero energy.
    geometric_valid = objects["geometric_energy_valid"][directional]
    valid_cones = ~(cones & ~geometric_valid[None, :, None]).any(axis=1)
    stored = accumulator["angular"]
    stored["reco_count"]["pending"].append(multiplicity.astype(np.uint16))
    stored["energy_ratio"]["pending"].append(np.where(valid_cones, energy_ratio, np.nan))
    detector_flags = objects['detector_flags'][directional]
    detector_matched = (cones[:, :, None, :] & detector_flags[None, :, :, None]).any(axis=1)
    stored['detector_matched']['pending'].append(detector_matched)
    # Event energy counts each Reco object once, even when Gen cones overlap.
    event_energy = np.zeros((1, len(ENERGY_RANGES), len(ANALYSIS_SELECTIONS), len(CUTS_DEG)))
    for range_index, (range_name, _, _) in enumerate(ENERGY_RANGES):
        for population_index, population in enumerate(ANALYSIS_SELECTIONS):
            selected = event["energy_selections"][range_name][population]
            union = cones[selected].any(axis=0)
            energy = (union * p4[:, 3, None]).sum(axis=0)
            valid = ~(union & ~geometric_valid[:, None]).any(axis=0)
            event_energy[0, range_index, population_index] = np.where(valid, energy, np.nan)
    stored["event_energy"]["pending"].append(event_energy)
    # Combine small event arrays to avoid one Python object per event.
    if len(stored["reco_count"]["pending"]) == 1024:
        for field in stored.values():
            field["blocks"].append(np.concatenate(field["pending"]))
            field["pending"].clear()

    for population in ANALYSIS_SELECTIONS:
        selected = event['analysis_selections'][population]
        values = accumulator['populations'][population]
        values['gen_count'] += int(selected.sum())
        values['angular_success'] += (multiplicity[selected] > 0).sum(axis=0)
        values['angular_energy_count'] += valid_cones[selected].sum(axis=0)
        ratios = np.where(valid_cones[selected], energy_ratio[selected], 0)
        values['angular_energy_sum'] += ratios.sum(axis=0)
        values['angular_energy_square_sum'] += (ratios ** 2).sum(axis=0)
    # Raw banks contain no Sim truth information; do not invent unresolved links.
    if 'gen_idx' not in objects:
        return
    own_origin = gen_idx[:, None] == objects["gen_idx"][directional][None, :]
    truth_cones = cones & own_origin[:, :, None]
    # Keep the established angle reference: all non-beam ISR, before fiducial selection.
    selected = event['all_gamma_category'] == 'nonbeam_isr'
    reference = accumulator['truth_reference']
    reference['success'] += int((objects['stats']['reco_count'][selected] > 0).sum())
    reference['recovered'] += truth_cones[selected].any(axis=1).sum(axis=0)


def angular_stats(accumulator, cut_degrees):
    """Extract the angle cut; event_energy axes are event, energy range, population."""
    index = int(np.flatnonzero(CUTS_DEG == cut_degrees).item())
    stats = {}
    for name, field in accumulator.pop("angular").items():
        if field["pending"]:
            field["blocks"].append(np.concatenate(field["pending"]))
        stats[name] = np.concatenate([block[..., index] for block in field["blocks"]])
    stats["energy_valid"] = np.isfinite(stats["energy_ratio"])
    return stats


def mean_and_sem(counts, sums, squares):
    """Mean including unmatched zeros; uncertainty is the standard error."""
    counts = np.broadcast_to(counts, np.shape(sums))
    mean = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    variance_sum = squares - np.divide(sums ** 2, counts, out=np.zeros_like(sums), where=counts > 0)
    sem = np.sqrt(np.divide(np.maximum(variance_sum, 0), counts * (counts - 1),
                            out=np.full_like(sums, np.nan), where=counts > 1))
    return mean, sem


def choose_cut(accumulators):
    """Smallest cut retaining 99% of pooled non-beam truth successes."""
    successes = 0
    recovered = np.zeros(len(CUTS_DEG), dtype=np.int64)
    for accumulator in accumulators:
        reference = accumulator['truth_reference']
        successes += reference['success']
        recovered += reference['recovered']
    eligible = (successes > 0) & (recovered >= 0.99 * successes)
    index = np.flatnonzero(eligible)[0] if eligible.any() else len(CUTS_DEG) - 1
    return float(CUTS_DEG[index])


def plot_angle_matching(output_root, sample, accumulator, cut_degrees):
    """All, ISR and no ISR, after excluding the beam generator component."""
    import matplotlib.pyplot as plt

    output_root.mkdir(parents=True, exist_ok=True)
    shown = CUTS_DEG <= 10
    for observable in ("efficiency", "energy_recovery"):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for gen_index, population in enumerate(ANALYSIS_SELECTIONS):
            values = accumulator["populations"][population]
            if observable == "efficiency":
                denominator = np.full(len(CUTS_DEG), values["gen_count"])
                numerator = values["angular_success"]
                mean = np.divide(numerator, denominator, out=np.full(len(CUTS_DEG), np.nan), where=denominator > 0)
                lower, upper = clopper_pearson(numerator, denominator)
                error = np.maximum([mean - lower, upper - mean], 0)
            else:
                mean, error = mean_and_sem(values["angular_energy_count"], values["angular_energy_sum"],
                                           values["angular_energy_square_sum"])
            ax.errorbar(CUTS_DEG[shown], mean[shown], yerr=error[..., shown], color=f"C{gen_index}",
                        marker="o", label=GEN_LABELS[gen_index], capsize=2, markersize=4)
        ax.axvline(cut_degrees, color="black", linestyle="--", linewidth=2)
        ax.text(0.02, 0.89, 'Gen fiducial', transform=ax.transAxes, va='top')
        ax.set_xlim(0, 10)
        if observable == "efficiency":
            ax.set_ylim(0, 1.35)
            ax.axhline(1, color='0.5', linestyle='--', linewidth=1.5)
        else:
            ax.set_ylim(0, ax.get_ylim()[1] * 1.65)
        ylabel = "Angular coverage" if observable == "efficiency" else r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$'
        finish(fig, ax, output_root / f"{observable}_vs_opening_angle.png", sample,
               "Opening-angle cut [deg]", ylabel)


def plot(study, sample, values, bins, accumulator, cut_degrees):
    from .detector_plots import detector_efficiencies
    from .results import matching_results

    stats = angular_stats(accumulator, cut_degrees)
    output = plot_path(study, '02_gamma', '02_angular_matching', sample, values['source'],
                       values['collection'], 'matching_validation')
    plot_angle_matching(output, sample, accumulator, cut_degrees)
    detector_efficiencies(study, sample, values, stats, '02_angular_matching')
    matching_results(study, sample, values, bins, stats, '02_angular_matching')
