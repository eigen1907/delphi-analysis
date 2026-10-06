"""Scan truth-associated opening angles, then compare independent cone matching."""

import numpy as np

from .data import ANALYSIS_SELECTIONS, GEN_LABELS

CUTS_DEG = np.array([0.25, 0.5, 1, 2, 3, 5, 10, 20])


def new_accumulator():
    """Store truth-retention curves and per-Gen cone counts/energy ratios."""
    return dict(populations={
        population: dict(
            gen_count=0, truth_success=0,
            recovered_truth_success=np.zeros(len(CUTS_DEG), dtype=np.int64),
            truth_energy_count=np.zeros(len(CUTS_DEG), dtype=np.int64),
            truth_energy_sum=np.zeros(len(CUTS_DEG)),
            truth_energy_square_sum=np.zeros(len(CUTS_DEG)),
        ) for population in ANALYSIS_SELECTIONS
    }, angular={name: dict(pending=[], blocks=[]) for name in ("reco_count", "energy_ratio")})


def accumulate(accumulator, event):
    """Every Gen photon gets its own cone; candidates can enter several cones."""
    gen_idx, gen_p4 = event["all_gamma_idx"], event["all_gamma_p4"]
    gen_unit = gen_p4[:, :3] / np.linalg.norm(gen_p4[:, :3], axis=1)[:, None]
    cut_cosines = np.cos(np.deg2rad(CUTS_DEG))
    objects = event["reco"]
    momentum = np.linalg.norm(objects["p4"][:, :3], axis=1)
    directional = momentum > 0
    p4 = objects["p4"][directional]
    reco_unit = p4[:, :3] / momentum[directional, None]
    cosine = np.clip(gen_unit @ reco_unit.T, -1, 1)
    cones = cosine[:, :, None] >= cut_cosines[None, None, :]
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
    # Combine small event arrays to avoid one Python object per event.
    if len(stored["reco_count"]["pending"]) == 1024:
        for field in stored.values():
            field["blocks"].append(np.concatenate(field["pending"]))
            field["pending"].clear()

    own_origin = gen_idx[:, None] == objects["gen_idx"][directional][None, :]
    truth_cones = cones & own_origin[:, :, None]
    truth_ratio = (truth_cones * reco_energy).sum(axis=1) / gen_p4[:, 3, None]
    energy_valid = objects["energy_valid"][directional]
    truth_valid = ~(truth_cones & ~energy_valid[None, :, None]).any(axis=1)
    truth_valid &= objects["stats"]["energy_valid"][:, None]
    for population in ANALYSIS_SELECTIONS:
        selected = event["analysis_selections"][population]
        values = accumulator["populations"][population]
        values["gen_count"] += int(selected.sum())
        values["truth_success"] += int((objects["stats"]["reco_count"][selected] > 0).sum())
        values["recovered_truth_success"] += truth_cones[selected].any(axis=1).sum(axis=0)
        values["truth_energy_count"] += truth_valid[selected].sum(axis=0)
        ratios = np.where(truth_valid[selected], truth_ratio[selected], 0)
        values["truth_energy_sum"] += ratios.sum(axis=0)
        values["truth_energy_square_sum"] += (ratios ** 2).sum(axis=0)


def angular_stats(accumulator, cut_degrees):
    """Extract the chosen cone cut, with no truth gate or unique assignment."""
    index = int(np.flatnonzero(CUTS_DEG == cut_degrees).item())
    stats = {}
    for name, field in accumulator.pop("angular").items():
        if field["pending"]:
            field["blocks"].append(np.concatenate(field["pending"]))
        stats[name] = np.concatenate([block[:, index] for block in field["blocks"]])
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
        values = accumulator["populations"]["gen_isr"]
        successes += values["truth_success"]
        recovered += values["recovered_truth_success"]
    eligible = (successes > 0) & (recovered >= 0.99 * successes)
    index = np.flatnonzero(eligible)[0] if eligible.any() else len(CUTS_DEG) - 1
    return float(CUTS_DEG[index])


def plot_angle_matching(output_root, sample, accumulator, cut_degrees):
    """All, ISR and no ISR, after excluding the beam generator component."""
    import matplotlib.pyplot as plt

    from .plot_isr_study import FIGURE_SIZE, clopper_pearson, finish

    output_root.mkdir(parents=True, exist_ok=True)
    shown = CUTS_DEG <= 10
    for observable in ("efficiency", "energy_recovery"):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for gen_index, population in enumerate(ANALYSIS_SELECTIONS):
            values = accumulator["populations"][population]
            if observable == "efficiency":
                denominator = np.full(len(CUTS_DEG), values["gen_count"])
                numerator = values["recovered_truth_success"]
                mean = np.divide(numerator, denominator, out=np.full(len(CUTS_DEG), np.nan), where=denominator > 0)
                lower, upper = clopper_pearson(numerator, denominator)
                error = np.maximum([mean - lower, upper - mean], 0)
            else:
                mean, error = mean_and_sem(values["truth_energy_count"], values["truth_energy_sum"],
                                           values["truth_energy_square_sum"])
            ax.errorbar(CUTS_DEG[shown], mean[shown], yerr=error[..., shown], color=f"C{gen_index}",
                        marker="o", label=GEN_LABELS[gen_index], capsize=2, markersize=4)
        ax.axvline(cut_degrees, color="black", linestyle="--", linewidth=2)
        ax.set_xlim(0, 10)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.65)
        ylabel = "Efficiency" if observable == "efficiency" else r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$'
        finish(fig, ax, output_root / f"{observable}_vs_opening_angle_{sample}.png", sample,
               "Opening-angle cut [deg]", ylabel)
