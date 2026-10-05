"""Opening-angle cross-check; stored truth remains the nominal association.

Each Gen ISR uses every candidate inside its cone. A candidate can therefore
appear in several cones; pair fractions are diagnostics, not a new truth match.
"""

import numpy as np

CUTS_DEG = np.array([0.25, 0.5, 1, 2, 3, 5, 10, 20])
CHANNELS = ("photon_conversion", "all_lineage")
POPULATIONS = ("all_ISR", "noBeamISR")
PAIR_CATEGORIES = ("own_ISR", "other_ISR", "non_ISR", "unresolved")


def new_accumulator():
    """Only scan counts and moments are retained, rather than every event."""
    result = dict(events=0, channels={})
    for channel in CHANNELS:
        populations = {}
        for population in POPULATIONS:
            values = dict(gen_count=0, truth_success=0,
                          pair_counts=np.zeros((4, len(CUTS_DEG)), dtype=np.int64))
            for field in ("geometric_success", "recovered_truth_success", "multiple_candidates",
                          "energy_count", "truth_energy_count"):
                values[field] = np.zeros(len(CUTS_DEG), dtype=np.int64)
            for field in ("energy_sum", "energy_square_sum", "truth_energy_sum",
                          "truth_energy_square_sum", "multiplicity_sum", "multiplicity_square_sum",
                          "truth_multiplicity_sum", "truth_multiplicity_square_sum"):
                values[field] = np.zeros(len(CUTS_DEG))
            populations[population] = values
        result["channels"][channel] = dict(candidates=0, zero_momentum_candidates=0,
                                            energy_invalid_candidates=0, populations=populations)
    return result


def accumulate(accumulator, event):
    """Accumulate all-ISR and non-beam ISR scans from one analyzed event."""
    accumulator["events"] += 1
    gen_idx, gen_p4 = event["gen_idx"], event["gen_p4"]
    gen_momentum = np.linalg.norm(gen_p4[:, :3], axis=1)
    assert np.all(gen_momentum > 0) and np.all(gen_p4[:, 3] > 0)
    gen_unit = gen_p4[:, :3] / gen_momentum[:, None]
    population_masks = (np.ones(len(gen_idx), dtype=bool), ~event["beam"])
    for channel in CHANNELS:
        objects = event["channels"][channel]
        p4, origins = objects["p4"], objects["gen_idx"]
        counts = accumulator["channels"][channel]
        momentum = np.linalg.norm(p4[:, :3], axis=1)
        directional = momentum > 0
        counts["candidates"] += len(p4)
        counts["zero_momentum_candidates"] += int(np.count_nonzero(~directional))
        counts["energy_invalid_candidates"] += int(np.count_nonzero(~objects["energy_valid"]))
        p4, origins = p4[directional], origins[directional]
        energy_valid = objects["energy_valid"][directional]
        reco_unit = p4[:, :3] / momentum[directional, None]
        cosine = np.clip(gen_unit @ reco_unit.T, -1, 1)
        cones = cosine[:, :, None] >= np.cos(np.deg2rad(CUTS_DEG))[None, None, :]
        own = gen_idx[:, None] == origins[None, :]
        other_isr = np.isin(origins, gen_idx)[None, :] & ~own
        non_isr = ((origins >= 0) & ~np.isin(origins, gen_idx))[None, :]
        unresolved = (origins < 0)[None, :]
        truth_cones = cones & own[:, :, None]
        multiplicity = cones.sum(axis=1)
        truth_multiplicity = truth_cones.sum(axis=1)
        energy_ratio = np.einsum("ijk,j->ik", cones, p4[:, 3]) / gen_p4[:, 3, None]
        truth_ratio = np.einsum("ijk,j->ik", truth_cones, p4[:, 3]) / gen_p4[:, 3, None]
        # Omit an energy-profile entry if its cone contains a representation
        # whose energy could double-count a saved parent/daughter structure.
        energy_rows = ~(cones & ~energy_valid[None, :, None]).any(axis=1)
        truth_energy_rows = ~(truth_cones & ~energy_valid[None, :, None]).any(axis=1)
        for population, selected in zip(POPULATIONS, population_masks, strict=True):
            values = counts["populations"][population]
            nominal = objects["stats"]["reco_count"][selected] > 0
            values["gen_count"] += int(selected.sum())
            values["truth_success"] += int(nominal.sum())
            values["geometric_success"] += cones[selected].any(axis=1).sum(axis=0)
            values["recovered_truth_success"] += truth_cones[selected].any(axis=1).sum(axis=0)
            values["multiple_candidates"] += (multiplicity[selected] > 1).sum(axis=0)
            for index, category in enumerate((own, other_isr, non_isr, unresolved)):
                values["pair_counts"][index] += (cones & category[:, :, None])[selected].sum(axis=(0, 1))
            for prefix, ratios, valid in (("energy", energy_ratio, energy_rows),
                                           ("truth_energy", truth_ratio, truth_energy_rows)):
                values[f"{prefix}_count"] += valid[selected].sum(axis=0)
                values[f"{prefix}_sum"] += np.where(valid[selected], ratios[selected], 0).sum(axis=0)
                values[f"{prefix}_square_sum"] += np.where(valid[selected], ratios[selected] ** 2, 0).sum(axis=0)
            values["multiplicity_sum"] += multiplicity[selected].sum(axis=0)
            values["multiplicity_square_sum"] += (multiplicity[selected] ** 2).sum(axis=0)
            values["truth_multiplicity_sum"] += truth_multiplicity[selected].sum(axis=0)
            values["truth_multiplicity_square_sum"] += (truth_multiplicity[selected] ** 2).sum(axis=0)


def mean_and_sem(counts, sums, squares):
    """Unweighted Gen-unit means, including zeros; SEM uses the sample variance."""
    counts = np.broadcast_to(counts, np.shape(sums))
    mean = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    variance_sum = squares - np.divide(sums ** 2, counts, out=np.zeros_like(sums), where=counts > 0)
    sem = np.sqrt(np.divide(np.maximum(variance_sum, 0), counts * (counts - 1),
                            out=np.full_like(sums, np.nan), where=counts > 1))
    return mean, sem


def summary(accumulator):
    """JSON-ready raw counts preserve the denominator of every scan quantity."""
    result = dict(events=accumulator["events"], cuts_degrees=CUTS_DEG.tolist(),
                  cone_rule="all candidates in each cone; a candidate may enter multiple Gen cones",
                  zero_momentum_rule="excluded only from geometric cones; retained in nominal truth success",
                  energy_rule="all Gen units including empty-cone zeros; omit cones containing energy-invalid objects",
                  channels={})
    for channel in CHANNELS:
        counts = accumulator["channels"][channel]
        populations = {}
        for population in POPULATIONS:
            values = counts["populations"][population]
            pairs = values["pair_counts"]
            assert np.all(values["recovered_truth_success"] <= values["truth_success"])
            entries = {field: value.tolist() if isinstance(value, np.ndarray) else value
                       for field, value in values.items() if field != "pair_counts"}
            entries["pair_counts"] = dict(zip(PAIR_CATEGORIES, pairs.tolist(), strict=True))
            assert np.array_equal(pairs.sum(axis=0), values["multiplicity_sum"])
            for prefix in ("energy", "truth_energy", "multiplicity", "truth_multiplicity"):
                denominator = values[f"{prefix}_count"] if prefix in ("energy", "truth_energy") else values["gen_count"]
                mean, sem = mean_and_sem(denominator, values[f"{prefix}_sum"], values[f"{prefix}_square_sum"])
                entries[f"{prefix}_mean"] = [float(value) if np.isfinite(value) else None for value in mean]
                entries[f"{prefix}_sem"] = [float(value) if np.isfinite(value) else None for value in sem]
            for field, numerator, denominator in (
                ("truth_recovery_given_truth_success", values["recovered_truth_success"], values["truth_success"]),
                ("truth_recovery_all_gen", values["recovered_truth_success"], values["gen_count"]),
                ("geometric_probability_all_gen", values["geometric_success"], values["gen_count"]),
            ):
                entries[field] = (numerator / denominator).tolist() if denominator else [None] * len(CUTS_DEG)
            populations[population] = entries
        result["channels"][channel] = {field: value for field, value in counts.items() if field != "populations"}
        result["channels"][channel]["populations"] = populations
    return result


def plot_angle_matching(output_root, sample, accumulator):
    """Four compact scan observables; energy also has a fully linear sibling."""
    import matplotlib.pyplot as plt

    from .plot_isr_study import FIGURE_SIZE, clopper_pearson, finish

    for number, channel in enumerate(CHANNELS, start=1):
        for population in POPULATIONS:
            values = accumulator["channels"][channel]["populations"][population]
            output = output_root / f"0{number}_{channel}" / population / sample
            output.mkdir(parents=True, exist_ok=True)
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for index, (numerator, denominator, label) in enumerate((
                (values["recovered_truth_success"], values["truth_success"], "Truth recovery"),
                (values["geometric_success"], values["gen_count"], "Any candidate"),
            )):
                denominators = np.full(len(CUTS_DEG), denominator)
                ratio = np.divide(numerator, denominators, out=np.full(len(CUTS_DEG), np.nan), where=denominators > 0)
                lower, upper = clopper_pearson(numerator, denominators)
                ax.errorbar(CUTS_DEG, ratio, yerr=[ratio - lower, upper - ratio], fmt='o-',
                            color=f'C{index}', label=label, capsize=2, markersize=4)
            ax.set_ylim(0, 1.2)
            ax.set_yticks(np.linspace(0, 1, 6))
            ax.set_xlim(0, CUTS_DEG[-1] * 1.05)
            finish(fig, ax, output / 'recovery_vs_angle.png', sample, 'Opening-angle cut [deg]', 'Probability')

            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            total = values["pair_counts"].sum(axis=0)
            fractions = np.divide(values["pair_counts"], total[None, :],
                                  out=np.zeros_like(values["pair_counts"], dtype=float), where=total[None, :] > 0)
            for index, label in enumerate(('Own ISR', 'Other ISR', 'Non-ISR', 'Unresolved')):
                ax.plot(CUTS_DEG, fractions[index], 'o-', label=label, color=f'C{index}', markersize=4)
            ax.set_ylim(0, 1.2)
            ax.set_yticks(np.linspace(0, 1, 6))
            ax.set_xlim(0, CUTS_DEG[-1] * 1.05)
            finish(fig, ax, output / 'pair_composition_vs_angle.png', sample, 'Opening-angle cut [deg]', 'Pair fraction')

            for linear in (False, True):
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                for index, (prefix, label) in enumerate((('energy', 'All candidates'), ('truth_energy', 'Own ISR'))):
                    mean, sem = mean_and_sem(values[f'{prefix}_count'], values[f'{prefix}_sum'], values[f'{prefix}_square_sum'])
                    ax.errorbar(CUTS_DEG, mean, yerr=sem, fmt='o-', color=f'C{index}', label=label, capsize=2, markersize=4)
                if not linear:
                    ax.set_yscale('log')
                else:
                    ax.set_ylim(bottom=0)
                suffix = '_linear' if linear else ''
                ax.set_xlim(0, CUTS_DEG[-1] * 1.05)
                finish(fig, ax, output / f'energy_recovery_vs_angle{suffix}.png', sample,
                       'Opening-angle cut [deg]', r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$')

            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for index, (prefix, label) in enumerate((('multiplicity', 'All candidates'), ('truth_multiplicity', 'Own ISR'))):
                mean, sem = mean_and_sem(values['gen_count'], values[f'{prefix}_sum'], values[f'{prefix}_square_sum'])
                ax.errorbar(CUTS_DEG, mean, yerr=sem, fmt='o-', label=label, color=f'C{index}', capsize=2, markersize=4)
            ax.set_ylim(bottom=0)
            ax.set_xlim(0, CUTS_DEG[-1] * 1.05)
            finish(fig, ax, output / 'multiplicity_vs_angle.png', sample, 'Opening-angle cut [deg]', r'$\langle N^{\rm reco}\rangle$ / gen')
