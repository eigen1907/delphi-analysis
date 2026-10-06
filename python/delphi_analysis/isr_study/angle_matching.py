"""Opening-angle truth validation and independent geometric matching."""

import numpy as np

CUTS_DEG = np.array([0.25, 0.5, 1, 2, 3, 5, 10, 20])
CHANNELS = ("reco_gamma", "reco_gamma_plus_conv")
RESULT_CHANNELS = (*CHANNELS, "reco_all")
POPULATIONS = ("gen_gamma_all", "gen_isr_all", "gen_isr_non_beam")
PAIR_CATEGORIES = ("own_gen", "other_ISR", "non_ISR", "unresolved")
GEN_LABELS = ("All gamma", "All ISR", "Non-beam ISR")
RECO_LABELS = (r"$\gamma$", r"$\gamma$ + conv")


def new_accumulator():
    """Keep scan moments and compact per-Gen angular observables."""
    result = dict(events=0, channels={}, angular={})
    for channel in RESULT_CHANNELS:
        result["angular"][channel] = {
            name: dict(pending=[], blocks=[])
            for name in ("reco_count", "energy_count", "energy_ratio")
        }
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
    """Scan the same global Reco candidates for all three Gen populations."""
    accumulator["events"] += 1
    gen_idx, gen_p4 = event["all_gamma_idx"], event["all_gamma_p4"]
    gen_momentum = np.linalg.norm(gen_p4[:, :3], axis=1)
    assert np.all(gen_momentum > 0) and np.all(gen_p4[:, 3] > 0)
    gen_unit = gen_p4[:, :3] / gen_momentum[:, None]
    isr_idx = gen_idx[event["selections"]["gen_isr_all"]]
    for channel in RESULT_CHANNELS:
        objects = event["reco"][channel]
        p4, origins = objects["p4"], objects["gen_idx"]
        momentum = np.linalg.norm(p4[:, :3], axis=1)
        directional = momentum > 0
        if channel in CHANNELS:
            counts = accumulator["channels"][channel]
            counts["candidates"] += len(p4)
            counts["zero_momentum_candidates"] += int(np.count_nonzero(~directional))
            counts["energy_invalid_candidates"] += int(np.count_nonzero(~objects["energy_valid"]))
        p4, origins = p4[directional], origins[directional]
        energy_valid = objects["energy_valid"][directional]
        geometric_energy_valid = objects["geometric_energy_valid"][directional]
        reco_unit = p4[:, :3] / momentum[directional, None]
        cosine = np.clip(gen_unit @ reco_unit.T, -1, 1)
        cones = cosine[:, :, None] >= np.cos(np.deg2rad(CUTS_DEG))[None, None, :]
        energy_ratio = np.einsum("ijk,j->ik", cones, p4[:, 3]) / gen_p4[:, 3, None]
        energy_rows = ~(cones & ~energy_valid[None, :, None]).any(axis=1)
        geometric_energy_rows = ~(cones & ~geometric_energy_valid[None, :, None]).any(axis=1)
        multiplicity = cones.sum(axis=1)
        raw_multiplicity = multiplicity
        if channel == "reco_all":
            # Stored Parts count as footprints; energy uses canonical objects.
            raw_p4 = objects["raw_p4"]
            raw_momentum = np.linalg.norm(raw_p4[:, :3], axis=1)
            raw_unit = raw_p4[raw_momentum > 0, :3] / raw_momentum[raw_momentum > 0, None]
            raw_cosine = np.clip(gen_unit @ raw_unit.T, -1, 1)
            raw_multiplicity = (raw_cosine[:, :, None] >= np.cos(np.deg2rad(CUTS_DEG))).sum(axis=1)
        stored = accumulator["angular"][channel]
        assert np.all(raw_multiplicity <= np.iinfo(np.uint16).max)
        assert np.all(multiplicity <= np.iinfo(np.uint16).max)
        for name, values in (("reco_count", raw_multiplicity.astype(np.uint16)),
                             ("energy_count", multiplicity.astype(np.uint16)),
                             ("energy_ratio", np.where(geometric_energy_rows, energy_ratio, np.nan))):
            stored[name]["pending"].append(values)
        # Consolidate small event arrays so the full scan stays compact.
        if len(stored["reco_count"]["pending"]) >= 1024:
            for field in stored.values():
                field["blocks"].append(np.concatenate(field["pending"]))
                field["pending"].clear()
        if channel == "reco_all":
            continue
        own = gen_idx[:, None] == origins[None, :]
        isr_origin = np.isin(origins, isr_idx)[None, :]
        other_isr = isr_origin & ~own
        non_isr = ((origins >= 0)[None, :] & ~isr_origin) & ~own
        unresolved = (origins < 0)[None, :]
        truth_cones = cones & own[:, :, None]
        truth_multiplicity = truth_cones.sum(axis=1)
        truth_ratio = np.einsum("ijk,j->ik", truth_cones, p4[:, 3]) / gen_p4[:, 3, None]
        # A saved parent/daughter ambiguity invalidates its energy, not its
        # presence. Empty cones remain zero-energy entries.
        truth_energy_rows = ~(truth_cones & ~energy_valid[None, :, None]).any(axis=1)
        truth_energy_rows &= objects["stats"]["energy_valid"][:, None]
        for population in POPULATIONS:
            selected = event["selections"][population]
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


def angular_stats(accumulator, cut_degrees):
    """Return per-Gen cone matching at the selected scan cut, with no truth gate.

    Arrays stay aligned with all stable Gen photons. A Reco candidate may enter
    several Gen cones; no unique assignment or angular truth repair is implied.
    Extracting the selected column releases the other scan-cut observations.
    """
    index = int(np.flatnonzero(CUTS_DEG == cut_degrees).item())
    result = {}
    for channel, fields in accumulator.pop("angular").items():
        stats = {}
        for name, field in fields.items():
            if field["pending"]:
                field["blocks"].append(np.concatenate(field["pending"]))
                field["pending"].clear()
            stats[name] = np.concatenate([block[:, index] for block in field["blocks"]])
            field["blocks"].clear()
        stats["raw_count"] = stats["reco_count"].copy()
        stats["energy_valid"] = np.isfinite(stats["energy_ratio"])
        result[channel] = stats
    return result


def mean_and_sem(counts, sums, squares):
    """Gen-unit means including zeros; SEM uses the sample variance."""
    counts = np.broadcast_to(counts, np.shape(sums))
    mean = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    variance_sum = squares - np.divide(sums ** 2, counts, out=np.zeros_like(sums), where=counts > 0)
    sem = np.sqrt(np.divide(np.maximum(variance_sum, 0), counts * (counts - 1),
                            out=np.full_like(sums, np.nan), where=counts > 1))
    return mean, sem


def choose_cut(accumulators):
    """Smallest cut retaining 99% of pooled non-beam truth successes in both channels."""
    pooled = {channel: dict(truth_success=0, recovered=np.zeros(len(CUTS_DEG), dtype=np.int64))
              for channel in CHANNELS}
    for accumulator in accumulators:
        for channel in CHANNELS:
            values = accumulator["channels"][channel]["populations"]["gen_isr_non_beam"]
            pooled[channel]["truth_success"] += values["truth_success"]
            pooled[channel]["recovered"] += values["recovered_truth_success"]
    eligible = np.ones(len(CUTS_DEG), dtype=bool)
    for values in pooled.values():
        eligible &= (values["truth_success"] > 0) & (values["recovered"] >= 0.99 * values["truth_success"])
    index = int(np.flatnonzero(eligible)[0]) if eligible.any() else len(CUTS_DEG) - 1
    channels = {}
    for channel, values in pooled.items():
        denominator = values["truth_success"]
        channels[channel] = dict(truth_success=denominator,
                                 recovered_truth_success=values["recovered"].tolist(),
                                 retained_fraction=(values["recovered"] / denominator).tolist() if denominator else [None] * len(CUTS_DEG))
    return dict(angle_degrees=float(CUTS_DEG[index]), target_truth_retention=0.99,
                target_met=bool(eligible.any()),
                selection_rule="smallest scan cut retaining >=99% of pooled non-beam ISR truth successes in both Photon channels",
                channels=channels)


def summary(accumulator):
    """Raw counts keep every denominator and cone-contamination diagnostic."""
    result = dict(events=accumulator["events"], cuts_degrees=CUTS_DEG.tolist(),
                  efficiency_rule="own truth-associated Reco inside cone / all selected Gen photons",
                  cone_rule="all candidates in each cone; a candidate may enter multiple Gen cones",
                  zero_momentum_rule="excluded only from geometric cones; retained in nominal truth success",
                  energy_rule="own truth-associated Reco energy / Gen energy, including empty-cone zeros; omit ambiguous nominal Gen energy and cones containing energy-invalid objects",
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


def plot_angle_matching(output_root, sample, accumulator, cut_degrees):
    """Six cases share each plot; color is Gen selection and line style is Reco definition."""
    import matplotlib.pyplot as plt

    from .plot_isr_study import FIGURE_SIZE, clopper_pearson, finish

    output_root.mkdir(parents=True, exist_ok=True)
    shown = CUTS_DEG <= 10
    for observable in ("efficiency", "energy_recovery"):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for gen_index, population in enumerate(POPULATIONS):
            for reco_index, channel in enumerate(CHANNELS):
                values = accumulator["channels"][channel]["populations"][population]
                label = f"{GEN_LABELS[gen_index]} / {RECO_LABELS[reco_index]}"
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
                            linestyle=("-", "--")[reco_index], marker=("o", "s")[reco_index],
                            label=label, capsize=2, markersize=4)
        ax.axvline(cut_degrees, color="black", linestyle="--", linewidth=2)
        ax.set_xlim(0, 10)
        ax.set_ylim(bottom=0)
        ylabel = "Efficiency" if observable == "efficiency" else r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$'
        # Reserve an interior legend band; the data determine the useful scale.
        lower, upper = ax.get_ylim()
        ax.set_ylim(lower, upper * 1.65)
        ax.legend(title=sample, ncol=2, loc="upper center", frameon=True,
                  framealpha=1, edgecolor="none")
        finish(fig, ax, output_root / f"{observable}_vs_opening_angle_{sample}.png", sample,
               "Opening-angle cut [deg]", ylabel)
