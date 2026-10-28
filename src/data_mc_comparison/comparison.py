"""Stored Raw and Native inputs in the final native two-lepton selections."""
from collections import Counter
import csv
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np

from common.plot_utils import SAMPLE_LABELS, clean_name, set_style
from .high_level_variables import HIGH_LEVEL_VARIABLES, UNDEFINED_FEATURES
from .variables import VARIABLES

CHANNELS = (("mumu", "01_mumu", "Zmumu", "C0"), ("ee", "02_ee", "Zee", "C1"))
LINK_COUNTS = dict(EmLayer_emShowerIdx="nEmShower", HadHit_hadShowerIdx="nHadShower",
                   VdAssocHit_tracRawIdx="nTracRaw")
FIELDS = ("entries", "finite", "zero", "negative", "nonfinite", "underflow", "overflow",
          "link_unlinked", "link_out_of_range")
LEVELS = (("low_level", "comparison", VARIABLES),
          ("high_level", "high_level_comparison", HIGH_LEVEL_VARIABLES))


def fill(hist, flows, key, selected, bins, bounds=None):
    """Keep finite source values and account for every stored or missing entry."""
    entries = np.asarray(ak.to_numpy(ak.flatten(
        ak.fill_none(selected, np.nan, axis=None), axis=None)), dtype=float)
    hist.setdefault(key, np.zeros(len(bins) - 1, dtype=np.int64))
    finite = entries[np.isfinite(entries)]
    hist[key] += np.histogram(finite, bins)[0]
    numbers = dict(entries=len(entries), finite=len(finite), zero=(finite == 0).sum(),
                   negative=(finite < 0).sum(), nonfinite=len(entries) - len(finite),
                   underflow=(finite < bins[0]).sum(), overflow=(finite > bins[-1]).sum())
    if bounds is not None:
        bounds = np.asarray(ak.to_numpy(ak.flatten(bounds, axis=None)))
        numbers["link_unlinked"] = (entries == -1).sum()
        numbers["link_out_of_range"] = ((entries < -1) | (entries >= bounds)).sum()
    for field, value in numbers.items():
        flows[(*key, field)] += int(value)
    assert hist[key].sum() + flows[(*key, "underflow")] + flows[(*key, "overflow")] \
        == flows[(*key, "finite")]


def accumulate(raw_tree, rows_by_channel, result):
    """Fill all stored objects, or one scalar per event; retain finite sentinels."""
    hist = result.setdefault("hist", {})
    flows = result.setdefault("flows", Counter())
    counts = result.setdefault("counts", Counter())
    rows_by_channel = {channel: np.asarray(rows, dtype=np.int64)
                       for channel, rows in rows_by_channel.items()}
    for channel, rows in rows_by_channel.items():
        if rows.ndim != 1 or len(np.unique(rows)) != len(rows) or np.any(
                (rows < 0) | (rows >= raw_tree.num_entries)):
            raise ValueError(f"Invalid or repeated selected raw rows for {channel}")
        counts[channel] += len(rows)
    branches = {branch for specs in VARIABLES.values() for branch, _, _, _ in specs}
    extra_names = {"EmShower_detector", *LINK_COUNTS.values()}
    missing = sorted(name for name in branches | extra_names if name not in raw_tree)
    if missing:
        raise KeyError(f"Required Raw comparison branches missing: {', '.join(missing)}")
    extra = {name: raw_tree[name].array(library="ak") for name in extra_names}
    for group, specs in VARIABLES.items():
        for branch, bins, _, _ in specs:
            values = extra[branch] if branch in extra else raw_tree[branch].array(library="ak")
            for channel, rows in rows_by_channel.items():
                selected = values[rows]
                if branch == "EmShower_nLayers":
                    selected = selected[extra["EmShower_detector"][rows] == 9]
                key = channel, group, branch
                bounds = None
                if branch in LINK_COUNTS:
                    bounds = ak.broadcast_arrays(selected, extra[LINK_COUNTS[branch]][rows])[1]
                fill(hist, flows, key, selected, bins, bounds)


def accumulate_native(native_tree, selections, result):
    """Fill two selected Part legs, all event vertices, and event scalars.

    Part, Dedx, Haid, Rich and Haidc use Part order. Track and lepton-ID
    banks require a unique explicit Part link; missing or ambiguous links
    contribute nonfinite entries rather than disappearing from the comparison.
    """
    hist = result.setdefault("hist", {})
    flows = result.setdefault("flows", Counter())
    counts = result.setdefault("counts", Counter())
    branches = {branch for specs in HIGH_LEVEL_VARIABLES.values() for branch, _, _, _ in specs}
    link_names = {f"{branch.split('_', 1)[0]}_partIdx" for branch in branches
                  if branch.startswith(("Trac_", "Muid_", "Elid_"))}
    missing = sorted(name for name in branches | link_names | {"nPart"} if name not in native_tree)
    if missing:
        raise KeyError(f"Required Native comparison branches missing: {', '.join(missing)}")
    n_part = native_tree["nPart"].array(library="np")
    selections = {channel: dict(events=np.asarray(selection["events"], dtype=np.int64),
                                parts=np.asarray(selection["parts"], dtype=np.int64))
                  for channel, selection in selections.items()}
    for channel, selection in selections.items():
        rows, parts = selection["events"], selection["parts"]
        if rows.ndim != 1 or len(np.unique(rows)) != len(rows) or np.any(
                (rows < 0) | (rows >= native_tree.num_entries)):
            raise ValueError(f"Invalid or repeated selected Native rows for {channel}")
        if parts.shape != (len(rows), 2) or np.any(parts < 0) or np.any(
                parts >= n_part[rows, None]) or np.any(parts[:, 0] == parts[:, 1]):
            raise ValueError(f"Invalid selected Native Part pair for {channel}")
        counts[channel] += len(rows)
    links = {name: native_tree[name].array(library="ak") for name in link_names}
    joins = {}
    for channel, selection in selections.items():
        rows, parts = selection["events"], selection["parts"]
        if not len(rows):
            continue
        for name, all_links in links.items():
            selected_links = all_links[rows]
            matches = selected_links[:, :, None] == ak.Array(parts[:, None, :])
            matches_per_leg = ak.sum(matches, axis=1)
            indices = ak.sum(ak.where(matches, ak.local_index(selected_links)[:, :, None], 0), axis=1)
            joins[channel, name] = ak.mask(indices, matches_per_leg == 1)
    for group, specs in HIGH_LEVEL_VARIABLES.items():
        for branch, bins, _, _ in specs:
            values = n_part if branch == "nPart" else native_tree[branch].array(library="ak")
            for channel, selection in selections.items():
                rows, parts = selection["events"], selection["parts"]
                selected = ak.Array(values[rows])
                if not len(rows):
                    selected = ak.Array([])
                elif branch.startswith(("Part_", "Dedx_", "Haid_", "Rich_", "Haidc_")):
                    assert np.array_equal(ak.to_numpy(ak.num(selected, axis=1)), n_part[rows]), \
                        f"Native {branch} does not follow Part order in {channel} selected events"
                    selected = selected[ak.from_regular(parts, axis=1)]
                elif branch.startswith(("Trac_", "Muid_", "Elid_")):
                    link_name = f"{branch.split('_', 1)[0]}_partIdx"
                    assert np.array_equal(ak.to_numpy(ak.num(selected, axis=1)),
                                          ak.to_numpy(ak.num(links[link_name][rows], axis=1))), \
                        f"Native {branch} and {link_name} lengths differ in {channel} selected events"
                    selected = selected[joins[channel, link_name]]
                elif not branch.startswith(("Vtx_", "Event_", "n")):
                    raise ValueError(f"No Native comparison scope defined for {branch}")
                fill(hist, flows, (channel, group, branch), selected, bins)


def scope(level, branch):
    if branch.startswith(("Event_", "n")):
        return "one scalar per selected event"
    if level == "low_level":
        return "all stored objects in selected events"
    if branch.startswith("Vtx_"):
        return "all Native vertices in selected events"
    return "two selected Native leptons"


def plot(results, output):
    """Write paired Data/matching-MC densities, ratios, and explicit storage flows."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    set_style()
    with (output / "unavailable_features.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("feature", "definition", "status"))
        writer.writerows((name, "", "Not computed: definition pending") for name in UNDEFINED_FEATURES)
    with (output / "histogram_flows.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("sample", "channel", "level", "group", "branch", "selected_events",
                         *FIELDS, "scope", "note"))
        for name, sample in results.items():
            for level, result_key, variables in LEVELS:
                result = sample[result_key]
                for group, specs in variables.items():
                    for branch, _, _, note in specs:
                        if level == "high_level" and branch.startswith(("Trac_", "Muid_", "Elid_")):
                            note = "; ".join(filter(None, (note,
                                "Missing or ambiguous Part links counted as nonfinite; stored zero tags retained")))
                        for channel, _, _, _ in CHANNELS:
                            key = channel, group, branch
                            writer.writerow((name, channel, level, group, branch, result["counts"][channel],
                                *(result["flows"][(*key, field)] for field in FIELDS), scope(level, branch), note))
    for channel, directory, mc_name, color in CHANNELS:
        if "Data" not in results or mc_name not in results:
            continue
        for level, result_key, variables in LEVELS:
            for group, specs in variables.items():
                folder = output / directory / level / group
                folder.mkdir(parents=True, exist_ok=True)
                for index, (branch, bins, xlabel, note) in enumerate(specs, 1):
                    key = channel, group, branch
                    bins = np.asarray(bins)
                    centers, widths = (bins[1:] + bins[:-1]) / 2, np.diff(bins)
                    fig, (ax, ratio_ax) = plt.subplots(2, 1, figsize=(12, 11), sharex=True,
                        gridspec_kw=dict(height_ratios=(3, 1), hspace=0.05))
                    densities, errors, raw_counts = {}, {}, {}
                    maximum = 0
                    for name, sample_color in (("Data", "black"), (mc_name, color)):
                        result = results[name][result_key]
                        counts = result["hist"].get(key, np.zeros(len(bins) - 1, dtype=int))
                        total = result["flows"][(*key, "finite")]
                        density = counts / (total * widths) if total else np.zeros_like(widths)
                        error = np.sqrt(counts) / (total * widths) if total else np.zeros_like(widths)
                        densities[name], errors[name], raw_counts[name] = density, error, counts
                        maximum = max(maximum, np.max(density + error))
                        label = SAMPLE_LABELS.get(name, name) + (" (0 finite entries)" if not total else "")
                        if name == "Data":
                            ax.errorbar(centers, density, yerr=error, fmt="o", color=sample_color,
                                        markersize=4, label=label)
                        else:
                            ax.stairs(density + error, bins, baseline=np.maximum(0, density - error),
                                      fill=True, color=sample_color, alpha=0.2, linewidth=0)
                            ax.stairs(density, bins, color=sample_color, label=label)
                    valid = densities[mc_name] > 0
                    relative_error = np.full(len(widths), np.nan)
                    relative_error[valid] = 1 / np.sqrt(raw_counts[mc_name][valid])
                    ratio_ax.stairs(1 + relative_error, bins, baseline=1 - relative_error,
                                    fill=True, color=color, alpha=0.25, linewidth=0)
                    valid &= results["Data"][result_key]["flows"][(*key, "finite")] > 0
                    ratio, ratio_error = (densities["Data"][valid] / densities[mc_name][valid],
                                          errors["Data"][valid] / densities[mc_name][valid])
                    ratio_ax.errorbar(centers[valid], ratio, yerr=ratio_error, fmt="o", color="black",
                                      markersize=4)
                    ratio_ax.axhline(1, color=color, linewidth=1)
                    ratio_ax.set_ylim(0, max(2, np.max(ratio + ratio_error) * 1.15) if len(ratio) else 2)
                    ratio_ax.set_ylabel("Data / MC")
                    ratio_ax.set_xlabel(xlabel)
                    ratio_ax.grid(axis="y", alpha=0.2)
                    ax.set_ylabel("Normalized entries / bin width")
                    ax.set_ylim(0, max(maximum, 0.01) * 1.5)
                    ax.set_xlim(bins[0], bins[-1])
                    ax.text(0.02, 0.97, SAMPLE_LABELS[mc_name], transform=ax.transAxes,
                            va="top", fontweight="bold")
                    ax.text(0.02, 0.88, branch, transform=ax.transAxes, va="top")
                    handles, labels = ax.get_legend_handles_labels()
                    order = sorted(range(len(labels)), key=lambda i: not labels[i].startswith("Data"))
                    ax.legend([handles[i] for i in order], [labels[i] for i in order],
                              loc="upper right", frameon=False)
                    ax.tick_params(labelbottom=False)
                    fig.subplots_adjust(left=0.14, right=0.96, bottom=0.13, top=0.90, hspace=0.06)
                    mh.label.exp_label(exp="DELPHI", llabel=" Private Work", rlabel="1994", ax=ax)
                    fig.savefig(folder / f"{index:02d}_{clean_name(branch)}.png", dpi=150)
                    plt.close(fig)
