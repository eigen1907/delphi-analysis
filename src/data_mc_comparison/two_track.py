"""Common two-track checks and separate dimuon/dielectron control selections."""
import argparse
from collections import Counter
import csv
from pathlib import Path
import sqlite3

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot
from delphi_metadata.db import get_db_path

from common.paths import DATA_ROOT, PLOT_ROOT
from common.plot_utils import SAMPLE_LABELS, sample_files, set_style
from common.samples import SAMPLE_NAMES, sample_name
from sample_preparation.run_quality import event_keys, QUALITY_MIN, QUALITY_MAX
from .selections import CHANNELS, DETECTORS, STAGES, select_pairs

REGIONS = ("All", "BB", "FF", "BF", "Other")
COMMON_STAGES = ("input", "tracking_quality", "two_charged", "two_os", "back_to_back")
BINS = dict(mass=np.arange(0, 142, 2), p=np.arange(0, 71, 1),
            cos_theta=np.linspace(-1, 1, 51),
            hpc_e_over_p=np.linspace(0, 2, 81))
LABELS = dict(p=r"$p$ [GeV]", cos_theta=r"$\cos\theta$",
              hpc_e_over_p=r"$E_{\rm HPC}/p$")
CHANNEL_DIRS = dict(mumu="01_mumu", ee="02_ee")
SAMPLE_COLORS = dict(Data="black", Zmumu="C0", Zee="C1", ZKK="C2", Zpipi="C3", Ztautau="C4")


def pair_tags(tree, collection, events, parts):
    """Join native PID rows to Part indices; missing/ambiguous links stay -1."""
    if not len(events):
        return np.empty((0, 2), dtype=int)
    rows = tree[f"{collection}_partIdx"].array()[events]
    tags = tree[f"{collection}_tag"].array()[events]
    same_part = parts[:, :, None] == rows[:, None, :]
    count = ak.to_numpy(ak.sum(same_part, axis=-1))
    values = ak.to_numpy(ak.sum(ak.where(same_part, tags[:, None, :], 0), axis=-1))
    return np.where(count == 1, values, -1).astype(int)


def read_pair(path, tracking=None):
    with uproot.open(path) as root:
        tree = root["Events"]
        if tracking is None:
            tracking = np.ones(tree.num_entries, dtype=bool)
        charge = tree["Part_charge"].array()
        charged = charge != 0
        two = ak.to_numpy(ak.sum(charged, axis=1) == 2)
        q = np.asarray(ak.to_numpy(charge[charged][two]), dtype=int).reshape(-1, 2)
        opposite = q[:, 0] * q[:, 1] == -1
        events = np.flatnonzero(two)[opposite]
        parts = np.asarray(ak.to_numpy(ak.local_index(charge)[charged][two]), dtype=int).reshape(-1, 2)[opposite]
        components = [tree[f"Part_fourMomentum.fCoordinates.f{axis}"].array() for axis in "XYZ"]
        charged_xyz = [ak.values_astype(values[charged], np.float64) for values in components]
        charged_p = np.sqrt(sum(values**2 for values in charged_xyz))
        with np.errstate(divide="ignore", invalid="ignore"):
            charged_cos_theta = charged_xyz[2] / charged_p
        charged_energy = ak.values_astype(tree["Part_fourMomentum.fCoordinates.fT"].array()[charged], np.float64)
        # With exactly two charged tracks, this is simply that selected pair.
        leading = ak.argsort(charged_p, axis=1, ascending=False)[:, :2]
        mass2 = ak.to_numpy(ak.sum(charged_energy[leading], axis=1)**2
                            - sum(ak.sum(values[leading], axis=1)**2 for values in charged_xyz))
        n_charged = ak.to_numpy(ak.sum(charged, axis=1))
        rankable = ak.to_numpy(ak.all(np.isfinite(charged_p), axis=1))
        charged_mass = np.sqrt(np.where((n_charged >= 2) & rankable & (mass2 >= 0), mass2, np.nan))
        momentum = np.stack([
            np.asarray(ak.to_numpy(values[charged][two]), dtype=float).reshape(-1, 2)[opposite]
            for values in components
        ], axis=-1)
        stats = Counter(input_events=tree.num_entries, two_charged=int(two.sum()), two_os_raw=len(events))
        p2 = np.sum(momentum**2, axis=2)
        usable = np.isfinite(momentum).all(axis=(1, 2)) & (p2 > 0).all(axis=1)
        stats["invalid_momentum"] = int((~usable).sum())
        valid = usable
        events, parts, momentum = events[valid], parts[valid], momentum[valid]
        usable_events = np.zeros(tree.num_entries, dtype=bool)
        usable_events[events] = True
        common, common_counts = {}, {}
        for stage, selected in (("input", np.ones(tree.num_entries, dtype=bool)),
                                 ("tracking_quality", tracking), ("two_charged", tracking & two),
                                 ("two_os", tracking & usable_events)):
            common_counts[stage] = int(selected.sum())
            common[stage] = dict(p=ak.to_numpy(ak.flatten(charged_p[selected])),
                                 cos_theta=ak.to_numpy(ak.flatten(charged_cos_theta[selected])),
                                 mass=charged_mass[selected])
        stats["no_leading_pair"] = int((n_charged < 2).sum())
        stats["unrankable_leading_pair"] = int((~rankable).sum())
        stats["negative_leading_pair_mass2"] = int(((n_charged >= 2) & rankable & (mass2 < 0)).sum())
        tags = {key: pair_tags(tree, collection, events, parts)
                for key, collection in (("muid", "Muid"), ("elid", "Elid"))}
        hpc_energy = np.asarray(ak.to_numpy(tree["Part_hpcTotalShowerEnergy"].array()[charged][two]),
                                dtype=float).reshape(-1, 2)[opposite][valid]
    # Count all charged Parts first: an unusable third track cannot promote a two-track event.
    stats["valid_two_os"] = len(events)
    for key, values in tags.items():
        stats[f"unknown_{key}_legs"] = int((values < 0).sum())
    assert stats["two_os_raw"] == len(events) + stats["invalid_momentum"]
    return dict(events=events, momentum=momentum, tags=tags, hpc_energy=hpc_energy,
                native_mass=charged_mass[events], common=common, common_counts=common_counts), stats


def observables(pair, channel):
    xyz = pair["momentum"]
    momentum = np.linalg.norm(xyz, axis=2)
    energy = np.sqrt(momentum**2 + CHANNELS[channel]["mass"]**2)
    mass2 = energy.sum(axis=1)**2 - np.sum(xyz.sum(axis=1)**2, axis=1)
    theta = np.degrees(np.arctan2(np.hypot(xyz[:, :, 0], xyz[:, :, 1]), xyz[:, :, 2]))
    hpc = ((44 <= theta) & (theta <= 88)) | ((92 <= theta) & (theta <= 136))
    return dict(mass=np.sqrt(np.maximum(mass2, 0)), p=momentum,
                cos_theta=xyz[:, :, 2] / momentum,
                hpc_e_over_p=pair["hpc_energy"] / momentum,
                hpc=hpc)


def fill_histogram(result, key, entries):
    bins = BINS[key[-1]]
    finite = entries[np.isfinite(entries)]
    result["hist"].setdefault(key, np.zeros(len(bins) - 1, dtype=int))
    result["hist"][key] += np.histogram(finite, bins)[0]
    for name, count in (("entries", len(entries)), ("nonfinite", len(entries) - len(finite)),
                        ("underflow", (finite < bins[0]).sum()), ("overflow", (finite > bins[-1]).sum())):
        result["flows"][(*key, name)] += int(count)


def accumulate(files, quality=None):
    result = dict(hist={}, counts=Counter(), flows=Counter(), diagnostics=Counter(),
                  inputs=[], data=quality is not None)
    for index, path in enumerate(files, 1):
        event_flags, tracking = None, None
        if quality is not None:
            keys = event_keys(path.parent)
            event_flags = np.asarray([quality.get(key, (-1,) * len(DETECTORS)) for key in keys],
                                     dtype=int).reshape(-1, len(DETECTORS))
            tracking = ((event_flags[:, :4] >= QUALITY_MIN) & (event_flags[:, :4] <= QUALITY_MAX)).all(axis=1)
        pair, stats = read_pair(path, tracking)
        flags = None if event_flags is None else event_flags[pair["events"]]
        if flags is not None:
            stats["missing_quality_events"] = int((event_flags[:, 0] == -1).sum())
            stats["missing_quality_two_os"] = int((flags[:, 0] == -1).sum())
        result["diagnostics"].update(stats)
        result["inputs"].append((path, stats))
        for stage, values in pair["common"].items():
            result["counts"]["common", "All", stage] += pair["common_counts"][stage]
            for coordinate, entries in values.items():
                fill_histogram(result, ("common", "All", stage, coordinate), entries)
        for channel in CHANNELS:
            values = observables(pair, channel)
            category, masks = select_pairs(pair["momentum"], values["mass"],
                                            pair["tags"], channel, flags)
            if channel == "mumu":
                common_values = dict(p=values["p"], cos_theta=values["cos_theta"], mass=pair["native_mass"])
                selected = masks["back_to_back"]
                result["counts"]["common", "All", "back_to_back"] += int(selected.sum())
                for coordinate, entries in common_values.items():
                    fill_histogram(result, ("common", "All", "back_to_back", coordinate),
                                   entries[selected].reshape(-1))
            for earlier, later in zip(STAGES[channel][:-1], STAGES[channel][1:], strict=True):
                assert np.all(~masks[later] | masks[earlier])
            for region in REGIONS:
                in_region = np.ones(len(category), dtype=bool) if region == "All" else category == region
                for stage, selected in masks.items():
                    selected = selected & in_region
                    result["counts"][channel, region, stage] += int(selected.sum())
                    for coordinate in BINS:
                        if coordinate == "hpc_e_over_p":
                            if channel != "ee":
                                continue
                            entries = values[coordinate][selected[:, None] & values["hpc"]]
                        else:
                            entries = values[coordinate][selected].reshape(-1)
                        fill_histogram(result, (channel, region, stage, coordinate), entries)
        if index % 25 == 0 or index == len(files):
            print(f"  {index}/{len(files)} files; {result['diagnostics']['input_events']:,} events", flush=True)
    for key, histogram in result["hist"].items():
        assert histogram.sum() + sum(result["flows"][(*key, name)] for name in
                                     ("underflow", "overflow", "nonfinite")) == result["flows"][(*key, "entries")]
    for channel in CHANNELS:
        for stage in STAGES[channel]:
            assert result["counts"][channel, "All", stage] == sum(
                result["counts"][channel, region, stage] for region in REGIONS[1:])
    assert all(result["counts"]["mumu", "All", stage] == result["counts"]["ee", "All", stage]
               for stage in ("two_os", "back_to_back"))
    counts = [result["counts"]["common", "All", stage] for stage in COMMON_STAGES]
    assert all(earlier >= later for earlier, later in zip(counts[:-1], counts[1:], strict=True))
    return result


def decorate(ax, text, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.text(0.02, 0.97, text, transform=ax.transAxes, va="top", fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda index: labels[index] != "Data")
    ax.legend([handles[index] for index in order], [labels[index] for index in order],
              loc="upper right", frameon=False)


def save_plot(fig, ax, path):
    # Place the experiment label after layout, including plots with long tick labels.
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    mh.label.exp_label(exp="DELPHI", llabel=" Private Work", rlabel="1994", ax=ax)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_common_cutflow(results, output):
    fig, ax = plt.subplots(figsize=(16, 9))
    x = np.arange(len(COMMON_STAGES))
    maximum, minimum = 0, np.inf
    for name, result in results.items():
        counts = np.array([result["counts"]["common", "All", stage] for stage in COMMON_STAGES], dtype=float)
        maximum, minimum = max(maximum, counts.max()), min(minimum, counts.min())
        ax.errorbar(x, counts, yerr=np.sqrt(counts), fmt="o-", capsize=3,
                    label=SAMPLE_LABELS.get(name, name), color=SAMPLE_COLORS[name])
    ax.set_xticks(x, ("No cut", "TPC + VD\nquality", "2 charged\ntracks", "Opposite\nsign",
                     "Backtoback"))
    ax.set_xlim(-0.3, len(x) - 0.7)
    ax.set_yscale("log")
    ax.set_ylim(max(0.5, minimum * 0.5), maximum * 10)
    ax.grid(axis="y", alpha=0.2)
    decorate(ax, "", "Cut", "Events")
    save_plot(fig, ax, output / "00_cutflow.png")


def plot_channel_cutflow(results, channel, output):
    fig, ax = plt.subplots(figsize=(16, 9))
    stages = STAGES[channel][1:]
    x = np.arange(len(stages) + 1)
    maximum, minimum = 0, np.inf
    for name, result in results.items():
        counts = np.array([result["diagnostics"]["input_events"]]
                          + [result["counts"][channel, "All", stage] for stage in stages], dtype=float)
        maximum = max(maximum, counts.max())
        minimum = min(minimum, counts.min())
        color = SAMPLE_COLORS[name]
        ax.errorbar(x, counts, yerr=np.sqrt(counts), fmt="o-", capsize=3,
                    label=SAMPLE_LABELS.get(name, name), color=color)
    lepton = "muon" if channel == "mumu" else "electron"
    labels = dict(back_to_back="Common", detector_quality="Detector\nquality", fiducial="Acceptance",
                  standard_id="Both " + lepton + " IDs", mass_window=cut_text(channel, "mass_window")[0])
    ax.set_xticks(x, ("No cut", *(labels[stage] for stage in stages)))
    ax.set_xlim(-0.3, len(x) - 0.7)
    ax.set_yscale("log")
    ax.set_ylim(max(0.5, minimum * 0.5), maximum * 10)
    ax.grid(axis="y", alpha=0.2)
    decorate(ax, CHANNELS[channel]["label"], "Cut", "Events")
    save_plot(fig, ax, output / "00_cutflow.png")


def draw_distribution(ax, results, channel, stage, coordinate):
    bins = BINS[coordinate]
    centers = (bins[:-1] + bins[1:]) / 2
    maximum = 0
    for name, result in results.items():
        key = channel, "All", stage, coordinate
        counts = result["hist"][key]
        total = result["flows"][(*key, "entries")] - result["flows"][(*key, "nonfinite")]
        color = SAMPLE_COLORS[name]
        label = SAMPLE_LABELS.get(name, name)
        if total == 0:
            ax.plot([], [], color=color, label=f"{label} (0 entries)")
            continue
        density = counts / (total * np.diff(bins))
        error = np.sqrt(counts) / (total * np.diff(bins))
        maximum = max(maximum, (density + error).max())
        if name == "Data":
            ax.errorbar(centers, density, yerr=error, fmt="o", color=color, markersize=5, label=label)
        else:
            ax.stairs(density, bins, color=color, label=label)
            shown = counts > 0
            ax.errorbar(centers[shown], density[shown], yerr=error[shown], fmt="none", color=color)
    return maximum


def cut_text(channel, stage):
    """Only the requirement added at this stage."""
    common = {
        "input": [],
        "two_charged": [r"$N_{\mathrm{charged\ tracks}}=2$"],
        "two_os": [r"$q_1q_2=-1$", r"Finite $|\mathbf{p}_i|>0$"],
        "tracking_quality": ["TPC + VD quality: 6–7"],
        "back_to_back": [r"$\pi-\Delta\phi<0.35$ rad"],
    }
    if channel == "common":
        return common[stage].copy()
    if stage == "back_to_back":
        return ["Common"]
    if stage == "mass_window":
        return [r"$70<" + CHANNELS[channel]["mass_label"][1:-1] + r"<110$ GeV"]
    if stage == "fiducial":
        return ([r"$\theta_i\in[20,42]\cup[53,88.5]$", r"$\quad\cup[91.5,127]\cup[138,160]$ deg"]
                if channel == "mumu" else
                [r"$\theta_i\in[20,35]\cup[44,88]$", r"$\quad\cup[92,136]\cup[145,160]$ deg"])
    if stage == "detector_quality":
        detectors = ("MUB b,d (barrel); MUF A/C (endcap)" if channel == "mumu" else
                     "HPC 0,1 (barrel); EMF A/C (forward)")
        return [detectors, "Quality: 6–7"]
    if stage == "standard_id":
        return (["Both: Muid_tag >= 0 and (tag & 4) != 0"] if channel == "mumu" else
                ["Both: Elid_tag >= 4"])


def plot_stage(results, channel, stage, coordinate, path, show_cut=True):
    """One axis per cumulative stage; annotate only its new cut."""
    fig, ax = plt.subplots(figsize=(12, 11))
    bins = BINS[coordinate]
    maximum = draw_distribution(ax, results, channel, stage, coordinate)
    lines = cut_text(channel, stage)
    context = "" if channel == "common" else CHANNELS[channel]["label"]
    if coordinate == "hpc_e_over_p":
        lines.append("HPC tracks only")
    if coordinate == "mass":
        xlabel = r"$m_{12}$ [GeV]" if channel == "common" else CHANNELS[channel]["mass_label"] + " [GeV]"
    else:
        xlabel = LABELS[coordinate]
    ylabel = "Normalized / GeV" if coordinate in ("mass", "p") else "Normalized"
    decorate(ax, context, xlabel, ylabel)
    if show_cut:
        ax.text(0.02, 0.90 if context else 0.97, "\n".join(lines), transform=ax.transAxes, va="top")
    ax.set_xlim(bins[0], bins[-1])
    headroom = 1.6 + 0.14 * len(lines)
    ax.set_ylim(0, max(maximum, 0.01) * headroom)
    if coordinate == "hpc_e_over_p":
        ax.set_yscale("log")
        ax.set_ylim(1e-4, max(maximum, 0.01) * 100)
    save_plot(fig, ax, path)


def write_counts(results, output):
    with (output / "inputs.csv").open("w") as stream:
        columns = ("input_events", "two_charged", "two_os_raw", "invalid_momentum",
                   "valid_two_os", "missing_quality_events", "missing_quality_two_os", "unknown_muid_legs", "unknown_elid_legs",
                   "no_leading_pair", "unrankable_leading_pair", "negative_leading_pair_mass2")
        writer = csv.writer(stream)
        writer.writerow(("sample", "file", *columns))
        for name, result in results.items():
            for path, stats in result["inputs"]:
                writer.writerow((name, path.resolve(), *(stats[column] for column in columns)))
    with (output / "cutflow.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(("sample", "channel", "region", "stage", "events", "mass_overflow_gt_140", "applied"))
        for name, result in results.items():
            for stage in COMMON_STAGES:
                applied = result["data"] or stage != "tracking_quality"
                writer.writerow((name, "common", "All", stage, result["counts"]["common", "All", stage], "", applied))
            for channel in CHANNELS:
                stages = STAGES[channel]
                print(name, channel, {stage: result["counts"][channel, "All", stage] for stage in stages}, flush=True)
                for region in REGIONS:
                    for stage in stages:
                        applied = result["data"] or stage not in ("tracking_quality", "detector_quality")
                        writer.writerow((name, channel, region, stage, result["counts"][channel, region, stage],
                            result["flows"][channel, region, stage, "mass", "overflow"], applied))
    with (output / "histogram_flows.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(("sample", "channel", "region", "stage", "observable", "entries", "underflow", "overflow", "nonfinite"))
        for name, result in results.items():
            for key in result["hist"]:
                writer.writerow((name, *key, *(result["flows"][(*key, field)] for field in
                    ("entries", "underflow", "overflow", "nonfinite"))))


def plot_results(results, output):
    """Flat numbered files: each distribution is after the named cumulative cut."""
    set_style()
    common = output / "00_common"
    common.mkdir(exist_ok=True)
    plot_common_cutflow(results, common)
    common_stages = (("01_all_charged_parts", "input"), ("02_tracking_quality", "tracking_quality"),
                     ("03_two_charged_tracks", "two_charged"), ("04_opposite_sign", "two_os"),
                     ("05_backtoback", "back_to_back"), ("06_selected", "back_to_back"))
    for prefix, stage in common_stages:
        for coordinate in ("p", "cos_theta", "mass"):
            name = "inv_mass" if coordinate == "mass" else coordinate
            plot_stage(results, "common", stage, coordinate, common / f"{prefix}_{name}.png",
                       show_cut=prefix != "06_selected")
    for channel, directory in CHANNEL_DIRS.items():
        sample = "Zmumu" if channel == "mumu" else "Zee"
        selected = {name: result for name, result in results.items() if name in ("Data", sample)}
        folder = output / directory
        folder.mkdir(exist_ok=True)
        plot_channel_cutflow(selected, channel, folder)
        stages = ((("01_common", "back_to_back"), ("02_angular_acceptance", "fiducial"),
                   ("03_muon_id", "standard_id"), ("04_detector_quality", "detector_quality"))
                  if channel == "mumu" else
                  (("01_common", "back_to_back"), ("02_angular_acceptance", "fiducial"),
                   ("03_electron_id", "standard_id"), ("04_detector_quality", "detector_quality")))
        stages += (("05_mass_window", "mass_window"), ("06_selected", "mass_window"))
        for prefix, stage in stages:
            coordinates = ["p", "cos_theta", "mass"]
            if channel == "ee" and stage in ("fiducial", "standard_id"):
                coordinates.append("hpc_e_over_p")
            for coordinate in coordinates:
                name = "inv_mass" if coordinate == "mass" else coordinate
                suffix = "_log" if coordinate == "hpc_e_over_p" else ""
                plot_stage(selected, channel, stage, coordinate, folder / f"{prefix}_{name}{suffix}.png",
                           show_cut=prefix != "06_selected")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", choices=SAMPLE_NAMES, default=["Zmumu", "Zee"],
                        help="MC templates (default: Zmumu Zee); Data selections always run both channels")
    parser.add_argument("--data", type=Path, default=DATA_ROOT / "20260828_florian/20261026_short94_c2")
    parser.add_argument("--mc", type=Path, default=DATA_ROOT / "20260828_florian")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or PLOT_ROOT / args.mc.name / "03_data_mc_comparison/two_track"
    output.mkdir(parents=True, exist_ok=True)
    data_files = []
    for path in sorted(args.data.glob("*/nanoaod.root")):
        raw, log = path.parent / "nanoaod_raw_sdst.root", path.parent / "conversion.log"
        if path.stat().st_size and raw.is_file() and raw.stat().st_size and log.is_file() \
                and "Native/raw physics keys: exact 1:1" in log.read_text():
            data_files.append(path)
    assert data_files, "No validated collision conversions found"
    print(f"Data snapshot: {len(data_files)} completed files", flush=True)
    with sqlite3.connect(f"file:{get_db_path()}?mode=ro", uri=True) as db:
        quality = {(run, sequence): flags for run, sequence, *flags in db.execute(
            f"SELECT run, cartridge, {', '.join(DETECTORS)} FROM run_quality WHERE year = 1994")}
    results = {"Data": accumulate(data_files, quality)}
    mc_files = {sample_name(dataset): paths for dataset, paths in sample_files(args.mc, "nanoaod.root")
                if dataset.endswith("_photosFSR") and sample_name(dataset) in SAMPLE_NAMES}
    for name in args.samples:
        print(name, "PHOTOS-FSR", flush=True)
        results[name] = accumulate(mc_files[name])
    write_counts(results, output)
    plot_results(results, output)
    print(f"Plots and cutflow: {output}")


if __name__ == "__main__":
    main()
