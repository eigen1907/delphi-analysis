"""Inspect 1994 Run Quality using the email's quality codes 6 or 7."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot
from delphi_metadata.db import get_db_path

from common.paths import DATA_ROOT, DATASET_ROOT, PLOT_ROOT
from common.plot_utils import FIGURE_SIZE, set_style


# Detector groups follow the official dmd runquality legend.
TRACKING = ("tpc_0", "tpc_1", "mvx_a", "mvx_c")
BARREL = ("mub_b", "mub_d")
FORWARD = ("muf_a", "muf_c")
DETECTORS = TRACKING + BARREL + FORWARD
QUALITY_MIN, QUALITY_MAX = 6, 7


def plot_quality():
    """Compare full-year metadata requirements, without selecting event angles."""
    with sqlite3.connect(f"file:{get_db_path()}?mode=ro", uri=True) as db:
        quality = db.execute(
            f"SELECT run, cartridge, {', '.join(DETECTORS)} FROM run_quality WHERE year = 1994"
        ).fetchall()
        lumi = {(run, sequence): value for run, sequence, value in db.execute("""
            SELECT run, cartridge, SUM(lumi_nb) FROM luminosity
            WHERE year = 1994 AND detector = 'STIC' GROUP BY run, cartridge
        """)}

    flags = np.asarray([row[2:] for row in quality])
    has_lumi = np.asarray([(row[0], row[1]) in lumi for row in quality])
    weights = np.asarray([lumi.get((row[0], row[1]), 0) for row in quality]) / 1000
    good = (flags >= QUALITY_MIN) & (flags <= QUALITY_MAX)
    # Columns are TPC/VD (4), MUB B/D (2), MUF A/C (2).
    tracking = good[:, :4].all(axis=1)
    barrel = good[:, 4:6].all(axis=1)
    forward = good[:, 6:8].all(axis=1)
    selections = (tracking, tracking & barrel, tracking & forward, tracking & barrel & forward)
    labels = ("Tracking\nonly", "Tracking\n+ MUB", "Tracking\n+ MUF", "Tracking\n+ MUB + MUF")
    totals = np.asarray([weights[selected].sum() for selected in selections])
    total_lumi = sum(lumi.values()) / 1000
    print("Full 1994 metadata; STIC luminosity; all required flags in [6,7].")
    print(f"No quality cut: {total_lumi:.6f} pb^-1, {len(lumi):,} luminosity keys")
    quality_keys = {(row[0], row[1]) for row in quality}
    lumi_only = lumi.keys() - quality_keys
    print(f"Lumi keys without quality: {len(lumi_only)} ({sum(lumi[key] for key in lumi_only)/1000:.6f} pb^-1); "
          f"quality keys without lumi: {len(quality_keys - lumi.keys())}")
    print(f"{'Requirement':26} {'Quality keys':>12} {'With lumi':>10} {'Lumi [pb^-1]':>14} {'Retained':>10}")
    for label, selected, value in zip(labels, selections, totals, strict=True):
        print(f"{label.replace(chr(10), ' '):26} {selected.sum():12,d} "
              f"{(selected & has_lumi).sum():10,d} {value:14.6f} {100*value/total_lumi:9.2f}%")
    print("Status 8/9 occurrences:", {name: int(np.isin(flags[:, i], [8, 9]).sum())
                                       for i, name in enumerate(DETECTORS)})

    output = PLOT_ROOT / "1994/run_quality"
    output.mkdir(parents=True, exist_ok=True)
    set_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.bar(0, total_lumi, width=0.6, color="0.7", alpha=0.75)
    ax.bar(range(1, 5), totals, width=0.6, color="C0", alpha=0.75, label="Quality 6 or 7")
    for i, value in enumerate(np.r_[total_lumi, totals]):
        ax.text(i, value - 1, f"{value:.2f}\n({100*value/total_lumi:.2f}%)",
                ha="center", va="top", fontsize=21)
    ax.set_xticks(range(5), ("No quality\ncut",) + labels, fontsize=20)
    ax.set_ylabel(r"Integrated luminosity [pb$^{-1}$]")
    ax.set_ylim(0, total_lumi * 1.4)
    ax.legend(loc="upper right", frameon=False, fontsize=20)
    mh.label.exp_label(exp="DELPHI", llabel=" Private Work", rlabel="1994", ax=ax)
    fig.tight_layout()
    fig.savefig(output / "luminosity_by_requirement.png", dpi=150)
    plt.close(fig)

    # Four disjoint categories within tracking-selected metadata.
    overlap = np.asarray([[weights[tracking & b & f].sum() for f in (~forward, forward)]
                          for b in (~barrel, barrel)])
    assert np.isclose(overlap.sum(), totals[0]), "Muon categories do not partition tracking luminosity"
    percent = 100 * overlap / totals[0]
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    image = ax.imshow(percent, cmap="Blues", vmin=0, vmax=100)
    for row, column in np.ndindex(2, 2):
        ax.text(column, row, f"{overlap[row,column]:.3f} pb$^{{-1}}$\n({percent[row,column]:.2f}%)",
                ha="center", va="center", fontsize=23,
                color="white" if percent[row,column] > 50 else "black")
    ax.set_xticks([0, 1], ["Fail", "Pass"])
    ax.set_yticks([0, 1], ["Fail", "Pass"])
    ax.set_xlabel("MUF quality")
    ax.set_ylabel("MUB quality")
    ax.tick_params(which="both", length=0)
    fig.colorbar(image, ax=ax, label="Tracking lumi [%]", fraction=0.046, pad=0.04)
    mh.label.exp_label(exp="DELPHI", llabel="  Private Work", rlabel="1994", ax=ax)
    fig.tight_layout()
    fig.savefig(output / "muon_quality_overlap.png", dpi=150)
    plt.close(fig)
    print(f"Plots: {output}")


def event_keys(directory):
    """Associate native events with raw fileSeq; no CSV or row-order join."""
    kinds = [line.split()[1] for line in (directory / "conversion.log").read_text().splitlines()
             if line.startswith("RAW_RECORD ")]
    with uproot.open(directory / "nanoaod_raw_sdst.root") as root:
        tree = root["Events"]
        assert len(kinds) == tree.num_entries, "Record log does not match raw entries"
        columns = [tree[name].array(library="np") for name in
                   ("Event_runNumber", "Event_fileSequenceNumber", "Event_eventNumber")]
    records = [(int(run), int(sequence), int(event))
               for kind, run, sequence, event in zip(kinds, *columns, strict=True) if kind == "DST"]
    lookup = {(run, event): (run, sequence) for run, sequence, event in records}
    assert len(lookup) == len(records), "Duplicate raw physics event keys"
    with uproot.open(directory / "nanoaod.root") as root:
        tree = root["Events"]
        runs = tree["Event_runNumber"].array(library="np")
        events = tree["Event_evtNumber"].array(library="np")
    native = list(zip(runs.tolist(), events.tolist(), strict=True))
    assert len(native) == len(set(native)), "Duplicate native event keys"
    assert set(native) == set(lookup), "Native/raw physics events differ"
    return [lookup[key] for key in native]


def check(input_root, output_path=None):
    database = get_db_path()
    with sqlite3.connect(f"file:{database}?mode=ro", uri=True) as db:
        quality = {(run, sequence): flags for run, sequence, *flags in db.execute(
            f"SELECT run, cartridge, {', '.join(DETECTORS)} FROM run_quality WHERE year = 1994")}
        # Same accepted-key query as Dietrich's RunSelection.
        conditions = " AND ".join(f"{name} BETWEEN {QUALITY_MIN} AND {QUALITY_MAX}" for name in DETECTORS)
        accepted = set(db.execute(
            "SELECT run, cartridge FROM run_quality WHERE year = 1994 AND " + conditions))

    print("Dietrich delphi-metadata v0.1.0; processing 94C2")
    print("All-region reference: TPC, VD, MUB and MUF must all be in [6,7].")
    print("This reference does not apply track-angle-dependent muon requirements.")
    print("No luminosity-existence cut; the default bad-HPC list applies only to 1997.")
    print(f"Metadata DB: {database.resolve()}")
    print(f"{'File':12} {'Input':>9} {'Pass':>9} {'Reject':>9} {'Pass %':>9}")
    totals = np.zeros(6, dtype=int)
    seen, kept = Counter(), Counter()
    for path in sorted(input_root.glob("*/nanoaod.root")):
        directory = path.parent
        keys = event_keys(directory)
        seen.update(keys)
        kept.update(key for key in keys if key in accepted)
        selected = np.asarray([key in accepted for key in keys])
        known = np.asarray([key in quality for key in keys])
        flags = np.asarray([quality.get(key, [-1] * len(DETECTORS)) for key in keys])
        remaining = known.copy()
        cutflow = [len(keys), int(remaining.sum())]
        for offset in range(0, len(DETECTORS), 2):
            remaining &= ((flags[:, offset:offset + 2] >= QUALITY_MIN)
                          & (flags[:, offset:offset + 2] <= QUALITY_MAX)).all(axis=1)
            cutflow.append(int(remaining.sum()))
        assert np.array_equal(remaining, selected), "Cutflow differs from Dietrich's SQL selection"
        totals += cutflow
        passed = int(selected.sum())
        print(f"{directory.name:12} {len(keys):9,d} {passed:9,d} {len(keys)-passed:9,d} {100*selected.mean():8.2f}%")
    assert totals[0] > 0, f"No collision events found under {input_root}"
    print(f"{'TOTAL':12} {totals[0]:9,d} {totals[-1]:9,d} {totals[0]-totals[-1]:9,d} {100*totals[-1]/totals[0]:8.2f}%")
    print("Sequential event cutflow:")
    for label, count in zip(("Input", "Known quality key", "TPC", "VD", "MUB", "MUF"), totals, strict=True):
        print(f"  {label:18}: {count:,}")
    all_keys = set(seen)
    missing = all_keys - quality.keys()
    print(f"Observed run/fileSeq keys: {len(all_keys)}; passing: {len(all_keys & accepted)}; missing: {len(missing)}")
    if missing:
        print("Missing quality keys:", sorted(missing))

    # Dietrich's writeReport schema, directly usable by runselection_summary.py.
    report = {}
    for (run, sequence), count in sorted(seen.items()):
        entry = report.setdefault(str(run), {
            "events_seen": 0, "events_accepted": 0, "cartridges": {},
        })
        entry["events_seen"] += count
        entry["events_accepted"] += kept[run, sequence]
        entry["cartridges"][str(sequence)] = {
            "events_seen": count, "events_accepted": kept[run, sequence],
        }
    output_path = output_path or DATASET_ROOT / input_root.name / "run_quality_runs.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Run report: {output_path}")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true",
                        help="Plot full-1994 metadata comparisons without reading collision files")
    parser.add_argument("--input", type=Path, default=DATA_ROOT / "20260828_florian/20261026_short94_c2")
    parser.add_argument("-o", "--output", type=Path,
                        help="default: output/dataset/<input>/run_quality_runs.json")
    args = parser.parse_args()
    if args.plot:
        plot_quality()
    else:
        check(args.input, args.output)


if __name__ == "__main__":
    main()
