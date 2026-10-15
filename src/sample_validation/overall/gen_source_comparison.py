from __future__ import annotations

import argparse
from common.dataset import SOURCE_TREES, dataset_files
from common.paths import DATASET_ROOT, PROJECT_ROOT
from common.plot_utils import default_plot_root, add_samples_argument, default_check_root
import csv
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from common.plot_utils import GEN_COMPARE_PAIR_STYLES, add_delphi_label, set_hist_yaxis


SOURCES = (
    ("nanoaod", "Event_evtNumber"),
    ("nanoaod_raw_sdst", "Event_eventNumber"),
)
PAIRS = (("nanoaod", "nanoaod_raw_sdst"),)

FIELDS = (
    ("nGenPart", r"N_{\mathrm{GenPart}}", (-0.5, 0.5), "nGenPart"),
    ("status", r"\mathrm{status}", (-0.5, 0.5), "GenPart_status"),
    ("pdgId", r"\mathrm{PDG~ID}", (-0.5, 0.5), "GenPart_pdgId"),
    ("parentIdx", r"\mathrm{parent~index}", (-0.5, 0.5), "GenPart_parentIdx"),
    ("firstChildIdx", r"\mathrm{first~child~index}", (-0.5, 0.5), "GenPart_firstChildIdx"),
    ("lastChildIdx", r"\mathrm{last~child~index}", (-0.5, 0.5), "GenPart_lastChildIdx"),
    ("mass", r"m_{\mathrm{GenPart}}", (-1e-6, 1e-6), "GenPart_mass"),
    ("vertex_x", r"x_{\mathrm{vtx}}", (-0.2, 0.2), "GenPart_vertex.fCoordinates.fX"),
    ("vertex_y", r"y_{\mathrm{vtx}}", (-0.002, 0.002), "GenPart_vertex.fCoordinates.fY"),
    ("vertex_z", r"z_{\mathrm{vtx}}", (-2.5, 2.5), "GenPart_vertex.fCoordinates.fZ"),
    (
        "vertex_t",
        r"t_{\mathrm{vtx}}",
        (-1e-6, 1e-6),
        {
            "nanoaod": "GenPart_vertex.fCoordinates.fT",
            "nanoaod_raw_sdst": "GenPart_productionTime",
            "nanoaod_raw_fadana": "GenPart_productionTime",
        },
    ),
    (
        "p4_px",
        r"p_x",
        (-1e-6, 1e-6),
        {
            "nanoaod": "GenPart_vector.fCoordinates.fX",
            "nanoaod_raw_sdst": "GenPart_fourMomentum.fCoordinates.fX",
            "nanoaod_raw_fadana": "GenPart_fourMomentum.fCoordinates.fX",
        },
    ),
    (
        "p4_py",
        r"p_y",
        (-1e-6, 1e-6),
        {
            "nanoaod": "GenPart_vector.fCoordinates.fY",
            "nanoaod_raw_sdst": "GenPart_fourMomentum.fCoordinates.fY",
            "nanoaod_raw_fadana": "GenPart_fourMomentum.fCoordinates.fY",
        },
    ),
    (
        "p4_pz",
        r"p_z",
        (-1e-6, 1e-6),
        {
            "nanoaod": "GenPart_vector.fCoordinates.fZ",
            "nanoaod_raw_sdst": "GenPart_fourMomentum.fCoordinates.fZ",
            "nanoaod_raw_fadana": "GenPart_fourMomentum.fCoordinates.fZ",
        },
    ),
    (
        "p4_e",
        r"E",
        (-1e-6, 1e-6),
        {
            "nanoaod": "GenPart_vector.fCoordinates.fT",
            "nanoaod_raw_sdst": "GenPart_fourMomentum.fCoordinates.fT",
            "nanoaod_raw_fadana": "GenPart_fourMomentum.fCoordinates.fT",
        },
    ),
    (
        "lifetime",
        r"\tau",
        (-1e-6, 1e-6),
        {
            "nanoaod": "GenPart_tau",
            "nanoaod_raw_sdst": "GenPart_properLifetime",
            "nanoaod_raw_fadana": "GenPart_properLifetime",
        },
    ),
)


def branch_name(branches, source: str) -> str:
    return branches[source] if isinstance(branches, dict) else branches


def as_list(value) -> list:
    return value if isinstance(value, list) else [value]


def plot_gen_compare(input_root: Path, output_root: Path, data_root: Path, samples: list[str] | tuple[str, ...] | None = None, split: str = "train") -> None:
    mh.style.use(mh.styles.CMS)
    files = dict(dataset_files(input_root, split, samples))
    samples = tuple(files)

    output_root = output_root / "00_overall" / "03_gen_source_comparison"
    data_root = data_root / "gen_compare"
    data_root.mkdir(parents=True, exist_ok=True)
    diff_rows = []

    for sample in samples:
        sample_output = output_root / sample
        sample_output.mkdir(parents=True, exist_ok=True)

        data = {}
        for source, event_branch in SOURCES:
            tree = uproot.open(files[sample][0])[SOURCE_TREES[source]]
            branches = ["job_id", "Event_runNumber", event_branch]
            branches += [branch_name(field_branches, source) for _, _, _, field_branches in FIELDS]
            arrays = {branch: tree[branch].array(library="ak") for branch in dict.fromkeys(branches)}
            keys = list(zip(ak.to_numpy(arrays["job_id"]).tolist(), ak.to_numpy(arrays["Event_runNumber"]).tolist(), ak.to_numpy(arrays[event_branch]).tolist(), strict=True))
            if len(set(keys)) != len(keys):
                raise ValueError(f"{sample} {source}: duplicate job/run/event keys")
            data[source] = {
                "arrays": arrays,
                "keys": keys,
                "index": {key: i for i, key in enumerate(keys)},
            }

        for field_name, x_label, (lo, hi), field_branches in FIELDS:
            values_by_pair = {}

            for left, right in PAIRS:
                left_data = data[left]
                right_data = data[right]
                keys = [key for key in left_data["keys"] if key in right_data["index"]]
                left_index = np.asarray([left_data["index"][key] for key in keys])
                right_index = np.asarray([right_data["index"][key] for key in keys])
                left_values = left_data["arrays"][branch_name(field_branches, left)][left_index]
                right_values = right_data["arrays"][branch_name(field_branches, right)][right_index]

                pair = f"{left} - {right}"
                pair_residuals = []
                pair_diff_rows = []

                for (job_id, run, event), left_event, right_event in zip(keys, ak.to_list(left_values), ak.to_list(right_values), strict=True):
                    left_list = as_list(left_event)
                    right_list = as_list(right_event)
                    residuals = [float(a) - float(b) for a, b in zip(left_list, right_list)]
                    n_nonzero = sum(abs(value) > 0 for value in residuals)
                    n_nonzero += abs(len(left_list) - len(right_list))
                    pair_residuals += residuals

                    if n_nonzero:
                        pair_diff_rows.append(
                            {
                                "sample": sample,
                                "field": field_name,
                                "pair": pair,
                                "job_id": int(job_id),
                                "run": int(run),
                                "event": int(event),
                                "n_left": len(left_list),
                                "n_right": len(right_list),
                                "n_nonzero": n_nonzero,
                                "max_abs": max([abs(value) for value in residuals], default=0.0),
                            }
                        )

                if pair_diff_rows:
                    events = ", ".join(str(row["event"]) for row in pair_diff_rows[:10])
                    suffix = " ..." if len(pair_diff_rows) > 10 else ""
                    print(f"[{sample}] {field_name} {pair}: {len(pair_diff_rows)} events differ ({events}{suffix})")

                values_by_pair[(left, right)] = np.asarray(pair_residuals, dtype=float)
                diff_rows += pair_diff_rows

            max_abs = max((float(np.max(np.abs(values))) for values in values_by_pair.values() if len(values)), default=0.0)
            if max_abs > max(abs(lo), abs(hi)):
                lo, hi = -1.05 * max_abs, 1.05 * max_abs
            bins = np.linspace(lo, hi, 101)

            scale = 1.0
            scale_label = ""
            width = float(np.max(np.abs(bins)))
            if width < 1e-2 or width >= 1e3:
                power = int(np.floor(np.log10(width)))
                scale = 10.0**power
                scale_label = rf"/10^{{{power}}}"

            fig, ax = plt.subplots(figsize=(12, 10))
            ax.plot([], [], color="none", label=f"{'Pair':<36} {'Entries':<7} {'Nonzero':<3}")

            counts_list = []
            for pair, values in values_by_pair.items():
                nonzero = int(np.count_nonzero(np.abs(values) > 0))
                pair_label = f"{pair[0]}, {pair[1]}"
                label = f"{pair_label:<36} {len(values):<7} {nonzero:<3}"
                color, hatch = GEN_COMPARE_PAIR_STYLES[pair]
                counts, _ = np.histogram(values / scale, bins=bins / scale)
                counts_list.append(counts)
                ax.hist(
                    values / scale,
                    bins=bins / scale,
                    histtype="stepfilled",
                    label=label,
                    facecolor="none",
                    edgecolor=color,
                    hatch=hatch,
                    linewidth=0.0,
                )
                ax.hist(
                    values / scale,
                    bins=bins / scale,
                    histtype="step",
                    label="_nolegend_",
                    color=color,
                    alpha=0.9,
                    linewidth=2,
                )

            set_hist_yaxis(ax, counts_list)
            ax.set_xlabel(rf"$\Delta {x_label}{scale_label}$")
            ax.set_ylabel("Entries")
            ax.ticklabel_format(axis="x", style="plain", useOffset=False)
            ax.legend(loc="upper right", prop={"family": "monospace", "size": 16})
            ax.grid(alpha=0.3)
            add_delphi_label(ax)
            fig.tight_layout()
            fig.savefig(sample_output / f"{field_name}.png", dpi=150)
            plt.close(fig)

        print(f"[{sample}] wrote {sample_output}")

    diff_fields = ("sample", "field", "pair", "job_id", "run", "event", "n_left", "n_right", "n_nonzero", "max_abs")
    with (data_root / "diff-events.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=diff_fields)
        writer.writeheader()
        writer.writerows(diff_rows)

    print(f"wrote {data_root / 'diff-events.csv'}")


def main() -> None:
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", type=Path,
                        default=DATASET_ROOT / "20261015_jongwon",
                        help="input dataset root")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="plot output root (default: output/plots/<campaign>/03_sample_validation/<split>)")
    parser.add_argument("--data-root", type=Path, help="diagnostic data root")
    parser.add_argument("--split", choices=("train", "val", "test"), default="train")
    args = parser.parse_args()

    input_root = args.input
    output_root = args.output or default_plot_root(PROJECT_ROOT, input_root, args.split)
    data_root = args.data_root or default_check_root(PROJECT_ROOT, input_root) / args.split
    plot_gen_compare(input_root, output_root, data_root, args.samples, args.split)


if __name__ == "__main__":
    main()
