from __future__ import annotations

import argparse
from common.dataset import SOURCE_TREES, dataset_files
from common.paths import DATASET_ROOT, PROJECT_ROOT
from common.plot_utils import default_plot_root, add_samples_argument
from pathlib import Path

import hist
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from common.plot_utils import (
    SOURCE_STYLES,
    add_delphi_label,
    build_branch_summary,
    histogram_bins,
    read_numeric_values,
    set_hist_yaxis,
    tree_branch_plot_path,
)


def plot_branches_compare(input_root: Path, output_root: Path, samples: list[str] | tuple[str, ...] | None = None, split: str = "train") -> None:
    mh.style.use(mh.styles.CMS)
    files = dict(dataset_files(input_root, split, samples))
    samples = tuple(files)
    output_dir = output_root / "00_overall" / "01_source_comparison"

    for sample in samples:
        summary = build_branch_summary(files[sample][0].parent, sample, split)
        trees = {}

        root_file = uproot.open(files[sample][0])
        for source, tree_name in SOURCE_TREES.items():
            trees[source] = root_file[tree_name]

        compare_branches = {}
        event_branches = summary.get("trees", {}).get("Events", {}).get("branches", {})
        for branch, branch_info in event_branches.items():
            sources = [
                source
                for source in SOURCE_TREES
                if source in branch_info["source"] and source in trees
            ]
            if len(sources) >= 2:
                compare_branches[branch] = sources

        n_plots = 0
        for branch, sources in sorted(compare_branches.items()):
            values_by_source = {}
            stats_by_source = {}
            for source in sources:
                dtype = trees[source].typenames()[branch]
                values = read_numeric_values(trees[source], branch, dtype)
                if not values.size:
                    continue
                values_by_source[source] = values
                stats_by_source[source] = {
                    "dtype": dtype,
                    "n_event": int(trees[source].num_entries),
                    "n_object": int(values.size),
                    "min": float(np.min(values)),
                    "q005": float(np.quantile(values, 0.005)),
                    "q25": float(np.quantile(values, 0.25)),
                    "q75": float(np.quantile(values, 0.75)),
                    "q995": float(np.quantile(values, 0.995)),
                    "max": float(np.max(values)),
                    "mean": float(np.mean(values)),
                    "std": float(np.std(values)),
                }

            if len(stats_by_source) < 2:
                continue

            nbins, lo, hi = histogram_bins(stats_by_source)

            fig, ax = plt.subplots(figsize=(12, 10))
            legend_header = f"{'Source':<18} {'Event':<8} {'Object':<8} {'UF':<6} {'OF':<6}"
            ax.plot([], [], color="none", label=legend_header)
            bin_edges = np.linspace(lo, hi, nbins + 1)
            counts_list = []
            density_list = []
            for source, branch_values in values_by_source.items():
                h = hist.Hist.new.Reg(nbins, lo, hi, name="value", label=branch).Double()
                h.fill(value=branch_values)
                stats = stats_by_source[source]
                underflow = int(np.count_nonzero(branch_values < lo))
                overflow = int(np.count_nonzero(branch_values > hi))
                counts, _ = np.histogram(branch_values, bins=bin_edges)
                counts_list.append(counts)
                density_list.append(counts / (np.sum(counts) * np.diff(bin_edges)) if np.sum(counts) else counts)
                label = (
                    f"{source:<18} "
                    f"{stats['n_event']:<8} "
                    f"{stats['n_object']:<8} "
                    f"{underflow:<6} "
                    f"{overflow:<6}"
                )
                color, hatch = SOURCE_STYLES[source]
                mh.histplot(
                    h,
                    ax=ax,
                    histtype="fill",
                    density=True,
                    label=label,
                    facecolor="none",
                    edgecolor=color,
                    hatch=hatch,
                    linewidth=0.0,
                )
                mh.histplot(
                    h,
                    ax=ax,
                    histtype="step",
                    density=True,
                    label="_nolegend_",
                    yerr=True,
                    color=color,
                    alpha=0.9,
                )

            set_hist_yaxis(ax, counts_list, density_list)
            ax.legend(loc="upper right", prop={"family": "monospace", "size": 16})
            ax.set_xlabel(branch)
            ax.set_ylabel("Normalized")
            add_delphi_label(ax)
            fig.tight_layout()

            plot_path = tree_branch_plot_path(output_dir / sample, "Events", branch)
            plot_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(plot_path, dpi=150)
            plt.close(fig)
            n_plots += 1

        print(f"[{sample}] plotted {n_plots} compare Events branches under {output_dir / sample}")


def main() -> None:
    parser = argparse.ArgumentParser()
    add_samples_argument(parser)
    parser.add_argument("-i", "--input", type=Path,
                        default=DATASET_ROOT / "20261015_jongwon",
                        help="input dataset root")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="plot output root (default: output/plots/<campaign>/03_sample_validation/<split>)")
    parser.add_argument("--split", choices=("train", "val", "test"), default="train")
    args = parser.parse_args()

    input_root = args.input
    output_root = args.output or default_plot_root(PROJECT_ROOT, input_root, args.split)
    plot_branches_compare(input_root, output_root, args.samples, args.split)


if __name__ == "__main__":
    main()
