from __future__ import annotations

from pathlib import Path

import hist
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from plot_utils import (
    COMPARE_FILE_GROUPS,
    SAMPLE_LABELS,
    SOURCE_STYLES,
    add_delphi_label,
    histogram_bins,
    load_branch_summary,
    read_numeric_values,
    resolve_file,
    set_hist_yaxis,
    tree_branch_plot_path,
)


def plot_branches_compare(input_root: Path, output_root: Path, branch_root: Path) -> None:
    mh.style.use(mh.styles.CMS)
    output_dir = output_root / "branches_compare"

    for sample in SAMPLE_LABELS:
        summary = load_branch_summary(branch_root, sample)
        sample_dir = input_root / sample
        trees = {}

        for source, file_names in COMPARE_FILE_GROUPS:
            path = resolve_file(sample_dir, file_names)
            trees[source] = uproot.open(path)["Events"]

        compare_branches = {}
        for branch, branch_info in summary["trees"]["Events"]["branches"].items():
            sources = [
                source
                for source, _ in COMPARE_FILE_GROUPS
                if source in branch_info["source"]
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
