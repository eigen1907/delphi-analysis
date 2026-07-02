from __future__ import annotations

from pathlib import Path

import hist
import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
import uproot

from plot_utils import (
    FILE_GROUPS,
    TREE_NAMES_BY_SOURCE,
    add_delphi_label,
    build_branch_summary,
    histogram_bins,
    read_numeric_values,
    resolve_file,
    resolve_samples,
    sample_styles,
    set_hist_yaxis,
    tree_branch_plot_path,
)


def plot_branches_all(input_root: Path, output_root: Path, samples: list[str] | tuple[str, ...] | None = None) -> None:
    mh.style.use(mh.styles.CMS)
    samples = resolve_samples(input_root, samples)
    styles = sample_styles(samples)
    output_root = output_root / "branches_all"
    summaries = {sample: build_branch_summary(input_root / sample, sample) for sample in samples}

    for source, file_names in FILE_GROUPS:
        for tree_name in TREE_NAMES_BY_SOURCE.get(source, ("Events",)):
            source_branches = set()
            for sample in samples:
                branches = summaries[sample].get("trees", {}).get(tree_name, {}).get("branches", {})
                for branch, branch_info in branches.items():
                    if source in branch_info["source"]:
                        source_branches.add(branch)

            if not source_branches:
                continue

            trees = {}
            typenames = {}
            for sample in samples:
                try:
                    path = resolve_file(input_root / sample, file_names)
                except FileNotFoundError:
                    continue
                root_file = uproot.open(path)
                if tree_name in root_file:
                    trees[sample] = root_file[tree_name]
                    typenames[sample] = trees[sample].typenames()

            n_plots = 0
            for branch in sorted(source_branches):
                values_by_sample = {}
                stats_by_sample = {}
                for sample in samples:
                    if sample not in trees or branch not in typenames[sample]:
                        continue
                    dtype = typenames[sample][branch]
                    values = read_numeric_values(trees[sample], branch, dtype)
                    if not values.size:
                        continue
                    values_by_sample[sample] = values
                    stats_by_sample[sample] = {
                        "dtype": dtype,
                        "n_event": int(trees[sample].num_entries),
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

                if not stats_by_sample:
                    continue

                nbins, lo, hi = histogram_bins(stats_by_sample)

                fig, ax = plt.subplots(figsize=(12, 10))
                legend_header = f"{'Sample':<8} {'Event':<8} {'Object':<8} {'UF':<6} {'OF':<6}"
                ax.plot([], [], color="none", label=legend_header)
                bin_edges = np.linspace(lo, hi, nbins + 1)
                counts_list = []
                density_list = []
                for sample, branch_values in values_by_sample.items():
                    h = hist.Hist.new.Reg(nbins, lo, hi, name="value", label=branch).Double()
                    h.fill(value=branch_values)
                    stats = stats_by_sample[sample]
                    underflow = int(np.count_nonzero(branch_values < lo))
                    overflow = int(np.count_nonzero(branch_values > hi))
                    counts, _ = np.histogram(branch_values, bins=bin_edges)
                    counts_list.append(counts)
                    density_list.append(counts / (np.sum(counts) * np.diff(bin_edges)) if np.sum(counts) else counts)
                    label = (
                        f"{sample:<8} "
                        f"{stats['n_event']:<8} "
                        f"{stats['n_object']:<8} "
                        f"{underflow:<6} "
                        f"{overflow:<6}"
                    )
                    color, hatch = styles[sample]
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
                ax.legend(loc="upper right", prop={"family": "monospace", "size": 18})
                ax.set_xlabel(branch)
                ax.set_ylabel("Normalized")
                add_delphi_label(ax)
                fig.tight_layout()

                plot_path = tree_branch_plot_path(output_root / source, tree_name, branch)
                plot_path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(plot_path, dpi=150)
                plt.close(fig)
                n_plots += 1

            print(f"[{source}/{tree_name}] plotted {n_plots} branches under {output_root / source / tree_name}")
