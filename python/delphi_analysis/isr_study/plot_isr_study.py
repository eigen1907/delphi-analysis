"""Read the five samples, then draw distributions, link checks, and efficiencies."""
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.stats import beta

from .angle_matching import accumulate, angular_stats, choose_cut, new_accumulator, plot_angle_matching
from .data import ANALYSIS_SELECTIONS, GEN_LABELS, GEN_SELECTIONS, LINK_STATUSES, RECO_LABELS, SAMPLES, read_sample

FIGURE_SIZE = (11, 9)
COVERAGE = 0.6826894921370859
ENERGY_EDGES = np.array([0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50])
MIN_GEN_ENERGY = 0.1
GEN_COMPONENTS = (('beam_isr', r'$\gamma_{\mathrm{BeamISR}}$'),
                  ('nonbeam_isr', r'$\gamma_{\mathrm{NonBeamISR}}$'),
                  ('fsr', r'$\gamma_{\mathrm{FSR}}$'), ('decayed', r'$\gamma_{\mathrm{Decayed}}$'))
GEN_STACKS = dict(gen_gamma=GEN_COMPONENTS, gen_isr=GEN_COMPONENTS[:2],
                  gen_no_isr=GEN_COMPONENTS[2:], gen_isr_non_beam=GEN_COMPONENTS[1:2])
XLABELS = dict(energy=r'$E_\gamma^{\rm gen}$ [GeV]', cos_theta=r'$\cos\theta^{\rm gen}$')
SAMPLE_LABELS = dict(Zee=r'$\boldsymbol{Z\to e^+e^-}$', Zmumu=r'$\boldsymbol{Z\to\mu^+\mu^-}$',
                     Ztautau=r'$\boldsymbol{Z\to\tau^+\tau^-}$', ZKK=r'$\boldsymbol{Z\to K^+K^-}$',
                     Zpipi=r'$\boldsymbol{Z\to\pi^+\pi^-}$', combined='All samples')


def clopper_pearson(numerator, denominator):
    alpha = 1 - COVERAGE
    lower, upper = np.zeros_like(denominator, dtype=float), np.ones_like(denominator, dtype=float)
    nonzero, not_all = numerator > 0, numerator < denominator
    lower[nonzero] = beta.ppf(alpha / 2, numerator[nonzero], denominator[nonzero] - numerator[nonzero] + 1)
    upper[not_all] = beta.ppf(1 - alpha / 2, numerator[not_all] + 1, denominator[not_all] - numerator[not_all])
    return lower, upper


def finish(fig, ax, path, sample, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ax.get_legend() is None:
        ax.legend(loc='upper right', frameon=False)
    ax.text(0.02, 0.97, SAMPLE_LABELS[sample], transform=ax.transAxes,
            va='top', fontweight='bold')
    ax.grid(alpha=0.2)
    mh.label.exp_label(exp='DELPHI', llabel='Simulation', rlabel='LEP 1 (91.2 GeV)', loc=0, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def energy_bins(entries, logarithmic=False):
    maximum = np.nextafter(entries.max(), np.inf)
    if logarithmic:
        positive = entries[entries > 0]
        low = positive.min() / 2
        count = max(1, int(np.ceil(6 * np.log10(maximum / low))))
        return np.r_[0, np.geomspace(low, maximum, count + 1)]
    return np.r_[ENERGY_EDGES[ENERGY_EDGES < maximum], maximum]


def bin_centers(ax, bins, logarithmic=False):
    ax.set_xlim(bins[0], bins[-1])
    if logarithmic:
        # The first interval contains zero only, rather than dropping zero-energy candidates.
        ax.set_xscale('symlog', linthresh=bins[1])
        centers = np.r_[0, np.sqrt(bins[1:-1] * bins[2:])]
        # Space the labels around zero; the bins and all entries are unchanged.
        ax.set_xticks(np.r_[0, ax.get_xticks()[2:]])
    else:
        centers = (bins[:-1] + bins[1:]) / 2
    return centers


def gen_spectrum(ax, entries, category, components, bins, events, logarithmic):
    centers = bin_centers(ax, bins, logarithmic)
    bottom = np.zeros(len(bins) - 1)
    for key, label in components:
        index = GEN_COMPONENTS.index((key, label))
        selected = entries[(category == key) & np.isfinite(entries)]
        counts = np.histogram(selected, bins=bins)[0]
        assert counts.sum() == len(selected), 'Spectrum must retain all finite entries'
        ax.stairs(bottom + counts / events, bins, baseline=bottom, fill=True,
                  label=f'Gen {label}', color=f'C{index}', linewidth=0.7)
        bottom += counts / events
    total = np.histogram(entries[np.isfinite(entries)], bins=bins)[0]
    assert np.allclose(bottom * events, total)
    shown = total > 0
    ax.errorbar(centers[shown], bottom[shown], yerr=np.sqrt(total[shown]) / events,
                fmt='none', ecolor='grey', capsize=2, elinewidth=0.8)


def reco_spectrum(ax, entries, bins, events, logarithmic, label):
    """One inclusive Reco curve, with no truth-origin selection."""
    centers = bin_centers(ax, bins, logarithmic)
    entries = entries[np.isfinite(entries)]
    counts = np.histogram(entries, bins=bins)[0]
    assert counts.sum() == len(entries)
    ax.stairs(counts / events, bins, fill=True, facecolor='none', edgecolor='black',
              hatch='///', linewidth=1.8, label=f'Reco {label}')
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / events, yerr=np.sqrt(counts[shown]) / events,
                fmt='none', ecolor='black', capsize=2, elinewidth=0.8)
    return counts / events


def sample_distributions(study, sample, values, bins):
    """Change only the Gen selection; the inclusive Reco curve stays the same."""
    gamma = values['gamma']
    candidates = values['reco_distributions']
    for gen_name in GEN_SELECTIONS:
        output = study / '01_sample_distribution' / gen_name
        output.mkdir(parents=True, exist_ok=True)
        gen_selected = values['selections'][gen_name]
        for coordinate in ('cos_theta', 'energy'):
            for log in (False, True):
                log_energy = log and coordinate == 'energy'
                edges = bins['energy_log'] if log_energy else bins[coordinate]
                gen_entries = gamma[coordinate]
                reco_entries = candidates[coordinate]
                if coordinate == 'energy' and not log:
                    # E >= 50 GeV goes into the labelled overflow bin, never out of the study.
                    last = np.nextafter(edges[-1], -np.inf)
                    gen_entries = np.minimum(gen_entries, last)
                    reco_entries = np.minimum(reco_entries, last)
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                gen_spectrum(ax, gen_entries[gen_selected], gamma['category'][gen_selected],
                             GEN_STACKS[gen_name], edges, values['events'], log_energy)
                reco_height = reco_spectrum(ax, reco_entries, edges, values['events'], log_energy, values['reco_label'])
                all_gen_counts = np.histogram(gen_entries[np.isfinite(gen_entries)], bins=edges)[0]
                # Use one y range for all four Gen selections.
                height = max(all_gen_counts.max() / values['events'], reco_height.max(), 1 / values['events'])
                if log:
                    ax.set_yscale('log')
                    ax.set_ylim(0.5 / values['events'], height * 1000)
                else:
                    ax.set_ylim(0, height * 1.8)
                if coordinate == 'energy' and not log:
                    ax.set_xticks([0, 10, 20, 30, 40, 55], ['0', '10', '20', '30', '40', r'$\geq 50$'])
                suffix = '_log' if log else ''
                xlabel = r'$E_\gamma$ [GeV]' if coordinate == 'energy' else r'$\cos\theta$'
                finish(fig, ax, output / f'{coordinate}_{sample}{suffix}.png', sample, xlabel,
                       r'$N_\gamma$ / bin per event')
        event_multiplicity(output, sample, values, GEN_STACKS[gen_name])


def event_multiplicity(output, sample, values, components):
    """Photon counts per event, including N_gamma = 0; overlapping curves are not stacked."""
    entries = [(values['gen_event_counts'][name], f'Gen {label}', f'C{GEN_COMPONENTS.index((name, label))}', '-')
               for name, label in components]
    entries.append((values['reco_event_counts'], f"Reco {values['reco_label']}", 'C4', '--'))
    bins = np.arange(-0.5, max(item[0].max() for item in entries) + 1.5)
    centers = (bins[:-1] + bins[1:]) / 2
    for log in (False, True):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        height = 0
        for numbers, label, color, line in entries:
            counts = np.histogram(numbers, bins=bins)[0]
            assert counts.sum() == len(numbers) == values['events']
            shown = counts > 0
            ax.stairs(counts, bins, label=label, color=color, linestyle=line)
            ax.errorbar(centers[shown], counts[shown], yerr=np.sqrt(counts[shown]),
                        fmt='none', color=color, capsize=2)
            height = max(height, (counts + np.sqrt(counts)).max())
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlim(bins[0], bins[-1])
        if log:
            ax.set_yscale('log')
            ax.set_ylim(0.5, height * 100)
        else:
            ax.set_ylim(0, height * 1.8)
        suffix = '_log' if log else ''
        finish(fig, ax, output / f'multiplicity_{sample}{suffix}.png',
               sample, r'$N_\gamma$', 'Events')


def truth_link_validation(study, sample, values):
    labels = ('Agree', 'Forward\nonly', 'Reverse\nonly', 'Conflict', 'No origin')
    output = study / '02_truth_link_matching_validation'
    output.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    x = np.arange(len(labels))
    maximum = 0
    for index, (gen_name, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
        counts = np.array([values['matching_validation'][gen_name][key] for key in LINK_STATUSES])
        positions = x + (index - 1) * 0.25
        ax.bar(positions, counts, width=0.25, color=f'C{index}', label=label)
        ax.errorbar(positions, counts, yerr=np.sqrt(counts), fmt='none', color=f'C{index}', capsize=2)
        maximum = max(maximum, counts.max())
    ax.set_xticks(x, labels)
    ax.set_yscale('log')
    ax.set_ylim(top=max(maximum * 100, 1))
    finish(fig, ax, output / f'truth_link_status_{sample}.png', sample,
           'Link status', f"# of linked {values['reco_label']}")


def profile_points(coordinate, response, bins):
    defined = np.isfinite(response)
    count = np.histogram(coordinate[defined], bins=bins)[0]
    sums = np.histogram(coordinate[defined], bins=bins, weights=response[defined])[0]
    squares = np.histogram(coordinate[defined], bins=bins, weights=response[defined] ** 2)[0]
    assert count.sum() == defined.sum()
    mean = np.divide(sums, count, out=np.full(len(count), np.nan), where=count > 0)
    error = np.full(len(count), np.nan)
    measured = count > 1
    error[measured] = np.sqrt(np.maximum(squares[measured] - sums[measured] ** 2 / count[measured], 0)
                              / (count[measured] * (count[measured] - 1)))
    return mean, error


def matched_energy_distributions(output, sample, values, stats):
    """Dimensionless loss and ratio for matched Gen photons above the energy cut."""
    populations = []
    for population, label in zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True):
        selected = (values['analysis_selections'][population]
                    & (values['gamma']['energy'] >= MIN_GEN_ENERGY) & (stats['reco_count'] > 0))
        defined = selected & np.isfinite(stats['energy_ratio'])
        populations.append((stats['energy_ratio'][defined], label))
        print(f'{sample} {output.parents[1].name}/{output.name} {population}: '
              f'matched E>={MIN_GEN_ENERGY}={selected.sum()}, undefined energy={(selected & ~defined).sum()}', flush=True)
    for name, bins, xlabel, tail_title in (
        ('residual', np.linspace(-2, 2, 41),
         r'$(E_\gamma^{\rm gen}-\sum E^{\rm reco})/E_\gamma^{\rm gen}$', 'Frac. (<−2)'),
        ('ratio', np.linspace(0, 3, 31),
         r'$\sum E^{\rm reco}/E_\gamma^{\rm gen}$', 'Frac. (>3)'),
    ):
        centers = (bins[:-1] + bins[1:]) / 2
        histograms = []
        for ratio, label in populations:
            entries = 1 - ratio if name == 'residual' else ratio
            in_range = entries[(entries >= bins[0]) & (entries <= bins[-1])]
            counts = np.histogram(in_range, bins=bins)[0]
            tail = len(entries) - len(in_range)
            assert counts.sum() + tail == len(entries)
            histograms.append((counts, label, in_range.mean(), tail / len(entries)))
            print(f'{sample} {output.parents[1].name}/{output.name} {label} {name}: '
                  f'outside [{bins[0]:g}, {bins[-1]:g}]={tail}/{len(entries)}', flush=True)
        maximum = max(counts.max() for counts, _, _, _ in histograms)
        for log in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for index, (counts, label, _, _) in enumerate(histograms):
                color = f'C{index}'
                ax.stairs(counts, bins, label=label, color=color)
                shown = counts > 0
                ax.errorbar(centers[shown], counts[shown], yerr=np.sqrt(counts[shown]),
                            fmt='none', color=color, capsize=2)
            # Three single-column legends; statistics use empty handles.
            groups = ax.legend(title='Gen', loc='upper right', bbox_to_anchor=(0.46, 0.96),
                               frameon=False, handleheight=1.6)
            ax.add_artist(groups)
            empty = [Line2D([], [], color='none') for _ in histograms]
            means = ax.legend(empty, [f'{mean:.2f}' for _, _, mean, _ in histograms],
                              title='Mean (shown)', loc='upper right', bbox_to_anchor=(0.75, 0.96),
                              frameon=False, handlelength=0, handletextpad=0, handleheight=1.6)
            ax.add_artist(means)
            ax.legend(empty, [f'{fraction:.4f}' for _, _, _, fraction in histograms],
                      title=tail_title, loc='upper right', bbox_to_anchor=(0.99, 0.96),
                      frameon=False, handlelength=0, handletextpad=0, handleheight=1.6)
            ax.text(0.02, 0.89, rf'$E_\gamma^{{\rm gen}}\geq{MIN_GEN_ENERGY}$ GeV',
                    transform=ax.transAxes, va='top')
            ax.set_xlim(bins[0], bins[-1])
            ax.set_xticks(np.arange(bins[0], bins[-1] + 0.5, 0.5))
            if log:
                ax.set_yscale('log')
                ax.set_ylim(0.5, 0.5 * (maximum / 0.5) ** 1.65)
            else:
                ax.set_ylim(0, maximum * 1.65)
            suffix = '_log' if log else ''
            finish(fig, ax, output / f'matched_energy_{name}_{sample}{suffix}.png',
                   sample, xlabel, 'Gen photons')


def gen_reco_energy_histograms(output, sample, values, stats):
    """One matched Gen photon per entry; sum every associated Reco energy."""
    gen_energy = values['gamma']['energy']
    gen_edges = np.r_[MIN_GEN_ENERGY, ENERGY_EDGES[ENERGY_EDGES > MIN_GEN_ENERGY]]
    reco_edges = np.r_[ENERGY_EDGES, 60]
    for population, label in zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True):
        selected = (values['analysis_selections'][population] & (gen_energy >= MIN_GEN_ENERGY)
                    & (stats['reco_count'] > 0) & np.isfinite(stats['energy_ratio']))
        energy = stats['energy_ratio'][selected] * gen_energy[selected]
        # Keep the full Reco tail in the labelled E >= 50 GeV bin.
        energy = np.minimum(energy, np.nextafter(reco_edges[-1], -np.inf))
        counts = np.histogram2d(gen_energy[selected], energy, bins=(gen_edges, reco_edges))[0]
        assert counts.sum() == selected.sum()
        for log in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            norm = LogNorm(vmin=1, vmax=counts.max()) if log else None
            mesh = ax.pcolormesh(gen_edges, reco_edges, np.ma.masked_equal(counts.T, 0),
                                 norm=norm, cmap='viridis', shading='flat')
            fig.colorbar(mesh, ax=ax, label='Gen photons')
            ax.plot([gen_edges[0], gen_edges[-1]], [gen_edges[0], gen_edges[-1]],
                    color='grey', linestyle='--', linewidth=1)
            # One population label, using the same empty-handle convention as the statistics.
            ax.legend([Line2D([], [], color='none')], [label], loc='upper right',
                      frameon=False, handlelength=0, handletextpad=0)
            ax.set_xlim(gen_edges[0], gen_edges[-1])
            # Leave a white band above the overflow bin for the sample and population.
            ax.set_ylim(reco_edges[0], reco_edges[-1] + 10)
            ax.set_xticks([MIN_GEN_ENERGY, 10, 20, 30, 40, 50])
            ax.set_yticks([0, 10, 20, 30, 40, 55], ['0', '10', '20', '30', '40', r'$\geq50$'])
            suffix = '_log' if log else ''
            finish(fig, ax, output / f'gen_reco_energy_2d_{population}_{sample}{suffix}.png',
                   sample, XLABELS['energy'], r'$\sum E^{\rm reco}$ [GeV]')


def detector_efficiency(study, samples):
    """Pooled regional efficiencies for both matching methods and all populations."""
    energy = np.concatenate([value['gamma']['energy'] for value in samples.values()])
    cosine = np.concatenate([value['gamma']['cos_theta'] for value in samples.values()])
    theta = np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))
    regions = {'All angles': np.ones(len(theta), dtype=bool),
               'HPC': (theta > 40) & (theta < 140),
               'FEMC': ((theta > 10) & (theta < 37)) | ((theta > 143) & (theta < 170)),
               'STIC': ((theta > 2) & (theta < 10)) | ((theta > 170) & (theta < 178))}
    x = np.arange(len(regions))
    for source, method in (('reco', 'truth_matching'), ('angular', 'angular_matching')):
        output = study / '05_detector_efficiency' / method
        output.mkdir(parents=True, exist_ok=True)
        matched = np.concatenate([value[source]['reco_count'] > 0 for value in samples.values()])
        for index, (population, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
            selected = np.concatenate([value['analysis_selections'][population] for value in samples.values()])
            selected &= energy >= MIN_GEN_ENERGY
            denominator = np.array([(selected & region).sum() for region in regions.values()])
            numerator = np.array([(selected & region & matched).sum() for region in regions.values()])
            efficiency = numerator / denominator
            lower, upper = clopper_pearson(numerator, denominator)
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            ax.bar(x, efficiency, color=f'C{index}', width=0.6, label=label)
            ax.errorbar(x, efficiency, yerr=[efficiency - lower, upper - efficiency],
                        fmt='none', color='black', capsize=4)
            for position, passed, total, high in zip(x, numerator, denominator, upper, strict=True):
                text = f'{100 * passed / total:.1f}%\n' + rf'$\frac{{{passed:,}}}{{{total:,}}}$'
                ax.text(position, high + 0.025, text, ha='center', va='bottom')
            ax.text(0.02, 0.89, rf'$E_\gamma^{{\rm gen}}\geq{MIN_GEN_ENERGY}$ GeV',
                    transform=ax.transAxes, va='top')
            ax.set_xticks(x, regions)
            ax.set_ylim(0, 1.35)
            ax.set_yticks(np.linspace(0, 1, 6))
            finish(fig, ax, output / f'efficiency_by_detector_{population}_combined.png',
                   'combined', 'Region', 'Efficiency')
            print(f'{study.name} {method} {population} E>={MIN_GEN_ENERGY} detector counts: '
                  f'{dict(zip(regions, zip(numerator, denominator), strict=True))}', flush=True)


def matching_results(study, sample, values, bins):
    gamma = values['gamma']
    for source, method in (('reco', 'truth_matching'), ('angular', 'angular_matching')):
        output = study / '04_matching_result' / method
        output.mkdir(parents=True, exist_ok=True)
        stats = values[source]
        matched_energy_distributions(output, sample, values, stats)
        gen_reco_energy_histograms(output, sample, values, stats)
        for coordinate in ('cos_theta', 'energy'):
            edges = bins[coordinate]
            centers = (edges[:-1] + edges[1:]) / 2
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for index, (gen_name, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
                color = f'C{index}'
                selected = values['analysis_selections'][gen_name]
                x = gamma[coordinate][selected]
                denominator = np.histogram(x, bins=edges)[0]
                assert denominator.sum() == selected.sum()
                matched = stats['reco_count'][selected] > 0
                numerator = np.histogram(x[matched], bins=edges)[0]
                assert numerator.sum() == matched.sum() and np.all(numerator <= denominator)
                shown = denominator > 0
                ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=shown)
                lower, upper = clopper_pearson(numerator, denominator)
                ax.stairs(ratio, edges, baseline=None, color=color, label=label)
                ax.errorbar(centers[shown], ratio[shown], yerr=np.maximum([ratio[shown] - lower[shown], upper[shown] - ratio[shown]], 0),
                            fmt='none', color=color, capsize=2)
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0, 1.5)
            ax.set_yticks(np.linspace(0, 1, 6))
            finish(fig, ax, output / f'efficiency_vs_{coordinate}_{sample}.png', sample, XLABELS[coordinate], 'Efficiency')
            profiles = []
            for index, (gen_name, label) in enumerate(zip(ANALYSIS_SELECTIONS, GEN_LABELS, strict=True)):
                selected = values['analysis_selections'][gen_name] & (gamma['energy'] >= MIN_GEN_ENERGY)
                x = gamma[coordinate][selected]
                mean, error = profile_points(x, stats['energy_ratio'][selected], edges)
                profiles.append((mean, error, f'C{index}', label))
            upper = max(np.nanmax(mean + np.nan_to_num(error, nan=0)) for mean, error, _, _ in profiles)
            for log_y in (False, True):
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                for mean, error, color, label in profiles:
                    ax.stairs(mean, edges, baseline=None, color=color, label=label)
                    ax.errorbar(centers, mean, yerr=error, fmt='none', color=color, capsize=2)
                ax.set_xlim(MIN_GEN_ENERGY if coordinate == 'energy' else edges[0], edges[-1])
                ax.text(0.02, 0.89, rf'$E_\gamma^{{\rm gen}}\geq{MIN_GEN_ENERGY}$ GeV',
                        transform=ax.transAxes, va='top')
                if log_y:
                    ax.set_yscale('log')
                    positive = np.concatenate([mean for mean, _, _, _ in profiles])
                    positive = positive[np.isfinite(positive) & (positive > 0)]
                    lower = positive.min() / 3
                    # Reserve the same upper fraction as in the linear plot for the legend.
                    ax.set_ylim(lower, lower * (upper / lower) ** 1.65)
                else:
                    ax.set_ylim(0, upper * 1.65)
                suffix = '_log' if log_y else ''
                finish(fig, ax, output / f'energy_response_vs_{coordinate}_{sample}{suffix}.png', sample,
                       XLABELS[coordinate], r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$')


def plot_isr_study(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    study = output_root / 'isr_study'
    study.mkdir(parents=True, exist_ok=True)
    # Read each event once: Gen selection, truth association, and angular scan.
    angles = {mode: {} for mode in RECO_LABELS}
    values = {mode: {} for mode in RECO_LABELS}
    for sample in SAMPLES:
        for mode in RECO_LABELS:
            angles[mode][sample] = new_accumulator()
        def scan(event):
            for mode, objects in event['reco'].items():
                accumulate(angles[mode][sample], event, objects)
        samples = read_sample(input_root, sample, scan)
        for mode, value in samples.items():
            values[mode][sample] = value
    # Identical bins make the two Reco definitions directly comparable.
    all_energy = np.concatenate([np.r_[value['gamma']['energy'], value['reco_distributions']['energy']]
                                 for samples in values.values() for value in samples.values()])
    bins = dict(energy=np.r_[ENERGY_EDGES, 60], energy_log=energy_bins(all_energy, True),
                cos_theta=np.linspace(-1, 1, 21))
    gen_bins = dict(bins, energy=energy_bins(np.concatenate([value['gamma']['energy'] for value in values['gamma'].values()])))
    for mode in RECO_LABELS:
        output = study / mode
        # Apply the same 99% pooled non-beam ISR truth-retention criterion to each mode.
        cut = choose_cut(angles[mode].values())
        for sample, value in values[mode].items():
            value['angular'] = angular_stats(angles[mode][sample], cut)
            sample_distributions(output, sample, value, bins)
            truth_link_validation(output, sample, value)
            plot_angle_matching(output / '03_angular_matching_validation', sample, angles[mode][sample], cut)
            matching_results(output, sample, value, gen_bins)
            print(f'plots: {mode} {sample}', flush=True)
        detector_efficiency(output, values[mode])
        print(f'ISR study: {output}; opening-angle reference {cut:g} deg', flush=True)
