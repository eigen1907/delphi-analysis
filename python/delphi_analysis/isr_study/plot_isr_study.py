"""Sample spectra, stored truth links, angular validation, and matching results."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from matplotlib.ticker import MaxNLocator
from scipy.stats import beta

from .angle_matching import accumulate, choose_cut, new_accumulator, plot_angle_matching, summary as angle_summary
from .data import GEN_SELECTIONS, LINK_STATUSES, RECO_DEFINITIONS, SAMPLES, read_sample

FIGURE_SIZE = (11, 9)
COVERAGE = 0.6826894921370859
ENERGY_EDGES = np.array([0, 0.1, 0.2, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50])
GEN_COMPONENTS = (('beam_isr', 'BeamISR'), ('nonbeam_isr', 'NonBeamISR'),
                  ('fsr', 'FSR'), ('decayed', 'Decayed'))
RECO_LABELS = ('Photon', 'Photon + conversion', 'All Part')
RESULT_DIRS = dict(gen_gamma_all='gen_all_gamma', gen_isr_all='gen_all_isr',
                   gen_isr_non_beam='gen_non_beam_isr')
XLABELS = dict(energy=r'$E_\gamma^{\rm gen}$ [GeV]', cos_theta=r'$\cos\theta^{\rm gen}$')


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
    legend = ax.get_legend()
    if legend is None:
        legend = ax.legend(title=sample, loc='best', frameon=True, framealpha=1, edgecolor='none')
    legend.get_title().set_fontweight('bold')
    legend.get_frame().set_facecolor('white')
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
    if logarithmic:
        # The first interval contains zero only, rather than dropping zero-energy candidates.
        ax.set_xscale('symlog', linthresh=bins[1])
        centers = np.r_[0, np.sqrt(bins[1:-1] * bins[2:])]
    else:
        centers = (bins[:-1] + bins[1:]) / 2
    ax.set_xlim(bins[0], bins[-1])
    if logarithmic:
        # Space the labels around zero; the bins and all entries are unchanged.
        ax.set_xticks(np.r_[0, ax.get_xticks()[2:]])
    return centers


def stacked_spectrum(ax, entries, category, components, bins, events, logarithmic):
    centers = bin_centers(ax, bins, logarithmic)
    bottom = np.zeros(len(bins) - 1)
    for key, label in components:
        index = {'beam_isr': 0, 'nonbeam_isr': 1, 'fsr': 2, 'decayed': 3,
                 'gamma': 0, 'gamma_conv': 1}[key]
        selected = entries[(category == key) & np.isfinite(entries)]
        counts = np.histogram(selected, bins=bins)[0]
        assert counts.sum() == len(selected), 'Spectrum must retain all finite entries'
        ax.stairs(bottom + counts / events, bins, baseline=bottom, fill=True, label=label,
                  facecolor=f'C{index}', edgecolor='black', hatch=('', '//', '\\\\', '..')[index], linewidth=0.7)
        bottom += counts / events
    total = np.histogram(entries[np.isfinite(entries)], bins=bins)[0]
    assert np.allclose(bottom * events, total)
    shown = total > 0
    ax.errorbar(centers[shown], bottom[shown], yerr=np.sqrt(total[shown]) / events,
                fmt='o', color='black', capsize=2, markersize=5)
    if logarithmic:
        ax.set_yscale('log')
        ax.set_ylim(top=max(bottom.max() * 100, 1 / events))
    else:
        ax.set_ylim(0, max(bottom.max() * 1.8, 1 / events))


def sample_distributions(study, sample, values, bins):
    gamma = values['gamma']
    candidates = values['reco_distributions']['reco_gamma_plus_conv']
    for gen_name in GEN_SELECTIONS:
        output = study / '01_sample_distribution' / gen_name
        output.mkdir(parents=True, exist_ok=True)
        gen_selected = values['selections'][gen_name]
        reco_selected = (np.ones(len(candidates['energy']), dtype=bool) if gen_name == 'gen_gamma_all' else
                         np.isin(candidates['origin_category'], ('beam_isr', 'nonbeam_isr')) if gen_name == 'gen_isr_all' else
                         candidates['origin_category'] == 'nonbeam_isr')
        components = dict(gen_gamma_all=GEN_COMPONENTS, gen_isr_all=GEN_COMPONENTS[:2],
                          gen_isr_non_beam=GEN_COMPONENTS[1:2])[gen_name]
        for observable in ('cos_theta', 'energy', 'energy_log'):
            coordinate = 'energy' if observable.startswith('energy') else 'cos_theta'
            logarithmic = observable == 'energy_log'
            gen_entries = gamma[coordinate][gen_selected]
            reco_entries = candidates[coordinate][reco_selected]
            if observable == 'energy':
                # Display all E >= 50 GeV in the labelled final bin; physics arrays stay unchanged.
                last = np.nextafter(bins['energy'][-1], -np.inf)
                gen_entries = np.minimum(gen_entries, last)
                reco_entries = np.minimum(reco_entries, last)
            fig, axes = plt.subplots(2, 1, figsize=(11, 12), sharex=True)
            stacked_spectrum(axes[0], gen_entries, gamma['category'][gen_selected],
                             components, bins[observable], values['events'], logarithmic)
            stacked_spectrum(axes[1], reco_entries, candidates['kind'][reco_selected],
                             (('gamma', r'$\gamma$'), ('gamma_conv', r'$\gamma_{\rm conv}$')),
                             bins[observable], values['events'], logarithmic)
            for ax, title in zip(axes, ('Gen', 'Reco'), strict=True):
                legend = ax.legend(title=f'{sample} · {title}', loc='upper center', ncol=2,
                                   frameon=True, framealpha=1, edgecolor='none')
                legend.get_title().set_fontweight('bold')
                ax.set_ylabel(r'$N_\gamma$ / bin per event')
                ax.grid(alpha=0.2)
            if observable == 'energy':
                axes[1].set_xticks([0, 10, 20, 30, 40, 55], ['0', '10', '20', '30', '40', r'$\geq 50$'])
            axes[1].set_xlabel(r'$E_\gamma$ [GeV]' if coordinate == 'energy' else r'$\cos\theta$')
            finish(fig, axes[0], output / f'{observable}_{sample}.png', sample, '',
                   r'$N_\gamma$ / bin per event')
    output = study / '01_sample_distribution' / 'multiplicity'
    output.mkdir(parents=True, exist_ok=True)
    for gen_name in GEN_SELECTIONS:
        counts = values['reco']['reco_all']['reco_count'][values['selections'][gen_name]]
        multiplicity_plot(output / f'{gen_name}_{sample}.png', sample, [(counts, 'All Part')],
                          r'$N_{\rm matched\ Part}$ / gen', 'Gen fraction / bin')
    for reco_name in RECO_DEFINITIONS[:2]:
        multiplicity_plot(output / f'{reco_name}_{sample}.png', sample,
                          [(values['reco_event_counts'][reco_name], RECO_LABELS[RECO_DEFINITIONS.index(reco_name)])],
                          r'$N^{\rm reco}$ / event', 'Event fraction / bin')


def multiplicity_plot(path, sample, entries, xlabel, ylabel, log_y=False):
    bins = np.arange(-0.5, max(values.max() for values, _ in entries) + 1.5)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    centers = (bins[:-1] + bins[1:]) / 2
    for index, (values, label) in enumerate(entries):
        counts = np.histogram(values, bins=bins)[0]
        assert counts.sum() == len(values)
        shown = counts > 0
        ax.stairs(counts / len(values), bins, label=label, color=f'C{index}')
        ax.errorbar(centers[shown], counts[shown] / len(values), yerr=np.sqrt(counts[shown]) / len(values),
                    fmt='o', color=f'C{index}', capsize=2, markersize=5)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlim(bins[0], bins[-1])
    if log_y:
        ax.set_yscale('log')
    else:
        ax.set_ylim(0, 1.4)
    finish(fig, ax, path, sample, xlabel, ylabel)


def truth_link_validation(study, sample, values):
    labels = ('Agree', 'Forward\nonly', 'Reverse\nonly', 'Conflict', 'No origin')
    for gen_name in GEN_SELECTIONS:
        for reco_name in RECO_DEFINITIONS:
            output = study / '02_truth_link_matching_validation' / gen_name / reco_name
            output.mkdir(parents=True, exist_ok=True)
            counts = np.array([values['matching_validation'][gen_name][reco_name][key] for key in LINK_STATUSES])
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            x = np.arange(len(labels))
            ax.bar(x, counts, color='C0', label=RECO_LABELS[RECO_DEFINITIONS.index(reco_name)])
            ax.errorbar(x, counts, yerr=np.sqrt(counts), fmt='o', color='black', capsize=2)
            ax.set_xticks(x, labels)
            ax.set_yscale('log')
            ax.set_ylim(top=max(counts.max() * 10, 1))
            finish(fig, ax, output / f'truth_link_status_{sample}.png', sample, 'Link status', 'Reco objects')


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


def matching_results(study, sample, values, bins):
    gamma = values['gamma']
    for gen_name in GEN_SELECTIONS:
        output = study / '04_matching_result' / RESULT_DIRS[gen_name]
        output.mkdir(parents=True, exist_ok=True)
        selected = values['selections'][gen_name]
        for coordinate in ('cos_theta', 'energy'):
            edges = bins[coordinate]
            centers = (edges[:-1] + edges[1:]) / 2
            x = gamma[coordinate][selected]
            denominator = np.histogram(x, bins=edges)[0]
            assert denominator.sum() == selected.sum()
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            for index, reco_name in enumerate(RECO_DEFINITIONS):
                matched = values['reco'][reco_name]['reco_count'][selected] > 0
                numerator = np.histogram(x[matched], bins=edges)[0]
                assert numerator.sum() == matched.sum() and np.all(numerator <= denominator)
                shown = denominator > 0
                ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=shown)
                lower, upper = clopper_pearson(numerator, denominator)
                ax.stairs(ratio, edges, baseline=None, color=f'C{index}', label=RECO_LABELS[index])
                ax.errorbar(centers[shown], ratio[shown], yerr=np.maximum([ratio[shown] - lower[shown], upper[shown] - ratio[shown]], 0),
                            fmt='o', color=f'C{index}', capsize=2, markersize=5)
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0, 1.3)
            ax.set_yticks(np.linspace(0, 1, 6))
            finish(fig, ax, output / f'efficiency_vs_{coordinate}_{sample}.png', sample, XLABELS[coordinate], 'Efficiency')
            for log_y in (False, True):
                fig, ax = plt.subplots(figsize=FIGURE_SIZE)
                means, upper_errors = [], []
                for index, reco_name in enumerate(RECO_DEFINITIONS):
                    mean, error = profile_points(x, values['reco'][reco_name]['energy_ratio'][selected], edges)
                    means.append(mean)
                    upper_errors.append(mean + np.nan_to_num(error, nan=0))
                    ax.stairs(mean, edges, baseline=None, color=f'C{index}', label=RECO_LABELS[index])
                    ax.plot(centers, mean, 'o', color=f'C{index}', markersize=5)
                    ax.errorbar(centers, mean, yerr=error, fmt='none', color=f'C{index}', capsize=2)
                ax.set_xlim(edges[0], edges[-1])
                if log_y:
                    ax.set_yscale('log')
                    positive = np.concatenate(means)
                    positive = positive[np.isfinite(positive) & (positive > 0)]
                    ax.set_ylim(positive.min() / 3, np.nanmax(upper_errors) * 3)
                else:
                    ax.set_ylim(bottom=0)
                suffix = '_logy' if log_y else ''
                finish(fig, ax, output / f'energy_response_vs_{coordinate}{suffix}_{sample}.png', sample,
                       XLABELS[coordinate], r'$\langle\sum E^{\rm reco}/E_\gamma^{\rm gen}\rangle$')
        entries = [(values['reco'][name]['reco_count'][selected], label)
                   for name, label in zip(RECO_DEFINITIONS, RECO_LABELS, strict=True)]
        for log_y in (False, True):
            suffix = '_logy' if log_y else ''
            multiplicity_plot(output / f'matched_reco_multiplicity{suffix}_{sample}.png', sample, entries,
                              r'$N^{\rm reco}$ / gen', 'Gen fraction / bin', log_y)


def matching_summary(gamma, stats, selected):
    count = stats['reco_count'][selected]
    ratio = stats['energy_ratio'][selected]
    energy = gamma['energy'][selected]
    defined = np.isfinite(ratio)
    numerator, denominator = int(np.count_nonzero(count)), len(count)
    lower, upper = clopper_pearson(np.array([numerator]), np.array([denominator]))
    matched_response = ratio[defined & (count > 0)]
    quantiles = np.quantile(matched_response, [0.16, 0.5, 0.84]).tolist() if len(matched_response) else None
    return dict(efficiency=dict(numerator=numerator, denominator=denominator, value=numerator / denominator,
                               lower_68_percent_cp=float(lower[0]), upper_68_percent_cp=float(upper[0])),
                raw_objects=int(stats['raw_count'][selected].sum()),
                canonical_objects=int(stats['energy_count'][selected].sum()),
                gen_with_multiple_reco=int(np.count_nonzero(count > 1)),
                response=dict(defined_gen=int(defined.sum()), undefined_gen=int((~defined).sum()),
                              mean_defined_gen=float(ratio[defined].mean()),
                              matched_quantiles_16_50_84=quantiles,
                              above_one=int(np.count_nonzero(ratio > 1)),
                              maximum=float(ratio[defined].max()),
                              ratio_of_defined_total_energies=float(np.sum(energy[defined] * ratio[defined]) / energy[defined].sum())))


def write_summary(study, values, angles, cut, bins):
    samples = {}
    for sample, value in values.items():
        matching, validation = {}, {}
        for gen_name in GEN_SELECTIONS:
            selected = value['selections'][gen_name]
            matching[gen_name] = {name: matching_summary(value['gamma'], value['reco'][name], selected)
                                  for name in RECO_DEFINITIONS}
            validation[gen_name] = {}
            for reco_name in RECO_DEFINITIONS:
                counts = value['matching_validation'][gen_name][reco_name]
                raw = value['unmatched_validation'][reco_name]
                outside = sum(raw.values()) - raw['no_origin'] - sum(counts.values())
                assert outside >= 0 and counts['no_origin'] == 0
                validation[gen_name][reco_name] = dict(counts=counts, known_outside_selected=outside)
        samples[sample] = dict(events=value['events'],
            gen_categories={key: int(np.count_nonzero(value['gamma']['category'] == key)) for key, _ in GEN_COMPONENTS},
            matching=matching, truth_link_validation=validation,
            all_reco_truth_status=value['unmatched_validation'],
            reco_counts={name: dict(candidates=int(counts.sum()), mean_per_event=float(counts.mean()),
                                    zero_momentum=int(np.count_nonzero(~np.isfinite(value['reco_distributions'][name]['cos_theta']))),
                                    above_50_GeV=int(np.count_nonzero(value['reco_distributions'][name]['energy'] >= 50)),
                                    maximum_energy_GeV=float(value['reco_distributions'][name]['energy'].max()))
                         for name, counts in value['reco_event_counts'].items()},
            diagnostics=value['diagnostics'], saved_link_validation=value['link_validation'],
            angular_matching=angle_summary(angles[sample]))
    result = dict(events=sum(value['events'] for value in values.values()),
                  gen_gamma=sum(len(value['gamma']['energy']) for value in values.values()),
                  energy_bins=bins['energy'].tolist(), energy_log_bins=bins['energy_log'].tolist(),
                  linear_spectrum_overflow_GeV=50,
                  gen_energy_bins=energy_bins(np.concatenate([v['gamma']['energy'] for v in values.values()])).tolist(),
                  chosen_opening_angle=cut, samples=samples)
    (study / 'study_summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    notes = [
        'Five PHOTOS-FSR samples; stable Gen gamma = PDG 22, status 1. No energy or fiducial cut.',
        'Decayed is the previous Others category; current photons originate from pi0, eta, or omega decays.',
        'Gen and Reco spectra are separate stacks normalized by event count, not bin width.',
        'All-gamma Reco spectra include every deduplicated Photon/conversion candidate, including unknown/non-photon origins.',
        'ISR Reco spectra require a unique selected Gen origin. They are truth-associated subsets.',
        'The three Gen multiplicities count raw associated Parts per Gen, including zero; Reco multiplicities count candidates per event.',
        'Nominal association is the conflict-free union of both stored directions, stopping at the first Gen anchor.',
        'Truth validation counts Reco objects with any known selected origin; conflicts are retained, No origin is zero by membership.',
        'Unassigned Reco counts are separate in all_reco_truth_status; known outside-selected counts are saved for each case.',
        'Photon/conversion and all-Part energy remove saved parent/daughter duplication; raw Part multiplicity remains separate.',
        'No association = zero response. Ambiguous energy = NaN, omitted from means with coverage saved.',
        'Response is associated reconstructed energy, not a measured fractional contribution from one Gen photon.',
        'Stage03 counts own truth-associated candidates inside each cone; unrelated/unresolved pairs remain JSON diagnostics.',
        'The dashed angle is chosen by pooled non-beam ISR truth retention; it is not applied to nominal stage04 matching.',
        'Zero-momentum Reco contributes energy and multiplicity, but not cos(theta) or angular cones; counts are saved.',
        'energy_log uses a logarithmic positive region plus a zero-only first interval; no zero-energy entry is dropped.',
        'Linear spectra put E >= 50 GeV into a labelled overflow bin. energy_log and matching retain original energies.',
        'Log-y response and multiplicity figures cannot display zeros; their linear counterparts retain zero values.',
        'Log-y profile limits follow positive means; error bars extending below zero reach the plot boundary.',
        'Count errors: sqrt(N). Efficiency errors: 68.27% Clopper-Pearson. Profile errors: SEM when estimable.',
        'All finite histogram entries and energy tails are retained. High-energy linear bins are deliberately wider.',
    ]
    (study / 'study_notes.txt').write_text('\n'.join(notes) + '\n')


def plot_isr_study(input_root: Path, output_root: Path):
    mh.style.use(mh.styles.CMS)
    study = output_root / 'isr_study'
    study.mkdir(parents=True, exist_ok=True)
    angles, values = {}, {}
    for sample in SAMPLES:
        angles[sample] = new_accumulator()
        values[sample] = read_sample(input_root, sample, lambda event: accumulate(angles[sample], event))
    all_energy = np.concatenate([np.r_[value['gamma']['energy'], value['reco_distributions']['reco_gamma_plus_conv']['energy']]
                                 for value in values.values()])
    bins = dict(energy=np.r_[ENERGY_EDGES, 60], energy_log=energy_bins(all_energy, True),
                cos_theta=np.linspace(-1, 1, 21))
    gen_bins = dict(bins, energy=energy_bins(np.concatenate([value['gamma']['energy'] for value in values.values()])))
    cut = choose_cut(angles.values())
    for sample, value in values.items():
        sample_distributions(study, sample, value, bins)
        truth_link_validation(study, sample, value)
        plot_angle_matching(study / '03_angular_matching_validation', sample, angles[sample], cut['angle_degrees'])
        matching_results(study, sample, value, gen_bins)
        print(f'plots: {sample}', flush=True)
    write_summary(study, values, angles, cut, bins)
    print(f'ISR study: {study}; opening-angle reference {cut["angle_degrees"]:g} deg', flush=True)
