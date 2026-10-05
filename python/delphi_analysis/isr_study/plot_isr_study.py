"""ISR link validation, photon/conversion reconstruction, and full Reco footprint."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as mh
import numpy as np
from matplotlib.ticker import MaxNLocator
from scipy.stats import beta

from .angle_matching import accumulate, new_accumulator, plot_angle_matching, summary as angle_summary
from .data import CHANNELS, SAMPLES, read_sample

FIGURE_SIZE = (11, 9)
COVERAGE = 0.6826894921370859
POPULATIONS = ('all_ISR', 'noBeamISR')
STAGES = {'photon_conversion': '01_photon_conversion', 'all_lineage': '02_all_lineage'}
XLABELS = {'energy': r'$E^{\rm gen}_\gamma$ [GeV]', 'pt': r'$p_T^{\rm gen}$ [GeV]',
           'cos_theta': r'$\cos\theta^{\rm gen}$'}
COMPOSITION = (('charged_e', r'Charged $e^\pm$'), ('other_charged', 'Other charged'),
               ('neutral_gamma', r'Neutral $\gamma$'), ('other_neutral', 'Other neutral'))


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
    legend = ax.legend(title=sample, loc='best', frameon=True, framealpha=1, edgecolor='none')
    legend.get_title().set_fontweight('bold')
    legend.get_frame().set_facecolor('white')
    ax.grid(alpha=0.2)
    mh.label.exp_label(exp='DELPHI', llabel='Simulation', rlabel='LEP 1 (91.2 GeV)', loc=0, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def log_bins(entries, bins_per_decade=3):
    positive = entries[entries > 0]
    low = positive.min() / 2
    high = np.nextafter(positive.max(), np.inf)
    count = max(1, int(np.ceil(bins_per_decade * np.log10(high / low))))
    return np.geomspace(low, high, count + 1)


def coordinate_bins(entries, coordinate, linear):
    if coordinate == 'cos_theta':
        return np.linspace(-1, 1, 21), 'linear'
    if linear:
        return np.arange(0, np.ceil(entries.max()) + 1), 'linear'
    bins = log_bins(entries)
    return (np.r_[0, bins], 'symlog') if coordinate == 'pt' else (bins, 'log')


def centers_and_scale(ax, bins, scale):
    centers = (bins[:-1] + bins[1:]) / 2
    if scale == 'log':
        centers = np.sqrt(bins[:-1] * bins[1:])
        ax.set_xscale('log')
    elif scale == 'symlog':
        centers = np.r_[0, np.sqrt(bins[1:-1] * bins[2:])]
        ax.set_xscale('symlog', linthresh=bins[1])
    ax.set_xlim(bins[0], bins[-1])
    return centers


def histogram(path, sample, entries, bins, normalization, xlabel, label, linear, scale='linear'):
    counts = np.histogram(entries, bins=bins)[0]
    assert counts.sum() == len(entries), 'Histogram must retain all entries, including tails'
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    centers = centers_and_scale(ax, bins, scale)
    ax.stairs(counts / normalization, bins, label=label, color='C0')
    shown = counts > 0
    ax.errorbar(centers[shown], counts[shown] / normalization,
                yerr=np.sqrt(counts[shown]) / normalization, fmt='o', color='C0', markersize=5, capsize=2)
    if linear:
        ax.set_ylim(0, 1.4 * (counts + np.sqrt(counts)).max() / normalization)
    else:
        ax.set_yscale('log')
        ax.set_ylim(top=10 * counts.max() / normalization)
    if np.issubdtype(entries.dtype, np.integer):
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    finish(fig, ax, path, sample, xlabel, 'Gen fraction / bin')


def efficiency(path, sample, coordinates, matched, bins, scale, label, xlabel):
    denominator = np.histogram(coordinates, bins=bins)[0]
    numerator = np.histogram(coordinates[matched], bins=bins)[0]
    assert denominator.sum() == len(coordinates) and numerator.sum() == matched.sum()
    assert np.all(numerator <= denominator)
    shown = denominator > 0
    ratio = np.divide(numerator, denominator, out=np.full(len(denominator), np.nan), where=shown)
    lower, upper = clopper_pearson(numerator, denominator)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    centers = centers_and_scale(ax, bins, scale)
    ax.stairs(ratio, bins, baseline=None, label=label, color='C0')
    ax.errorbar(centers[shown], ratio[shown], yerr=[ratio[shown] - lower[shown], upper[shown] - ratio[shown]],
                fmt='o', color='C0', markersize=5, capsize=2)
    ax.set_ylim(0, 1.25)
    ax.set_yticks(np.linspace(0, 1, 6))
    finish(fig, ax, path, sample, xlabel, 'Efficiency')


def profile(path, sample, coordinates, values, bins, scale, xlabel, ylabel, label):
    defined = np.isfinite(values)
    counts = np.histogram(coordinates[defined], bins=bins)[0]
    sums = np.histogram(coordinates[defined], bins=bins, weights=values[defined])[0]
    squares = np.histogram(coordinates[defined], bins=bins, weights=values[defined] ** 2)[0]
    assert counts.sum() == defined.sum()
    mean = np.divide(sums, counts, out=np.full(len(counts), np.nan), where=counts > 0)
    measured = counts > 1
    sem = np.sqrt(np.maximum(squares[measured] - sums[measured] ** 2 / counts[measured], 0)
                  / (counts[measured] * (counts[measured] - 1)))
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    centers = centers_and_scale(ax, bins, scale)
    ax.stairs(mean, bins, baseline=None, color='C0', label=label)
    ax.plot(centers[counts > 0], mean[counts > 0], 'o', color='C0', markersize=5)
    ax.errorbar(centers[measured], mean[measured], yerr=sem, fmt='none', color='C0', capsize=2)
    ax.set_ylim(bottom=0)
    finish(fig, ax, path, sample, xlabel, ylabel)


def gen_spectra(study, sample, values, common):
    output = study / '01_photon_conversion' / 'gen' / sample
    output.mkdir(parents=True, exist_ok=True)
    photons = values['photons']
    for coordinate in XLABELS:
        for linear in (False, True):
            bins, scale = coordinate_bins(common[coordinate], coordinate, linear)
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            centers = centers_and_scale(ax, bins, scale)
            bottom = np.zeros(len(bins) - 1)
            for selected, label, color, hatch in ((~photons['beam'], 'Non-beam ISR', 'C0', ''),
                                                  (photons['beam'], 'Beam ISR', 'C1', '//')):
                counts = np.histogram(photons[coordinate][selected], bins=bins)[0]
                ax.stairs(bottom + counts / values['events'], bins, baseline=bottom, fill=True,
                          label=label, facecolor=color, edgecolor='black', hatch=hatch, linewidth=0.7)
                bottom += counts / values['events']
            total = np.histogram(photons[coordinate], bins=bins)[0]
            assert total.sum() == len(photons['energy'])
            shown = total > 0
            ax.errorbar(centers[shown], bottom[shown], yerr=np.sqrt(total[shown]) / values['events'],
                        fmt='o', color='black', markersize=5, capsize=2)
            if linear:
                ax.set_ylim(0, bottom.max() * 1.4)
            else:
                ax.set_yscale('log')
                ax.set_ylim(top=bottom.max() * 10)
            suffix = '_linear' if linear else ''
            finish(fig, ax, output / f'gen_{coordinate}{suffix}.png', sample,
                   XLABELS[coordinate], r'$N_\gamma$ / bin per event')


def link_validation(study, sample, values):
    output = study / '00_link_validation' / sample
    output.mkdir(parents=True, exist_ok=True)
    keys = ('agreement', 'forward_only', 'reverse_only', 'conflict', 'no_origin')
    labels = ('Agree', 'Forward\nonly', 'Reverse\nonly', 'Conflict', 'No origin')
    for name, filename, label in (('all_parts', 'reco_link_status', 'All Part'),
                                  ('ISR_involved', 'reco_link_status_isr', 'ISR-related Part')):
        counts = np.asarray([values['link_validation'][name].get(key, 0) for key in keys])
        for linear in (False, True):
            fig, ax = plt.subplots(figsize=FIGURE_SIZE)
            ax.bar(np.arange(len(keys)), counts, color='C0', label=label)
            ax.errorbar(np.arange(len(keys)), counts, yerr=np.sqrt(counts), fmt='o', color='black', capsize=2)
            ax.set_xticks(np.arange(len(keys)), labels)
            if linear:
                ax.set_ylim(0, counts.max() * 1.4)
            else:
                ax.set_yscale('log')
                ax.set_ylim(top=counts.max() * 10)
            suffix = '_linear' if linear else ''
            finish(fig, ax, output / f'{filename}{suffix}.png', sample, 'Link status', 'Reco Parts')


def composition_plot(path, sample, photons, raw_count, composition, selected, coordinate, bins, scale):
    gen_counts = np.histogram(photons[coordinate][selected], bins=bins)[0]
    records = np.ones(len(composition['category']), dtype=bool) if selected.all() else ~composition['beam']
    record_x = composition['gen_' + coordinate]
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    centers = centers_and_scale(ax, bins, scale)
    bottom = np.zeros(len(bins) - 1)
    total = np.zeros(len(bins) - 1, dtype=int)
    for index, (category, label) in enumerate(COMPOSITION):
        counts = np.histogram(record_x[records & (composition['category'] == category)], bins=bins)[0]
        mean = np.divide(counts, gen_counts, out=np.zeros(len(counts)), where=gen_counts > 0)
        ax.stairs(bottom + mean, bins, baseline=bottom, fill=True, label=label,
                  facecolor=f'C{index}', edgecolor='black', hatch=('', '//', '\\\\', '..')[index], linewidth=0.7)
        bottom += mean
        total += counts
    sums = np.histogram(photons[coordinate][selected], bins=bins, weights=raw_count[selected])[0]
    squares = np.histogram(photons[coordinate][selected], bins=bins, weights=raw_count[selected] ** 2)[0]
    assert np.array_equal(total, sums), 'Composition must retain every associated Part'
    measured = gen_counts > 1
    sem = np.sqrt(np.maximum(squares[measured] - sums[measured] ** 2 / gen_counts[measured], 0)
                  / (gen_counts[measured] * (gen_counts[measured] - 1)))
    ax.errorbar(centers[measured], bottom[measured], yerr=sem, fmt='o', color='black', markersize=5, capsize=2)
    ax.set_ylim(0, max(0.1, bottom.max() * 1.5))
    finish(fig, ax, path, sample, XLABELS[coordinate], r'$\langle N_{\rm Part}\rangle$ / gen')


def channel_plots(study, sample, values, channel, population, common, response_ranges):
    output = study / STAGES[channel] / population / sample
    output.mkdir(parents=True, exist_ok=True)
    photons, stats = values['photons'], values['channels'][channel]
    selected = np.ones(len(photons['energy']), dtype=bool) if population == 'all_ISR' else ~photons['beam']
    n_gen = int(selected.sum())
    label = 'Photon + conversion' if channel == 'photon_conversion' else 'Any Part'
    for coordinate in XLABELS:
        for linear in (False, True):
            if coordinate == 'cos_theta' and linear:
                continue  # These efficiency and profile plots already have fully linear axes.
            bins, scale = coordinate_bins(common[coordinate], coordinate, linear)
            suffix = '_linear' if linear else ''
            efficiency(output / f'efficiency_vs_{coordinate}{suffix}.png', sample,
                       photons[coordinate][selected], stats['reco_count'][selected] > 0, bins, scale, label, XLABELS[coordinate])
            profile(output / f'energy_response_vs_{coordinate}{suffix}.png', sample,
                    photons[coordinate][selected], stats['energy_ratio'][selected], bins, scale,
                    XLABELS[coordinate], r'$\langle E^{\rm reco}/E^{\rm gen}\rangle$', 'Defined response')
            if channel == 'all_lineage' and coordinate in ('pt', 'cos_theta'):
                composition_plot(output / f'composition_vs_{coordinate}{suffix}.png', sample,
                                 photons, stats['raw_count'], values['composition'], selected, coordinate, bins, scale)
    for linear in (False, True):
        suffix = '_linear' if linear else ''
        counts = stats['reco_count'][selected]
        bins = np.arange(-0.5, response_ranges['multiplicity'] + 1.5)
        histogram(output / f'reco_multiplicity{suffix}.png', sample, counts, bins, n_gen,
                  r'$N^{\rm reco}$ / gen', label, linear)
        fields = (('energy_ratio', 'energy_response', 'Defined response'),
                  ('leading_ratio', 'leading_energy_response', 'Matched, defined response'))
        if channel == 'all_lineage':
            fields += (('raw_energy_ratio', 'raw_activity_response', 'Raw Part activity'),)
        for field, filename, legend in fields:
            entries = stats[field][selected]
            entries = entries[np.isfinite(entries)]
            all_entries = response_ranges[field]
            bins = np.linspace(0, np.nextafter(all_entries.max(), np.inf), 81) if linear else np.r_[0, log_bins(all_entries, 6)]
            scale = 'linear' if linear else 'symlog'
            normalization = len(entries) if field == 'leading_ratio' else n_gen
            histogram(output / f'{filename}{suffix}.png', sample, entries, bins, normalization,
                      r'$E^{\rm reco}/E^{\rm gen}$', legend, linear, scale)


def response_summary(stats, selected):
    matched = stats['reco_count'][selected] > 0
    ratio = stats['energy_ratio'][selected]
    defined = np.isfinite(ratio)
    matched_defined = ratio[defined & matched]
    return dict(defined_gen=int(defined.sum()), undefined_gen=int((~defined).sum()),
                mean_defined_gen=float(ratio[defined].mean()),
                matched_quantiles_16_50_84=np.quantile(matched_defined, [0.16, 0.5, 0.84]).tolist(),
                mean_matched_defined=float(matched_defined.mean()),
                maximum=float(ratio[defined].max()), above_one=int(np.count_nonzero(ratio > 1)),
                raw_mean_all_gen=float(stats['raw_energy_ratio'][selected].mean()))


def write_summary(study, values, angles):
    results = {}
    for sample, value in values.items():
        photons = value['photons']
        channels = {}
        for channel in CHANNELS:
            stats = value['channels'][channel]
            channels[channel] = {}
            for population in POPULATIONS:
                selected = np.ones(len(photons['energy']), dtype=bool) if population == 'all_ISR' else ~photons['beam']
                denominator = int(selected.sum())
                numerator = int(np.count_nonzero(stats['reco_count'][selected]))
                lower, upper = clopper_pearson(np.array([numerator]), np.array([denominator]))
                valid = selected & stats['energy_valid']
                energy = photons['energy']
                channels[channel][population] = dict(
                    efficiency=dict(numerator=numerator, denominator=denominator, value=numerator / denominator,
                                    lower_68_percent_cp=float(lower[0]), upper_68_percent_cp=float(upper[0])),
                    raw_objects=int(stats['raw_count'][selected].sum()),
                    canonical_objects=int(stats['energy_count'][selected].sum()),
                    gen_with_multiple_reco=int(np.count_nonzero(stats['reco_count'][selected] > 1)),
                    response=response_summary(stats, selected),
                    ratio_of_defined_total_energies=float(np.sum(energy[valid] * stats['energy_ratio'][valid]) / energy[valid].sum()))
        results[sample] = dict(events=value['events'], gen_ISR=len(photons['energy']),
            gen_Beam_ISR=int(photons['beam'].sum()), link_validation=value['link_validation'],
            diagnostics=value['diagnostics'], channels=channels, angle_matching=angle_summary(angles[sample]))
    result = dict(events=sum(value['events'] for value in values.values()),
                  gen_ISR=sum(len(value['photons']['energy']) for value in values.values()), samples=results)
    (study / 'study_summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    notes = [
        'ISR study: five PHOTOS-FSR samples; one stable Gen ISR photon per analysis unit.',
        'Nominal truth = unique Gen origin in the union of both saved association directions.',
        'Conflicts include known non-ISR/non-photon Gen anchors; no angular fallback.',
        'PhotonConv direct index is not used; parent/daughter saved structure and origins are checked.',
        'Raw Part counts/composition include composite parents and daughters.',
        'Energy response removes saved duplicate representations; contradictory structural groups are undefined.',
        'Unmatched defined response = 0. Undefined response is omitted from profiles, not changed to zero.',
        'Response profiles use defined Gen units, with SEM; singleton errors are undefined.',
        'Count errors sqrt(N); efficiencies use 68.27% Clopper-Pearson intervals.',
        'The pT logarithmic first bin contains zero only; Beam ISR remains visible.',
        'No energy, angle, PID, lock, or Sim-terminal cut; all histogram tails are retained.',
        'Part/SimPart species codes are DELPHI mass hypotheses, not Gen PDG IDs.',
        'Missing truth is unresolved, not proof of failed reconstruction or a fake candidate.',
        'Energy response is associated activity, not a decomposition of ISR energy contributions.',
        'Angle scans are per-Gen cones: one Reco object may enter cones of multiple Gen photons.',
        'Unresolved angular pairs are separate from known non-ISR/other-ISR contamination.',
        'Geometric energy profiles exclude cones containing ambiguous energy representations; coverage is saved.',
        'An unresolved retained parent overlapping retained daughters is invalid only for geometric energy.',
        'Angle truth recovery is conditional on nominal truth success; any-candidate probability uses all Gen units.',
        'No best opening-angle cut is selected. Previous outputs and manual audit remain preserved.',
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
    common = {name: np.concatenate([value['photons'][name] for value in values.values()]) for name in XLABELS}
    for sample, value in values.items():
        link_validation(study, sample, value)
        gen_spectra(study, sample, value, common)
    for channel in CHANNELS:
        ranges = {field: np.concatenate([value['channels'][channel][field][np.isfinite(value['channels'][channel][field])]
                                         for value in values.values()])
                  for field in ('energy_ratio', 'leading_ratio', 'raw_energy_ratio')}
        ranges['multiplicity'] = max(value['channels'][channel]['reco_count'].max() for value in values.values())
        for population in POPULATIONS:
            for sample, value in values.items():
                channel_plots(study, sample, value, channel, population, common, ranges)
                print(f'plots: {STAGES[channel]} {population} {sample}', flush=True)
    for sample in SAMPLES:
        plot_angle_matching(study / '03_angle_matching', sample, angles[sample])
    write_summary(study, values, angles)
    print(f'ISR study: {study}', flush=True)
