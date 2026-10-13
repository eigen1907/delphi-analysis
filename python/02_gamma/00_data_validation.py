"""Inspect inclusive Reco truth links before any Gen selection."""
import matplotlib.pyplot as plt
import numpy as np

from common.plot_utils import FIGURE_SIZE, finish, plot_path
from .data import LINK_STATUSES

def data_validation(study, sample, values):
    """Inclusive link status includes unresolved candidates, which have no Gen population."""
    native = values['source'] == 'nanoaod'
    labels = (('Agree', 'Forward\nonly', 'Reverse\nonly', 'Conflict', 'No origin') if native
              else ('Stored', 'Usable', 'Invalid', 'Charged\nPA'))
    fields = LINK_STATUSES if native else ('stored', 'usable', 'invalid', 'charged_PA')
    counts = np.array([(values['link_counts'] if native else values['branch_counts'])[key] for key in fields])
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    x = np.arange(len(labels))
    ax.bar(x, counts, color='C0')
    ax.errorbar(x, counts, yerr=np.sqrt(counts), fmt='none', color='C0', capsize=2)
    ax.set_xticks(x, labels)
    ax.set_yscale('log')
    ax.set_ylim(0.5, max(counts.max() * 100, 1))
    output = plot_path(study, '02_gamma', '00_data_validation', sample, values['source'],
                       values['collection']) / 'branch_checks'
    name = 'reco_link_status' if native else 'cluster_validity'
    finish(fig, ax, output / f'{name}_{sample}.png', sample,
           'Link status' if native else 'Cluster status', 'Reco candidates')


def source_comparison(root, sample, native, raw):
    """Compare Gen records only after joining the same job, run and event."""
    first = dict(zip(native['event_keys'], native['gen_signatures'], strict=True))
    second = dict(zip(raw['event_keys'], raw['gen_signatures'], strict=True))
    assert len(first) == native['events'] and len(second) == raw['events'], 'Duplicate event IDs within a job'
    common = first.keys() & second.keys()
    equal = sum(first[key] == second[key] for key in common)
    counts = np.array([len(common), equal, len(common) - equal,
                       len(first.keys() - second.keys()), len(second.keys() - first.keys())])
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    x = np.arange(len(counts))
    ax.bar(x, counts, color='C0')
    ax.set_xticks(x, ('Paired', 'Gen equal', 'Gen different', 'Only NanoAOD', 'Only raw'))
    ax.set_ylim(0, max(counts.max(), 1) * 1.5)
    output = root / '02_gamma' / '00_data_validation' / sample / 'source_comparison'
    finish(fig, ax, output / f'gen_correspondence_{sample}.png', sample, 'Source comparison', 'Events')
    print(f'{sample}: paired={len(common)}, exact Gen agreement={equal}; '
          f'empty Gen records excluded: NanoAOD={native["metadata_excluded"]}, raw={raw["metadata_excluded"]}', flush=True)
    # Reco collections are different saved representations: compare inclusive
    # spectra, rather than assuming row-by-row cluster/object correspondence.
    for coordinate, bins, xlabel in (
            ('energy', np.array([0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 60]), r'$E$ [GeV]'),
            ('cos_theta', np.linspace(-1, 1, 101), r'$\cos\theta$')):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for index, (source, values) in enumerate((('NanoAOD', native), ('raw SDST', raw))):
            entries = values['reco_distributions'][coordinate]
            entries = entries[np.isfinite(entries)]
            if coordinate == 'energy':
                entries = np.minimum(entries, np.nextafter(bins[-1], -np.inf))
            counts = np.histogram(entries, bins)[0]
            assert counts.sum() == len(entries)
            ax.stairs(counts / values['events'], bins, color=f'C{index}', label=f"{source}: {values['reco_label']}")
        ax.set_yscale('log')
        ax.set_ylim(0.5 / min(native['events'], raw['events']), ax.get_ylim()[1] * 100)
        if coordinate == 'energy':
            ax.set_xticks([0, 10, 20, 30, 40, 55], ['0', '10', '20', '30', '40', r'$\geq50$'])
        finish(fig, ax, output / f'reco_{coordinate}_{sample}_log.png', sample, xlabel, 'Objects / bin per event')
