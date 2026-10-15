"""Check saved, selected and usable track counts."""
import numpy as np
import matplotlib.pyplot as plt

from common.plot_utils import FIGURE_SIZE, finish, plot_path
from .results import STAGES


def plot(root, dataset, source, collection, values):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    labels = ("Events", "Gen charged", "Saved tracks", "Usable tracks", "Metadata")
    counts = np.array([values['counts'][key] for key in ("events", "stable_charged", "raw_tracks", "valid_tracks", "metadata_excluded")])
    ax.bar(np.arange(5), counts, color='C0')
    ax.errorbar(np.arange(5), counts, yerr=np.sqrt(counts), fmt='none', color='black', capsize=3)
    ax.set_xticks(np.arange(5), labels)
    ax.set_ylim(0, counts.max(initial=0) * 1.7)
    output = plot_path(root, '01_track', STAGES[0], dataset, source, collection) / 'branch_checks'
    finish(fig, ax, output / 'object_counts.png', dataset, '', 'Count')


def source_comparison(root, dataset, raw, normal):
    """Join unique event IDs, then compare Gen records and independent Reco views."""
    output = root / '01_track' / STAGES[0] / dataset / 'source_comparison'
    raw_ids = set(raw['event_records']) - raw['duplicate_keys']
    normal_ids = set(normal['event_records']) - normal['duplicate_keys']
    shared = raw_ids & normal_ids
    agree = sum(raw['event_records'][key] == normal['event_records'][key] for key in shared)
    counts = np.array([len(shared), agree, len(shared) - agree, len(raw_ids - normal_ids),
                       len(normal_ids - raw_ids), len(raw['duplicate_keys']), len(normal['duplicate_keys'])])
    labels = ('Shared\nIDs', 'Gen\nagree', 'Gen\ndiffer', 'Raw\nonly', 'Normal\nonly',
              'Raw\nambiguous', 'Normal\nambiguous')
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.bar(np.arange(len(counts)), counts, color='C0')
    ax.errorbar(np.arange(len(counts)), counts, yerr=np.sqrt(counts), fmt='none', color='black', capsize=3)
    ax.set_xticks(np.arange(len(counts)), labels)
    ax.tick_params(axis='x', labelsize=16)
    ax.set_yscale('log')
    ax.set_ylim(0.5, counts.max(initial=1) * 100)
    finish(fig, ax, output / 'event_gen_correspondence.png', dataset, '', 'Events')
    # These are different Reco definitions, not arrays expected to be identical.
    from .results import histogram
    for name, bins, xlabel in (('p', np.r_[np.linspace(0, 60, 31), 70], r'$p$ [GeV]'),
                               ('cos_theta', np.linspace(-1, 1, 101), r'$\cos\theta$')):
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        for index, (values, label) in enumerate(((raw, 'TracRaw'), (normal, 'Trac → Part'))):
            entries = values['reco'][name]
            if name == 'p':
                entries = np.minimum(entries, np.nextafter(70, -np.inf))
            histogram(ax, entries, bins, label, f'C{index}', values['counts']['events'])
        if name == 'p':
            ax.set_xticks([0, 20, 40, 65], ['0', '20', '40', r'$\geq60$'])
        ax.set_ylim(0, ax.get_ylim()[1] * 1.7)
        finish(fig, ax, output / f'reco_{name}.png', dataset, xlabel, 'Tracks / bin per event')
    print(f'Track {dataset} source comparison: {dict(zip(labels, counts, strict=True))}', flush=True)
