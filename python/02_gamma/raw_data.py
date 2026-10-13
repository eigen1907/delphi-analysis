"""Raw EM clusters: directions and energies, without a Sim truth association."""
import numpy as np

from .data import gen_event

GEN_P4 = tuple(f'GenPart_fourMomentum.fCoordinates.f{axis}' for axis in 'XYZT')
EM_POSITION = tuple(f'EmShower_position.fCoordinates.f{axis}' for axis in 'XYZ')
STIC_MOMENTUM = tuple(f'Stic_directionFromMain.fCoordinates.f{axis}' for axis in 'XYZ')
RECO_LABELS = dict(EmShower='EM shower', Stic='STIC', EmShower_Stic='EM shower + STIC')
BRANCHES = ('GenPart_pdgId', 'GenPart_status', 'GenPart_parentIdx', *GEN_P4,
            'Event_simPVX', 'Event_simPVY', 'Event_simPVZ', 'Event_simPVErrorFlag',
            'EmShower_paIdx', 'EmShower_energy', *EM_POSITION,
            'Stic_paIdx', 'Stic_energyFromMain', *STIC_MOMENTUM,
            'TracRaw_paIdx', 'TracRaw_charge')


def candidates(energy, vectors, pa_indices, charged_pa):
    """Keep every usable cluster; sharing a PA does not make two showers duplicates."""
    energy, vectors = np.asarray(energy), np.asarray(vectors)
    length = np.linalg.norm(vectors, axis=1)
    valid = np.isfinite(energy) & (energy > 0) & np.isfinite(vectors).all(axis=1) & (length > 0)
    p4 = np.column_stack((vectors[valid] * (energy[valid] / length[valid])[:, None], energy[valid]))
    counts = dict(stored=len(energy), usable=int(valid.sum()), invalid=int((~valid).sum()),
                  charged_PA=int(np.isin(np.asarray(pa_indices)[valid], list(charged_pa)).sum()))
    return dict(p4=p4, geometric_energy_valid=np.ones(len(p4), dtype=bool), branch_counts=counts)


def analyze_event(raw, sample):
    event = gen_event(raw, sample, GEN_P4)
    assert raw['Event_simPVErrorFlag'] == 0, 'EM direction needs a valid simulated primary vertex'
    pv = np.array([raw[f'Event_simPV{axis}'] for axis in 'XYZ'])
    positions = np.column_stack([raw[name] for name in EM_POSITION])
    stic_momentum = np.column_stack([raw[name] for name in STIC_MOMENTUM])
    charged_pa = {pa for pa, charge in zip(raw['TracRaw_paIdx'], raw['TracRaw_charge'], strict=True) if charge != 0}
    em = candidates(raw['EmShower_energy'], positions - pv, raw['EmShower_paIdx'], charged_pa)
    stic = candidates(raw['Stic_energyFromMain'], stic_momentum, raw['Stic_paIdx'], charged_pa)
    # No cross-collection PA overlap was found in the full 100-file sample audit.
    # Multiple EmShower rows per PA are distinct saved clusters and remain separate.
    combined = dict(p4=np.concatenate((em['p4'], stic['p4'])),
                    geometric_energy_valid=np.r_[em['geometric_energy_valid'], stic['geometric_energy_valid']],
                    branch_counts={key: em['branch_counts'][key] + stic['branch_counts'][key]
                                   for key in em['branch_counts']})
    event['reco'] = dict(EmShower=em, Stic=stic, EmShower_Stic=combined)
    return event
