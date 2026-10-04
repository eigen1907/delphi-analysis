#!/usr/bin/env python3
"""Dump ten deliberately selected events for a manual photon-truth audit."""
from pathlib import Path

import numpy as np
import uproot

from delphi_analysis.photons.data import (
    BRANCHES, GEN_P4, RECO_P4, SAMPLES, SIM_GAMMA_CODE, TARGET_PDG,
    first_gen_ancestor, photon_origin,
)

PROJECT = Path(__file__).resolve().parents[2]
INPUT = PROJECT / 'data/20260828_florian'
OUTPUT = PROJECT / 'plots/20260828_florian/photon_study/manual_audit'
SIM_P4 = tuple(f'SimPart_fourMomentum.fCoordinates.f{axis}' for axis in 'XYZT')
AUDIT_BRANCHES = (*BRANCHES, *SIM_P4, 'Event_runNumber', 'Event_evtNumber',
                  'GenPart_simIdx', 'SimPart_decayVtxIdx', 'SimPart_charge', 'SimPart_partIdx')
CASES = ('direct ISR reco', 'descendant ISR reco', 'nonprimary photon anchor',
         'split ISR', 'ISR shower with secondary gamma')


def event_info(raw, sample):
    gen_p4 = np.column_stack([raw[name] for name in GEN_P4])
    sim_p4 = np.column_stack([raw[name] for name in SIM_P4])
    reco_p4 = np.column_stack([raw[name] for name in RECO_P4])
    pdg, status, gen_parent = (np.asarray(raw[name]) for name in
                               ('GenPart_pdgId', 'GenPart_status', 'GenPart_parentIdx'))
    daughters = np.flatnonzero((np.abs(pdg) == TARGET_PDG[sample]) & (status == (21 if sample == 'Ztautau' else 1)))
    assert len(daughters) == 2 and gen_parent[daughters[0]] == gen_parent[daughters[1]]
    hard_parent = gen_parent[daughters[0]]
    sim_gen = np.asarray(raw['SimPart_genIdx'])
    sim_origin, sim_decay = raw['SimPart_originVtxIdx'], raw['SimPart_decayVtxIdx']
    incoming = raw['SimVtx_incomingIdx']
    parents = [incoming[v] if v >= 0 else -1 for v in sim_origin]
    children = [[] for _ in sim_gen]
    for s, parent in enumerate(parents):
        if parent >= 0:
            children[parent].append(s)
            assert sim_decay[parent] == sim_origin[s], 'Nonreciprocal Sim vertex edge'
    origins, depths = np.asarray([
        first_gen_ancestor(s, sim_gen, sim_origin, incoming) for s in range(len(sim_gen))
    ], dtype=int).reshape(-1, 2).T
    assert np.all(origins < len(pdg)), 'Sim origin exceeds Gen collection'
    reco_sims = [raw['Part_simIdx'][p] if p >= 0 else -1 for p in raw['Photon_partIdx']]
    reco_at_sim = [[] for _ in sim_gen]
    for r, sim in enumerate(reco_sims):
        if sim >= 0:
            reco_at_sim[sim].append(r)
    photons = {}
    for g in np.flatnonzero((pdg == 22) & (status == 1)).tolist():
        origin = photon_origin(g, pdg, gen_parent, hard_parent, TARGET_PDG[sample])
        p = gen_parent[g]
        beam = (origin == 'isr' and gen_p4[g, 0] == gen_p4[g, 1] == 0
                and p >= 0 and abs(pdg[p]) == 11 and gen_parent[p] < 0)
        anchors = np.flatnonzero(sim_gen == g).tolist()
        # Independent forward traversal; another Gen anchor starts a different lineage.
        lineage, pending = set(), anchors.copy()
        while pending:
            s = pending.pop()
            if s in lineage:
                continue
            lineage.add(s)
            pending.extend(c for c in children[s] if sim_gen[c] < 0 or sim_gen[c] == g)
        assert lineage == set(np.flatnonzero(origins == g)), 'Forward/reverse origin disagreement'
        forward = raw['GenPart_simIdx'][g]
        if forward >= 0:
            assert sim_gen[forward] == g and forward in anchors
        secondary = sorted(s for s in lineage if depths[s] > 0)
        gamma = [s for s in secondary if raw['SimPart_pdgId'][s] == SIM_GAMMA_CODE]
        recos = sorted(r for s in lineage for r in reco_at_sim[s])
        assert len(recos) == len(set(recos))
        leading = max(recos, key=lambda r: reco_p4[r, 3]) if recos else None
        photons[g] = dict(origin='Beam ISR' if beam else 'Non-beam ISR' if origin == 'isr' else 'FSR' if origin == 'fsr' else 'Others',
                          topology='no_sim' if not lineage else 'shower' if secondary else 'direct',
                          anchors=anchors, lineage=lineage, secondary=secondary, gamma=gamma,
                          recos=recos, leading=leading)
    return dict(raw=raw, gen_p4=gen_p4, sim_p4=sim_p4, reco_p4=reco_p4, photons=photons,
                parents=parents, children=children, origins=origins, depths=depths,
                reco_sims=reco_sims, reco_at_sim=reco_at_sim)


def focus_photons(info, reason):
    photons = info['photons']
    if reason == 'direct ISR reco':
        return [g for g, p in photons.items() if p['origin'] == 'Non-beam ISR' and p['topology'] == 'direct' and p['recos']]
    if reason == 'descendant ISR reco':
        return [g for g, p in photons.items() if p['origin'] == 'Non-beam ISR'
                and any(info['depths'][info['reco_sims'][r]] > 0 for r in p['recos'])]
    if reason == 'nonprimary photon anchor':
        return [g for g, p in photons.items() if any(info['parents'][s] >= 0 for s in p['anchors'])]
    if reason == 'split ISR':
        return [g for g, p in photons.items() if p['origin'] == 'Non-beam ISR' and len(p['recos']) > 1]
    return [g for g, p in photons.items() if p['origin'] == 'Non-beam ISR' and p['gamma']]


def p4_text(p4):
    return '(' + ', '.join(f'{name}={value:.9g}' for name, value in zip(('px', 'py', 'pz', 'E'), p4, strict=True)) + ') GeV'


def dump_event(info, sample, path, entry, reason, focus, filename):
    raw, photons = info['raw'], info['photons']
    lines = [f'Sample: {sample} photosFSR', f'ROOT: {path.relative_to(PROJECT)}',
             f'Entry: {entry} (zero-based)', f'Event_runNumber: {raw["Event_runNumber"]}',
             f'Event_evtNumber: {raw["Event_evtNumber"]}', f'Selection: {reason}',
             f'Focus Gen indices: {focus}', '', 'GEN RECORD (all particles; indices are event-local)']
    for g, p4 in enumerate(info['gen_p4']):
        lines.append(f'Gen[{g}] pdg={raw["GenPart_pdgId"][g]} status={raw["GenPart_status"][g]} '
                     f'parent={raw["GenPart_parentIdx"][g]} simIdx={raw["GenPart_simIdx"][g]} {p4_text(p4)}')
    lines += ['', 'STABLE GEN PHOTONS (same definition/units as the photon study)']
    for g, p in photons.items():
        chain, ancestor = [], g
        while ancestor >= 0:
            assert ancestor not in chain, 'Cycle in Gen parents'
            chain.append(ancestor)
            ancestor = raw['GenPart_parentIdx'][ancestor]
        ancestry = ' <- '.join(f'Gen[{i}](pdg={raw["GenPart_pdgId"][i]},status={raw["GenPart_status"][i]})' for i in chain)
        energy = info['gen_p4'][g, 3]
        cosine = info['gen_p4'][g, 2] / np.linalg.norm(info['gen_p4'][g, :3])
        lines += [f'Gen[{g}] {p["origin"]}; topology={p["topology"]}; E={energy:.9g} GeV; cosTheta={cosine:.9g}',
                  f'  Gen ancestry: {ancestry}',
                  f'  anchors={p["anchors"]}; Sim lineage={sorted(p["lineage"])}',
                  f'  secondary Sim={p["secondary"]}; secondary Sim gamma={p["gamma"]}',
                  f'  associated RecoPhoton={p["recos"]}; efficiency_success={int(bool(p["recos"]))} (per Gen photon)']
        if p['recos']:
            leading = p['leading']
            total = info['reco_p4'][p['recos'], 3].sum()
            unit_gen = info['gen_p4'][g, :3] / np.linalg.norm(info['gen_p4'][g, :3])
            unit_reco = info['reco_p4'][leading, :3] / np.linalg.norm(info['reco_p4'][leading, :3])
            angle = np.arccos(np.clip(unit_gen @ unit_reco, -1, 1))
            lines.append(f'  leading=RecoPhoton[{leading}]; leading/gen={info["reco_p4"][leading, 3]/energy:.9g}; '
                         f'summed/gen={total/energy:.9g}; openingAngle={angle:.9g} rad')
    lines += ['', 'FULL STORED SIM FOREST (no species, energy, angle, or depth filtering)',
              'Indentation = physical Sim parentage. firstGen/depth = analysis origin and parent steps.',
              'A nested Gen anchor changes firstGen; descendants below it belong to that Gen lineage.']
    visited = set()

    def walk(s, prefix):
        assert s not in visited, 'Sim cycle or duplicate tree node'
        visited.add(s)
        parent = info['parents'][s]
        boundary = ' [GEN ANCHOR]' if raw['SimPart_genIdx'][s] >= 0 else ''
        lines.append(f'{prefix}Sim[{s}] mass_code={raw["SimPart_pdgId"][s]} Q={raw["SimPart_charge"][s]} '
                     f'genIdx={raw["SimPart_genIdx"][s]} firstGen={info["origins"][s]} depth={info["depths"][s]} '
                     f'parentSim={parent} originVtx={raw["SimPart_originVtxIdx"][s]} '
                     f'decayVtx={raw["SimPart_decayVtxIdx"][s]} partIdx={raw["SimPart_partIdx"][s]}{boundary}')
        lines.append(f'{prefix}  {p4_text(info["sim_p4"][s])}')
        for r in info['reco_at_sim'][s]:
            g = info['origins'][s]
            leading = g in photons and photons[g]['leading'] == r
            lines.append(f'{prefix}  -> RecoPhoton[{r}] Photon_partIdx={raw["Photon_partIdx"][r]} '
                         f'firstGen={g}' + (' [LEADING]' if leading else '') + f' {p4_text(info["reco_p4"][r])}')
        for child in info['children'][s]:
            walk(child, prefix + '  ')

    for s, parent in enumerate(info['parents']):
        if parent < 0:
            walk(s, '')
    assert visited == set(range(len(info['parents']))), 'Sim forest does not cover all saved nodes'
    lines += ['', 'RECO ASSOCIATIONS (all Photon candidates, including those without truth)']
    for r, s in enumerate(info['reco_sims']):
        g = info['origins'][s] if s >= 0 else -1
        direct = raw['SimPart_genIdx'][s] if s >= 0 else -1
        label = photons[g]['origin'] if g in photons else 'unassociated' if g < 0 else 'other Gen origin (not stable photon)'
        lines.append(f'RecoPhoton[{r}] Part[{raw["Photon_partIdx"][r]}] -> Sim[{s}] '
                     f'directGen={direct} firstGen={g} ({label}) {p4_text(info["reco_p4"][r])}')
    lines += ['', 'AUDIT CHECKS: PASS',
              '- Sim parent/decay vertices are reciprocal.',
              '- Every saved Sim appears exactly once in the full forest; no cycles.',
              '- Forward child traversal agrees with first-Gen-ancestor assignment for every stable photon.',
              '- Stable-photon Gen forward links agree with direct Sim Gen anchors.',
              '- Every associated reco has one origin; multiple reco count as one Gen success.']
    (OUTPUT / filename).write_text('\n'.join(lines) + '\n')


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    index = [
        'Manual Gen -> Sim -> Reco audit: 10 events, two per photosFSR sample.',
        'Selection is deliberately representative, not random or an efficiency estimate.',
        'For each sample: entry 0 of the first sorted ROOT file, then the first distinct event meeting its target criterion.',
        'Stable Gen photon: PDG 22, status 1. ISR/Beam definitions reuse the photon study.',
        'SimPart_pdgId is a DELPHI mass code (21=photon), not a PDG ID.',
        'Secondary Sim gamma excludes directly Gen-linked anchors. Truth association includes all Sim species.',
        'The full physical Sim forest includes other Gen origins; firstGen marks analysis lineage boundaries.',
        'No energy/angular cuts or angular matching. All saved Gen/Sim/Reco records of selected events are printed.',
        'Gen -> Sim means saved lineage. No transport, energy-deposit, or interaction-position claim is made.',
        'The full photon study has no no_sim photons; no missing-lineage example is fabricated.', '',
    ]
    number = 0
    for sample, target in zip(SAMPLES, CASES, strict=True):
        paths = sorted((INPUT / f'20260828_100kTest_{sample}_photosFSR/final_root').glob('job_*/nanoaod.root'))
        found = False
        for path in paths:
            with uproot.open(path) as root_file:
                tree = root_file['Events']
                branches = {name: tree[name].array(library='ak').to_list() for name in AUDIT_BRANCHES}
                for entry in range(tree.num_entries):
                    baseline = path == paths[0] and entry == 0
                    info = event_info({name: values[entry] for name, values in branches.items()}, sample)
                    focus = list(info['photons']) if baseline else focus_photons(info, target)
                    if not baseline and not focus:
                        continue
                    reason = 'baseline: first entry, all stable photons' if baseline else target
                    number += 1
                    filename = f'{number:02d}_{sample}_{path.parent.name}_entry_{entry}.txt'
                    dump_event(info, sample, path, entry, reason, focus, filename)
                    index.append(f'{filename}: {reason}; focus Gen={focus}; run={info["raw"]["Event_runNumber"]}, event={info["raw"]["Event_evtNumber"]}; PASS')
                    print(index[-1], flush=True)
                    if not baseline:
                        found = True
                        break
            if found:
                break
        assert found, f'No event found for {sample}: {target}'
    assert number == 10
    (OUTPUT / 'README.txt').write_text('\n'.join(index) + '\n')
    print(f'audit: {OUTPUT}')


if __name__ == '__main__':
    main()
