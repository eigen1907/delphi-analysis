#!/usr/bin/env python3
"""Preserve representative truth trees and audit ISR Part/conversion links."""
from pathlib import Path

import numpy as np
import uproot

from delphi_analysis.isr_study.data import (
    BRANCHES, GEN_P4, RECO_P4, SAMPLES, SIM_GAMMA_CODE, TARGET_PDG,
    first_gen_ancestor, photon_origin, analyze_event,
)

PROJECT = Path(__file__).resolve().parents[2]
INPUT = PROJECT / 'data/20260828_florian'
OUTPUT = PROJECT / 'plots/20260828_florian/isr_study/02_truth_link_matching_validation/manual_audit'
SIM_P4 = tuple(f'SimPart_fourMomentum.fCoordinates.f{axis}' for axis in 'XYZT')
PART_P4 = tuple(f'Part_fourMomentum.fCoordinates.f{axis}' for axis in 'XYZT')
AUDIT_BRANCHES = tuple(dict.fromkeys((*BRANCHES, *SIM_P4, 'Event_runNumber', 'Event_evtNumber',
                  'GenPart_simIdx', 'SimPart_decayVtxIdx', 'SimPart_charge', 'SimPart_partIdx', 'Part_lock')))
CASES = ('direct ISR reco', 'descendant ISR reco', 'nonprimary photon anchor',
         'split ISR', 'ISR shower with secondary gamma')
# Fixed cases from the preceding full-sample diagnostics (zero-based entries).
LINK_CASES = (
    ('Zee', 0, 17, 'conversion: two agreeing daughters'),
    ('Zee', 0, 54, 'conversion: parent-only forward truth'),
    ('Zee', 0, 4131, 'reverse-only ISR Photon associations'),
    ('Ztautau', 12, 2730, 'conversion: contradictory ISR/pion anchors'),
    ('Zpipi', 10, 4492, 'conversion: contradictory ISR/pion anchors'),
    ('Zpipi', 15, 483, 'conversion: contradictory ISR/pion anchors'),
    ('Zpipi', 14, 3870, 'suspicious Beam ISR association'),
    ('Zpipi', 17, 2420, 'suspicious Beam ISR association'),
)


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


def part_audit(info, sample):
    """Independently check both Part-link directions against the nominal union."""
    raw = info['raw']
    result = analyze_event(raw, sample)
    reverse = [set() for _ in raw['Part_simIdx']]
    for sim, part in enumerate(raw['SimPart_partIdx']):
        if part >= 0 and info['origins'][sim] >= 0:
            reverse[part].add(int(info['origins'][sim]))
    lines = ['', 'PART LINK VALIDATION (all Gen anchors, including non-ISR)',
             'forward = Part->Sim ancestry; reverse = Sim ancestry->Part.',
             'A unique union anchor is accepted; contradictory anchors remain unresolved.']
    for part, sim in enumerate(raw['Part_simIdx']):
        forward = int(info['origins'][sim]) if sim >= 0 else -1
        known = reverse[part] | ({forward} if forward >= 0 else set())
        expected = next(iter(known)) if len(known) == 1 else -1
        assert set(result['part_origin_sets'][part]) == known
        assert result['part_origins'][part] == expected, 'Incorrect conflict-free Part union'
        p4 = np.asarray([raw[name][part] for name in PART_P4])
        lines.append(f'Part[{part}] mass_code={raw["Part_pdgId"][part]} Q={raw["Part_charge"][part]} '
                     f'lock={raw["Part_lock"][part]} simIdx={sim} forwardGen={forward} '
                     f'reverseGen={sorted(reverse[part])} unionGen={expected} '
                     f'originVtx={raw["Part_originVtxIdx"][part]} decayVtx={raw["Part_decayVtxIdx"][part]} '
                     f'{p4_text(p4)}')
    lines += ['', 'PHOTON CONVERSION ASSOCIATIONS',
              'The saved neutral conversion parent and known daughter/parent anchors must agree.']
    for conv, origin in enumerate(result['conversion_origins']):
        parent = result['conversion_parents'][conv]
        daughters = result['conversion_daughters'][conv]
        members = [p for p in (parent, *daughters) if p >= 0]
        known = set().union(*(result['part_origin_sets'][p] for p in members))
        if len(known) > 1:
            assert origin < 0, 'Contradictory conversion was accepted'
        lines.append(f'PhotonConv[{conv}] simPhotonIdx={raw["PhotonConv_simPhotonIdx"][conv]} '
                     f'parentPart={parent} daughterParts={sorted(daughters)} '
                     f'knownGen={sorted(known)} acceptedGen={origin}')
    lines += ['', 'NOMINAL ISR OBSERVABLES (one row per Gen ISR)']
    for row, gen in enumerate(result['gen_idx']):
        parts = np.flatnonzero(result['part_origins'] == gen).tolist()
        conversions = np.flatnonzero(result['conversion_origins'] == gen).tolist()
        lines.append(f'Gen[{gen}] unionParts={parts} acceptedConversions={conversions}')
        for name, channel in result['channels'].items():
            stats = channel['stats']
            lines.append(f'  {name}: Reco={stats["reco_count"][row]} raw={stats["raw_count"][row]} '
                         f'canonical={stats["energy_count"][row]} energyValid={stats["energy_valid"][row]} '
                         f'energy/gen={stats["energy_ratio"][row]:.9g} '
                         f'rawEnergy/gen={stats["raw_energy_ratio"][row]:.9g}')
    if 'contradictory ISR/pion' in info.get('reason', ''):
        assert result['conversion_origins'][0] < 0, 'Known contradictory ISR conversion was accepted'
    if info.get('reason') == 'conversion: two agreeing daughters':
        assert result['conversion_origins'][0] == 17 and result['conversion_parents'][0] == 4
    if info.get('reason') == 'conversion: parent-only forward truth':
        assert result['conversion_origins'][0] == 8 and result['conversion_parents'][0] == 6
    if info.get('reason') == 'reverse-only ISR Photon associations':
        assert np.all(result['forward_origins'][[2, 3]] == -1)
        assert np.all(result['part_origins'][[2, 3]] == 13)
    return result, lines


def dump_event(info, sample, path, entry, reason, focus, filename):
    raw, photons = info['raw'], info['photons']
    info['reason'] = reason
    _, part_lines = part_audit(info, sample)
    lines = [f'Sample: {sample} photosFSR', f'ROOT: {path.relative_to(PROJECT)}',
             f'Entry: {entry} (zero-based)', f'Event_runNumber: {raw["Event_runNumber"]}',
             f'Event_evtNumber: {raw["Event_evtNumber"]}', f'Selection: {reason}',
             f'Focus Gen indices: {focus}', '', 'GEN RECORD (all particles; indices are event-local)']
    for g, p4 in enumerate(info['gen_p4']):
        lines.append(f'Gen[{g}] pdg={raw["GenPart_pdgId"][g]} status={raw["GenPart_status"][g]} '
                     f'parent={raw["GenPart_parentIdx"][g]} simIdx={raw["GenPart_simIdx"][g]} {p4_text(p4)}')
    lines += ['', 'STABLE GEN PHOTONS (one-direction Photon associations retained for comparison)',
              'Nominal ISR results use the conflict-free Part union printed below.']
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
    lines += part_lines
    lines += ['', 'AUDIT CHECKS: PASS',
              '- Sim parent/decay vertices are reciprocal.',
              '- Every saved Sim appears exactly once in the full forest; no cycles.',
              '- Forward child traversal agrees with first-Gen-ancestor assignment for every stable photon.',
              '- Stable-photon Gen forward links agree with direct Sim Gen anchors.',
              '- Every associated reco has one origin; multiple reco count as one Gen success.',
              '- Independent forward/reverse Part origin sets agree with the nominal conflict-free union.',
              '- Conflicting known conversion anchors are not accepted.']
    (OUTPUT / filename).write_text('\n'.join(lines) + '\n')


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    index = [
        'Manual Gen -> Sim -> Reco audit: 10 preserved cases plus 8 fixed link/conversion checks.',
        'Selection is deliberately representative, not random or an efficiency estimate.',
        'For each sample: entry 0 of the first sorted ROOT file, then the first distinct event meeting its target criterion.',
        'Stable Gen photon: PDG 22, status 1. ISR/Beam definitions reuse isr_study.',
        'SimPart_pdgId is a DELPHI mass code (21=photon), not a PDG ID.',
        'Secondary Sim gamma excludes directly Gen-linked anchors. Truth association includes all Sim species.',
        'The full physical Sim forest includes other Gen origins; firstGen marks analysis lineage boundaries.',
        'No energy/angular cuts or angular matching. All saved Gen/Sim/Reco records of selected events are printed.',
        'Gen -> Sim means saved lineage. No transport, energy-deposit, or interaction-position claim is made.',
        'The old Photon association is printed for comparison; Part sections validate the new nominal union.',
        'No unknown truth is changed to an unrelated origin or a zero-energy observation.', '',
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
    for sample, job, entry, reason in LINK_CASES:
        path = INPUT / f'20260828_100kTest_{sample}_photosFSR/final_root/job_{job}/nanoaod.root'
        with uproot.open(path) as root_file:
            raw = {name: root_file['Events'][name].array(entry_start=entry, entry_stop=entry + 1,
                                                       library='ak').to_list()[0] for name in AUDIT_BRANCHES}
        info = event_info(raw, sample)
        focus = [g for g, p in info['photons'].items() if 'ISR' in p['origin']]
        number += 1
        filename = f'{number:02d}_{sample}_job_{job}_entry_{entry}.txt'
        dump_event(info, sample, path, entry, reason, focus, filename)
        index.append(f'{filename}: {reason}; focus Gen={focus}; PASS')
        print(index[-1], flush=True)
    assert number == 18
    (OUTPUT / 'README.txt').write_text('\n'.join(index) + '\n')
    print(f'audit: {OUTPUT}')


if __name__ == '__main__':
    main()
