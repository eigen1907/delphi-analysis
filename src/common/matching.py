"""Small matching functions shared by the photon and charged-track studies."""

import numpy as np


def first_gen_ancestor(sim, sim_gens, sim_vertices, vertex_incoming):
    """Follow Sim mothers and stop at the first saved Gen anchor."""
    seen = set()
    while 0 <= sim < len(sim_gens):
        if sim in seen:
            raise ValueError("Cycle in Sim ancestry")
        seen.add(sim)
        if sim_gens[sim] >= 0:
            return int(sim_gens[sim])
        vertex = sim_vertices[sim]
        sim = vertex_incoming[vertex] if 0 <= vertex < len(vertex_incoming) else -1
    return -1


def descendants(index, children):
    """Saved descendants, without inferring relationships from geometry."""
    found = set()
    pending = list(children[index])
    while pending:
        child = pending.pop()
        if child == index or child in found:
            raise ValueError("Cycle or duplicate edge in saved ancestry")
        found.add(child)
        pending.extend(children[child])
    return found


def link_status(forward, reverse):
    if len(forward | reverse) > 1:
        return "conflict"
    if forward and reverse:
        return "agreement"
    if forward:
        return "forward_only"
    if reverse:
        return "reverse_only"
    return "no_origin"


def resolve_part_origins(raw):
    """Compare both saved truth-link directions before applying any Gen cuts.

    Forward: Part -> Sim -> first Gen anchor.
    Reverse: Gen-anchored Sim -> descendants -> Part, stopping at new anchors.
    A unique origin from either direction is usable; conflicting origins stay -1.
    """
    sim_gens = np.asarray(raw["SimPart_genIdx"], dtype=int)
    vertices = raw["SimPart_originVtxIdx"]
    incoming = raw["SimVtx_incomingIdx"]
    n_gen, n_part = len(raw["GenPart_pdgId"]), len(raw["Part_simIdx"])
    children = [[] for _ in sim_gens]
    for sim, vertex in enumerate(vertices):
        mother = incoming[vertex] if 0 <= vertex < len(incoming) else -1
        if 0 <= mother < len(sim_gens):
            children[mother].append(sim)

    reverse = [set() for _ in range(n_part)]
    for anchor in np.flatnonzero((sim_gens >= 0) & (sim_gens < n_gen)):
        gen, pending, seen = int(sim_gens[anchor]), [int(anchor)], set()
        while pending:
            sim = pending.pop()
            if sim in seen:
                raise ValueError("Cycle in downward Sim ancestry")
            seen.add(sim)
            if sim_gens[sim] >= 0 and sim_gens[sim] != gen:
                continue
            part = raw["SimPart_partIdx"][sim]
            if 0 <= part < n_part:
                reverse[part].add(gen)
            pending.extend(children[sim])

    forward = np.full(n_part, -1, dtype=int)
    resolved = {}
    for part, sim in enumerate(raw["Part_simIdx"]):
        if 0 <= sim < len(sim_gens):
            if sim not in resolved:
                resolved[sim] = first_gen_ancestor(sim, sim_gens, vertices, incoming)
            gen = resolved[sim]
            if 0 <= gen < n_gen:
                forward[part] = gen
    known = [origins | ({int(forward[part])} if forward[part] >= 0 else set())
             for part, origins in enumerate(reverse)]
    origin = np.asarray([next(iter(origins)) if len(origins) == 1 else -1
                         for origins in known], dtype=int)
    return dict(forward=forward, reverse=reverse, known=known, origin=origin)


def angle_matrix(gen_xyz, reco_xyz):
    """3D opening angles in radians; undefined directions have infinite distance."""
    gen_xyz, reco_xyz = np.asarray(gen_xyz), np.asarray(reco_xyz)
    gen_p = np.linalg.norm(gen_xyz, axis=1)
    reco_p = np.linalg.norm(reco_xyz, axis=1)
    gen_unit = np.divide(gen_xyz, gen_p[:, None], out=np.zeros_like(gen_xyz, dtype=float),
                         where=gen_p[:, None] > 0)
    reco_unit = np.divide(reco_xyz, reco_p[:, None], out=np.zeros_like(reco_xyz, dtype=float),
                          where=reco_p[:, None] > 0)
    angles = np.arccos(np.clip(gen_unit @ reco_unit.T, -1, 1))
    valid = (np.isfinite(gen_xyz).all(axis=1) & (gen_p > 0))[:, None]
    valid = valid & (np.isfinite(reco_xyz).all(axis=1) & (reco_p > 0))[None, :]
    return np.where(valid, angles, np.inf)


def match_tracks(gen_xyz, reco_xyz, gen_charge, reco_charge, cut_radians):
    """Same-charge, smallest-angle-first one-to-one track matching."""
    angles = angle_matrix(gen_xyz, reco_xyz)
    eligible = (angles < cut_radians) & (np.asarray(gen_charge)[:, None] == np.asarray(reco_charge)[None, :])
    candidates = np.argwhere(eligible)
    ordered = sorted(candidates, key=lambda pair: angles[tuple(pair)])
    used_gen, used_reco, pairs = set(), set(), []
    for gen, reco in ordered:
        if gen not in used_gen and reco not in used_reco:
            used_gen.add(gen)
            used_reco.add(reco)
            pairs.append((gen, reco))
    return np.asarray(pairs, dtype=int).reshape(-1, 2).T
