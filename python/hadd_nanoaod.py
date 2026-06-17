from __future__ import annotations

from pathlib import Path

import awkward as ak
import uproot

from utils import META_TREE_NAME, TREE_NAME


def read_tree(path: Path, tree_name: str) -> dict[str, ak.Array]:
    with uproot.open(path) as root_file:
        return root_file[tree_name].arrays(library="ak", how=dict)


def concatenate(chunks: list[dict[str, ak.Array]]) -> dict[str, ak.Array]:
    branches = list(chunks[0])
    return {branch: ak.concatenate([chunk[branch] for chunk in chunks], axis=0) for branch in branches}


def has_tree(path: Path, tree_name: str) -> bool:
    with uproot.open(path) as root_file:
        return tree_name in root_file or f"{tree_name};1" in root_file


def n_entries(arrays: dict[str, ak.Array]) -> int:
    return len(next(iter(arrays.values()))) if arrays else 0


def hadd_nanoaods(inputs: list[Path], output: Path) -> dict[str, int]:
    events = concatenate([read_tree(path, TREE_NAME) for path in inputs])
    metas = [read_tree(path, META_TREE_NAME) for path in inputs if has_tree(path, META_TREE_NAME)]

    output.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output) as root_file:
        root_file[TREE_NAME] = events
        if metas:
            root_file[META_TREE_NAME] = concatenate(metas)

    return {"files": len(inputs), "events": n_entries(events), "meta_entries": n_entries(concatenate(metas)) if metas else 0}
