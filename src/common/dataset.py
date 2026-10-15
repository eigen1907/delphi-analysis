"""Paths shared by paired analysis-level and raw SDST datasets."""
from pathlib import Path

from .samples import sample_name


SOURCE_TREES = {
    "nanoaod": "analysisLevel/Events",
    "nanoaod_raw_sdst": "raw_sdst/Events",
}
SPLITS = ("train", "val", "test")


def dataset_files(input_root: Path, split: str = "train", samples=None):
    """Return one paired split file per sample, from a campaign or sample path."""
    input_root = Path(input_root)
    if split not in SPLITS:
        raise ValueError(f"Unknown split: {split}")
    if input_root.is_file():
        paths = [input_root]
    elif (input_root / f"{split}.root").is_file():
        paths = [input_root / f"{split}.root"]
    else:
        paths = sorted(input_root.glob(f"*/{split}.root"))
    requested = set(samples or ())
    result = [(path.parent.name, [path]) for path in paths
              if not requested or path.parent.name in requested
              or sample_name(path.parent) in requested]
    if not result:
        raise FileNotFoundError(f"No {split} datasets found under {input_root}")
    return result
