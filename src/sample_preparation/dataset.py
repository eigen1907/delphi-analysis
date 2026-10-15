"""Prepare analysis chunks, merged ROOT files and classification datasets."""
from __future__ import annotations

import argparse
from pathlib import Path

from common.paths import DATA_ROOT, DATASET_ROOT, PROJECT_ROOT
from .chunks import prepare_chunks
from .merge import merge_dataset


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, default=DATA_ROOT / "20261015_jongwon",
                        help="raw campaign root; prepared chunks for --stage merge/classification")
    parser.add_argument("-o", "--output", type=Path,
                        help="campaign output root (default: output/dataset/<campaign>)")
    parser.add_argument("--stage", choices=("all", "chunks", "merge", "classification"), default="all")
    parser.add_argument("--samples", nargs="+", help="sample names or full sample directory names")
    parser.add_argument("--backend", choices=("uproot", "hadd"), default="uproot")
    parser.add_argument("--check", type=Path, help="event filtering summary CSV")
    parser.add_argument("--feature-config", type=Path,
                        default=PROJECT_ROOT / "config/pid/features/pid.json")
    parser.add_argument("--fractions", nargs=3, type=float, default=(0.6, 0.2, 0.2),
                        metavar=("TRAIN", "VAL", "TEST"))
    parser.add_argument("--seed", type=int, default=1907)
    args = parser.parse_args()

    campaign = args.input.parent.name if args.input.name == "chunks" else args.input.name
    output_root = args.output or DATASET_ROOT / campaign
    chunks_root = output_root / "chunks" if args.stage in ("all", "chunks") else args.input
    if args.stage in ("all", "chunks"):
        prepare_chunks(args.input, chunks_root, samples=args.samples, check_path=args.check)
    if args.stage in ("all", "merge"):
        merge_dataset(chunks_root, output_root / "merged", samples=args.samples, backend=args.backend)
    if args.stage in ("all", "classification"):
        from event_classification.dataset import prepare_dataset
        import json
        with args.feature_config.open() as source:
            feature_set = json.load(source)["name"]
        prepare_dataset(chunks_root, output_root / "classification" / feature_set,
                        feature_config_path=args.feature_config,
                        fractions=tuple(args.fractions), seed=args.seed)


if __name__ == "__main__":
    main()
