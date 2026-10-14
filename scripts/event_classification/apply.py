#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import uproot

from event_classification.bdt import TREE_NAME, extract_features

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def predict(input_path: Path, output_path: Path, bundle: dict) -> None:
    feature_names = list(bundle["features"])
    feature_config = bundle["feature_config"]

    label_items = sorted(bundle["labels"].items(), key=lambda item: item[1])
    class_names = [name for name, _ in label_items]

    with uproot.open(input_path) as input_file:
        if TREE_NAME not in input_file:
            raise KeyError(f"Missing tree {TREE_NAME} in {input_path}")
        tree = input_file[TREE_NAME]
        features = extract_features(tree, feature_config)
        if list(features) != feature_names:
            raise ValueError("Model and extracted feature lists differ")
        feature_matrix = np.column_stack([features[name] for name in feature_names])
        predicted_label = np.asarray(bundle["model"].predict(feature_matrix), dtype=np.int8)
        probabilities = np.asarray(bundle["model"].predict_proba(feature_matrix), dtype=np.float32)
        probability_column = {
            int(label): column
            for column, label in enumerate(bundle["model"].classes_)
        }
        output = {
            "run_number": np.asarray(tree["Event_runNumber"].array(library="np"), dtype=np.int32),
            "event_number": np.asarray(tree["Event_eventNumber"].array(library="np"), dtype=np.int32),
            "predicted_label": predicted_label,
            **{
                f"probability_{class_name}": probabilities[:, probability_column[label]]
                for class_name, label in label_items
            },
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output_path) as output_file:
        output_file[TREE_NAME] = output

    counts = np.bincount(predicted_label, minlength=len(class_names))
    for class_name, count in zip(class_names, counts, strict=True):
        print(f"{class_name}: {count}")
    print(f"predictions: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=Path,
                        default=PROJECT_ROOT / "data/202606xx_jongwon/dataset/20260606_100kTest",
                        help="dataset directory or one raw SDST ROOT file")
    parser.add_argument("-m", "--model", type=Path,
                        default=PROJECT_ROOT / "plots/event_classification/20260606_100kTest_pid_standard/bdt.joblib",
                        help="trained bdt.joblib")
    parser.add_argument("-o", "--output", type=Path,
                        default=PROJECT_ROOT / "data/202606xx_jongwon/ml/20260606_100kTest_pid_standard",
                        help="prediction directory, or ROOT file for a single input")
    args = parser.parse_args()
    bundle = joblib.load(args.model)
    if args.input.suffix == ".root":
        output = args.output if args.output.suffix == ".root" else args.output / f"prediction_{args.input.parent.name}.root"
        predict(args.input, output, bundle)
    else:
        for sample in bundle["labels"]:
            predict(args.input / sample / "nanoaod_raw_sdst.root",
                    args.output / f"prediction_{sample}.root", bundle)


if __name__ == "__main__":
    main()
