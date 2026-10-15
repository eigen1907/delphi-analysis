#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import uproot

from event_classification.bdt import TREE_NAME, extract_features

from common.paths import DATASET_ROOT, MODEL_ROOT


DEFAULT_INPUT = DATASET_ROOT / "20261015_jongwon" / "merged"
DEFAULT_MODEL = MODEL_ROOT / "20261015_jongwon_pid_standard" / "bdt.joblib"
DEFAULT_OUTPUT = (
    DATASET_ROOT / "20261015_jongwon" / "predictions" / "20261015_jongwon_pid_standard"
)


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


def apply_model(input_path: Path, model_path: Path, output_path: Path) -> None:
    bundle = joblib.load(model_path)
    if input_path.suffix == ".root":
        output = (
            output_path if output_path.suffix == ".root"
            else output_path / f"prediction_{input_path.parent.name}.root"
        )
        predict(input_path, output, bundle)
    else:
        for sample in bundle["labels"]:
            predict(
                input_path / sample / "nanoaod_raw_sdst.root",
                output_path / f"prediction_{sample}.root",
                bundle,
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply a saved event classifier to RAW SDST data")
    parser.add_argument("-i", "--input", type=Path, default=DEFAULT_INPUT,
                        help="merged dataset directory or one RAW SDST ROOT file")
    parser.add_argument("-m", "--model", type=Path, default=DEFAULT_MODEL,
                        help="trained bdt.joblib")
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT,
                        help="prediction directory, or ROOT file for a single input")
    args = parser.parse_args()
    apply_model(args.input, args.model, args.output)


if __name__ == "__main__":
    main()
