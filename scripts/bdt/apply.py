#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import uproot


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[2]))

from bdt import TREE_NAME, extract_features


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path, help="input raw SDST ROOT file")
    parser.add_argument("-m", "--model", required=True, type=Path, help="trained bdt.joblib")
    parser.add_argument("-o", "--output", required=True, type=Path, help="prediction ROOT file")
    args = parser.parse_args()

    try:
        import joblib
    except ImportError as error:
        raise RuntimeError(
            "joblib and XGBoost are required; update the environment from environment.yml"
        ) from error

    bundle = joblib.load(args.model)
    feature_names = list(bundle["features"])
    feature_config = bundle["feature_config"]

    label_items = sorted(bundle["labels"].items(), key=lambda item: item[1])
    class_names = [name for name, _ in label_items]

    with uproot.open(args.input) as input_file:
        if TREE_NAME not in input_file:
            raise KeyError(f"Missing tree {TREE_NAME} in {args.input}")
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

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(args.output) as output_file:
        output_file[TREE_NAME] = output

    counts = np.bincount(predicted_label, minlength=len(class_names))
    for class_name, count in zip(class_names, counts, strict=True):
        print(f"{class_name}: {count}")
    print(f"predictions: {args.output}")


if __name__ == "__main__":
    main()
