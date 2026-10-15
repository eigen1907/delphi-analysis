"""Train an event classifier using a separate validation split for early stopping."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
from xgboost import XGBClassifier

from common.paths import DATASET_ROOT, MODEL_ROOT, PROJECT_ROOT
from event_classification.bdt import DEFAULT_CONFIG_PATH
from event_classification.dataset import feature_metadata, load_metadata, load_dataset_split


DEFAULT_INPUT = DATASET_ROOT / "20261015_jongwon"
DEFAULT_OUTPUT = MODEL_ROOT / "20261015_jongwon_detector_pid_default"
DEFAULT_HYPERPARAMETERS = PROJECT_ROOT / "config" / "pid" / "hyperparameters" / "default.json"


def load_hyperparameters(path: Path) -> tuple[str, dict]:
    with path.open() as input_file:
        config = json.load(input_file)

    name = config.get("name")
    parameters = config.get("parameters")
    if not isinstance(name, str) or not isinstance(parameters, dict):
        raise ValueError("Hyperparameter config needs string name and object parameters")
    reserved = {"num_class", "random_state"}
    conflicts = reserved & set(parameters)
    if conflicts:
        raise ValueError(
            f"Hyperparameter config cannot set runtime values: {sorted(conflicts)}"
        )
    return name, parameters


def train_model(
    input_dir: Path,
    output_dir: Path,
    hyperparameter_path: Path = DEFAULT_HYPERPARAMETERS,
    seed: int = 1907,
    feature_config_path: Path = DEFAULT_CONFIG_PATH,
    ncpu: int = 32,
) -> dict:
    metadata = feature_metadata(load_metadata(input_dir), feature_config_path)
    selected_features = list(metadata["features"])
    hyperparameter_name, model_parameters = load_hyperparameters(hyperparameter_path)
    model_parameters["n_jobs"] = ncpu
    train_features, train_labels = load_dataset_split(input_dir, "train", metadata)
    val_features, val_labels = load_dataset_split(input_dir, "val", metadata)

    model = XGBClassifier(
        num_class=len(metadata["labels"]),
        random_state=seed,
        **model_parameters,
    )
    model.fit(
        train_features,
        train_labels,
        eval_set=[(val_features, val_labels)],
        verbose=False,
    )

    bundle = {
        "model": model,
        "feature_set": metadata.get("feature_set", "detector_pid"),
        "features": selected_features,
        "labels": metadata["labels"],
        "feature_config": metadata["feature_config"],
        "hyperparameter_name": hyperparameter_name,
        "hyperparameters": model_parameters,
        "tree": metadata.get("tree"),
        "trees": metadata.get("trees"),
        "format_version": metadata["format_version"],
        "training_dataset": str(input_dir.resolve()),
        "training_campaign": metadata.get("campaign"),
        "dataset_seed": metadata.get("seed"),
        "feature_extraction": {
            "version": 2 if metadata["format_version"] == 6 else 1,
            "analysisLevel_features": metadata["feature_config"].get("analysisLevel_features", []),
            "haid_track_association": (
                "Use unique zero-based PA indices when analysisLevel.Part_paIdx exists; "
                "otherwise leave the configured track-associated HaID feature slots missing. "
                "Configured event-level HaID counts are retained."
            ),
        },
        "seed": seed,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    model_path = output_dir / "bdt.joblib"
    joblib.dump(bundle, model_path)
    print(f"model: {model_path}")
    return bundle


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, default=DEFAULT_INPUT,
                        help="BDT dataset directory")
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT,
                        help="trained model directory")
    parser.add_argument("--seed", type=int, default=1907, help="model seed (default: 1907)")
    parser.add_argument("--hyperparameters", type=Path, default=DEFAULT_HYPERPARAMETERS,
                        help="XGBoost hyperparameter JSON config")
    parser.add_argument("--feature-config", type=Path, default=DEFAULT_CONFIG_PATH,
                        help="RAW features; paired datasets also add five reconstructed event features")
    parser.add_argument("--ncpu", type=int, default=32, help="XGBoost CPU threads (default: 32)")
    args = parser.parse_args()
    try:
        if args.ncpu < 1:
            raise ValueError("--ncpu must be positive")
        train_model(args.input, args.output, args.hyperparameters, args.seed,
                    args.feature_config, args.ncpu)
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
