#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import joblib
import matplotlib
import numpy as np
import uproot
from xgboost import XGBClassifier

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]


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


def load_split(
    path: Path,
    feature_names: list[str],
    tree_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    with uproot.open(path) as root_file:
        tree = root_file[tree_name]
        features = [
            np.asarray(tree[name].array(library="np"), dtype=np.float32)
            for name in feature_names
        ]
        labels = np.asarray(tree["label"].array(library="np"), dtype=np.int8)
    return np.column_stack(features), labels


def evaluate(model, features, labels, class_names):
    from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, f1_score

    prediction = model.predict(features)
    return {
        "n_events": int(len(labels)),
        "accuracy": float(accuracy_score(labels, prediction)),
        "macro_f1": float(f1_score(labels, prediction, average="macro")),
        "classification_report": classification_report(
            labels,
            prediction,
            labels=np.arange(len(class_names)),
            target_names=class_names,
            output_dict=True,
            zero_division=0,
        ),
        "confusion_matrix": confusion_matrix(
            labels,
            prediction,
            labels=np.arange(len(class_names)),
        ).tolist(),
    }


def plot_confusion_matrix(matrix, class_names, output_path: Path) -> None:
    values = np.asarray(matrix, dtype=float)
    row_sum = values.sum(axis=1, keepdims=True)
    normalized = np.divide(values, row_sum, out=np.zeros_like(values), where=row_sum > 0)

    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(normalized, vmin=0.0, vmax=1.0, cmap="Blues")
    fig.colorbar(image, ax=ax, label="Fraction")
    ax.set_xticks(range(len(class_names)), class_names)
    ax.set_yticks(range(len(class_names)), class_names)
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")

    for row in range(len(class_names)):
        for column in range(len(class_names)):
            value = normalized[row, column]
            color = "white" if value > 0.5 else "black"
            ax.text(column, row, f"{value:.3f}", ha="center", va="center", color=color)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def write_feature_importance(
    model,
    feature_names: list[str],
    output_dir: Path,
) -> None:
    importance = np.asarray(model.feature_importances_, dtype=float)
    order = np.argsort(importance)[::-1]

    csv_path = output_dir / "feature_importance.csv"
    with csv_path.open("w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(("feature", "gain_importance"))
        for index in order:
            writer.writerow((feature_names[index], importance[index]))

    fig, ax = plt.subplots(figsize=(10, 8))
    selected = order[:20][::-1]
    ax.barh(
        [feature_names[index] for index in selected],
        importance[selected],
    )
    ax.set_xlabel("XGBoost gain importance")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / "feature_importance.png", dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path, help="BDT dataset directory")
    parser.add_argument("-o", "--output", required=True, type=Path, help="model output directory")
    parser.add_argument("--seed", type=int, default=1907, help="model seed (default: 1907)")
    parser.add_argument(
        "--hyperparameters",
        type=Path,
        default=PROJECT_ROOT
        / "config"
        / "bdt"
        / "hyperparameters"
        / "standard.json",
        help="XGBoost hyperparameter JSON config",
    )
    args = parser.parse_args()

    with (args.input / "metadata.json").open() as metadata_file:
        metadata = json.load(metadata_file)
    if metadata.get("format_version") != 5:
        parser.error(
            "BDT dataset uses an older feature format; rerun scripts/bdt/prepare.py"
        )

    feature_names = list(metadata["features"])
    feature_set = metadata.get("feature_set", "pid")
    feature_config = metadata["feature_config"]
    selected_tree = metadata["tree"]
    label_items = sorted(metadata["labels"].items(), key=lambda item: item[1])
    class_names = [name for name, _ in label_items]
    try:
        hyperparameter_name, model_parameters = load_hyperparameters(
            args.hyperparameters
        )
    except ValueError as error:
        parser.error(str(error))

    train_features, train_labels = load_split(
        args.input / "train.root", feature_names, selected_tree
    )
    val_features, val_labels = load_split(
        args.input / "val.root", feature_names, selected_tree
    )

    model = XGBClassifier(
        num_class=len(class_names),
        random_state=args.seed,
        **model_parameters,
    )
    model.fit(
        train_features,
        train_labels,
        eval_set=[(val_features, val_labels)],
        verbose=False,
    )

    args.output.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": model,
            "feature_set": feature_set,
            "features": feature_names,
            "labels": metadata["labels"],
            "feature_config": feature_config,
            "hyperparameter_name": hyperparameter_name,
            "hyperparameters": model_parameters,
        },
        args.output / "bdt.joblib",
    )

    train_metrics = evaluate(model, train_features, train_labels, class_names)
    val_metrics = evaluate(model, val_features, val_labels, class_names)
    del train_features, train_labels, val_features, val_labels
    test_features, test_labels = load_split(
        args.input / "test.root", feature_names, selected_tree
    )
    test_metrics = evaluate(model, test_features, test_labels, class_names)

    metrics = {
        "model": {
            "type": type(model).__name__,
            "feature_set": feature_set,
            "feature_count": len(feature_names),
            "hyperparameter_name": hyperparameter_name,
            "parameters": model.get_params(),
            "best_iteration": int(model.best_iteration),
            "boosting_rounds_used": int(model.best_iteration + 1),
            "seed": args.seed,
        },
        "train": train_metrics,
        "val": val_metrics,
        "test": test_metrics,
    }
    with (args.output / "metrics.json").open("w") as output_file:
        json.dump(metrics, output_file, indent=2, sort_keys=True)
        output_file.write("\n")

    for split in ("val", "test"):
        plot_confusion_matrix(
            metrics[split]["confusion_matrix"],
            class_names,
            args.output / f"confusion_matrix_{split}.png",
        )

    write_feature_importance(
        model,
        feature_names,
        args.output,
    )

    for split in ("train", "val", "test"):
        print(
            f"{split}: accuracy={metrics[split]['accuracy']:.4f}, "
            f"macro_f1={metrics[split]['macro_f1']:.4f}"
        )
    print(f"model and diagnostics: {args.output}")


if __name__ == "__main__":
    main()
