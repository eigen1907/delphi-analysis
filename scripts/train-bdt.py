#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import matplotlib
import numpy as np
import uproot

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))

from bdt_features import FEATURE_NAMES, TREE_NAME


def load_split(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with uproot.open(path) as root_file:
        tree = root_file[TREE_NAME]
        features = [
            np.asarray(tree[name].array(library="np"), dtype=np.float32)
            for name in FEATURE_NAMES
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
    features: np.ndarray,
    labels: np.ndarray,
    output_dir: Path,
    seed: int,
    max_events: int,
) -> None:
    from sklearn.inspection import permutation_importance

    if len(labels) > max_events:
        rng = np.random.default_rng(seed)
        indices = rng.choice(len(labels), size=max_events, replace=False)
        features = features[indices]
        labels = labels[indices]

    result = permutation_importance(
        model,
        features,
        labels,
        scoring="f1_macro",
        n_repeats=3,
        random_state=seed,
        n_jobs=1,
    )
    order = np.argsort(result.importances_mean)[::-1]

    csv_path = output_dir / "feature_importance.csv"
    with csv_path.open("w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(("feature", "importance_mean", "importance_std"))
        for index in order:
            writer.writerow(
                (
                    FEATURE_NAMES[index],
                    result.importances_mean[index],
                    result.importances_std[index],
                )
            )

    fig, ax = plt.subplots(figsize=(10, 8))
    selected = order[:20][::-1]
    ax.barh(
        [FEATURE_NAMES[index] for index in selected],
        result.importances_mean[selected],
        xerr=result.importances_std[selected],
    )
    ax.set_xlabel("Permutation importance (macro F1)")
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
        "--importance-events",
        type=int,
        default=20_000,
        help="maximum validation events for permutation importance (default: 20000)",
    )
    args = parser.parse_args()

    try:
        import joblib
        from sklearn.ensemble import HistGradientBoostingClassifier
    except ImportError as error:
        raise RuntimeError(
            "scikit-learn is required; recreate or update the environment from environment.yml"
        ) from error

    with (args.input / "metadata.json").open() as metadata_file:
        metadata = json.load(metadata_file)

    if tuple(metadata["features"]) != FEATURE_NAMES:
        raise ValueError("Prepared dataset feature list differs from bdt_features.FEATURE_NAMES")

    label_items = sorted(metadata["labels"].items(), key=lambda item: item[1])
    class_names = [name for name, _ in label_items]

    train_features, train_labels = load_split(args.input / "train.root")
    val_features, val_labels = load_split(args.input / "val.root")
    test_features, test_labels = load_split(args.input / "test.root")

    model = HistGradientBoostingClassifier(
        learning_rate=0.05,
        max_iter=300,
        max_leaf_nodes=31,
        l2_regularization=1.0,
        early_stopping=False,
        random_state=args.seed,
    )
    model.fit(train_features, train_labels)

    args.output.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": model,
            "features": list(FEATURE_NAMES),
            "labels": metadata["labels"],
        },
        args.output / "bdt.joblib",
    )

    metrics = {
        "model": {
            "type": type(model).__name__,
            "learning_rate": model.learning_rate,
            "max_iter": model.max_iter,
            "max_leaf_nodes": model.max_leaf_nodes,
            "l2_regularization": model.l2_regularization,
            "seed": args.seed,
        },
        "train": evaluate(model, train_features, train_labels, class_names),
        "val": evaluate(model, val_features, val_labels, class_names),
        "test": evaluate(model, test_features, test_labels, class_names),
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
        val_features,
        val_labels,
        args.output,
        args.seed,
        args.importance_events,
    )

    for split in ("train", "val", "test"):
        print(
            f"{split}: accuracy={metrics[split]['accuracy']:.4f}, "
            f"macro_f1={metrics[split]['macro_f1']:.4f}"
        )
    print(f"model and diagnostics: {args.output}")


if __name__ == "__main__":
    main()
