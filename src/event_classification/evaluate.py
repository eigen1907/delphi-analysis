"""Evaluate a saved event classifier and write its metrics and diagnostic plots."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import joblib
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common.paths import DATASET_ROOT, MODEL_ROOT, PLOT_ROOT
from event_classification.dataset import load_metadata, load_split


DEFAULT_INPUT = DATASET_ROOT / "20261015_jongwon" / "classification" / "pid"
DEFAULT_MODEL = MODEL_ROOT / "20261015_jongwon_pid_standard" / "bdt.joblib"
DEFAULT_OUTPUT = PLOT_ROOT / "event_classification" / "20261015_jongwon_pid_standard"


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
    csv_dir: Path,
) -> None:
    importance = np.asarray(model.feature_importances_, dtype=float)
    order = np.argsort(importance)[::-1]

    csv_path = csv_dir / "feature_importance.csv"
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


def evaluate_model(input_dir: Path, model_path: Path, output_dir: Path) -> dict:
    metadata = load_metadata(input_dir)
    bundle = joblib.load(model_path)
    for key in ("features", "labels", "feature_config"):
        if metadata[key] != bundle[key]:
            raise ValueError(f"Model and evaluation dataset have different {key}")

    model = bundle["model"]
    selected_features = list(bundle["features"])
    label_items = sorted(bundle["labels"].items(), key=lambda item: item[1])
    class_names = [name for name, _ in label_items]
    best_iteration = int(getattr(model, "best_iteration", model.get_booster().num_boosted_rounds() - 1))
    metrics = {
        "model": {
            "type": type(model).__name__,
            "feature_set": bundle.get("feature_set", "pid"),
            "feature_count": len(selected_features),
            "hyperparameter_name": bundle["hyperparameter_name"],
            "parameters": model.get_params(),
            "best_iteration": best_iteration,
            "boosting_rounds_used": best_iteration + 1,
            "seed": bundle.get("seed", model.get_params().get("random_state")),
        },
    }
    for split in ("train", "val", "test"):
        features, labels = load_split(
            input_dir / f"{split}.root", selected_features, metadata["tree"]
        )
        metrics[split] = evaluate(model, features, labels, class_names)
        del features, labels
        print(
            f"{split}: accuracy={metrics[split]['accuracy']:.4f}, "
            f"macro_f1={metrics[split]['macro_f1']:.4f}"
        )

    metrics_path = model_path.parent / "metrics.json"
    with metrics_path.open("w") as output_file:
        json.dump(metrics, output_file, indent=2, sort_keys=True)
        output_file.write("\n")

    output_dir.mkdir(parents=True, exist_ok=True)
    for split in ("val", "test"):
        plot_confusion_matrix(
            metrics[split]["confusion_matrix"],
            class_names,
            output_dir / f"confusion_matrix_{split}.png",
        )
    write_feature_importance(model, selected_features, output_dir, model_path.parent)
    print(f"metrics: {metrics_path}")
    print(f"feature importance: {model_path.parent / 'feature_importance.csv'}")
    print(f"diagnostic plots: {output_dir}")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, default=DEFAULT_INPUT,
                        help="BDT dataset directory with train/val/test ROOT files")
    parser.add_argument("-m", "--model", type=Path, default=DEFAULT_MODEL,
                        help="trained bdt.joblib; metrics and importance CSV are saved beside it")
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT,
                        help="directory for diagnostic plots")
    args = parser.parse_args()
    try:
        evaluate_model(args.input, args.model, args.output)
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
