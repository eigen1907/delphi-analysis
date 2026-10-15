"""Project locations shared by the command-line programs."""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = PROJECT_ROOT / "data"
DATASET_ROOT = PROJECT_ROOT / "output" / "dataset"
PLOT_ROOT = PROJECT_ROOT / "output" / "plots"
MODEL_ROOT = PROJECT_ROOT / "output" / "models"
