"""Resolve the physics channel from campaign directory names."""
from pathlib import Path

SAMPLE_NAMES = ("Zee", "Zmumu", "Ztautau", "ZKK", "Zpipi")


def sample_name(path: str | Path) -> str:
    name = Path(path).name
    for token in name.split("_"):
        if token in SAMPLE_NAMES:
            return token
    return name.rsplit("_", 1)[-1]
