from __future__ import annotations

from pathlib import Path

import mplhep as mh

from plot_isr_photons import SAMPLES, plot_distribution, plot_efficiency, plot_multiplicity, read_sample
from plot_utils import sample_styles


def plot_fsr_photons(input_root: Path, output_root: Path) -> None:
    mh.style.use(mh.styles.CMS)
    study_dir = output_root / "fsr_photons"
    gen_dir = study_dir / "01_gen"
    matching_dir = study_dir / "03_matching"
    efficiency_dir = study_dir / "04_efficiency"
    for directory in (gen_dir, matching_dir, efficiency_dir):
        directory.mkdir(parents=True, exist_ok=True)

    styles = sample_styles(SAMPLES)
    values_by_sample = {sample: read_sample(input_root, sample) for sample in SAMPLES}
    for kind, directory in (("fsr", gen_dir), ("fsr_matched_reco", matching_dir)):
        plot_multiplicity(directory, kind, values_by_sample, styles)
        for name in ("energy", "cos_theta", "phi"):
            plot_distribution(directory, kind, name, values_by_sample, styles)
    for name in ("energy", "cos_theta"):
        plot_efficiency(efficiency_dir, "fsr", name, values_by_sample, styles)
    print(f"plots: {study_dir}")
