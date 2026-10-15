"""Storage checks stay separate from the saved-angle coverage study."""
from collections import Counter

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import uproot

from common.plot_utils import FIGURE_SIZE, clopper_pearson, finish, plot_path
from .data import FIELDS, SOURCES
from .plot_helpers import histogram

FIELD_LABELS = (r"$Q$", r"$\theta_g$", r"$\sigma_g$", r"$N_g$", r"$N_{\mathrm{exp},g}$", r"$F_g$",
                r"$\theta_l$", r"$\sigma_l$", r"$N_l$", r"$N_{\mathrm{exp},l}$", r"$F_l$")


def storage(values, files, source, output, sample):
    with uproot.open(files[0]) as root:
        types = root["Events"].typenames()
    normal, subnormal, zero = [], [], []
    for name, entries in values["stored"].items():
        entries = np.asarray(entries)
        finite = np.isfinite(entries)
        small = finite & (entries != 0) & (abs(entries) < np.finfo(np.float32).tiny)
        integer_like = finite & (entries != 0) & (entries == np.round(entries))
        print(f"{sample} {source} {FIELDS[source][name]}: {types[FIELDS[source][name]]}; "
              f"entries={len(entries)}, nonfinite={(~finite).sum()}, subnormal={small.sum()}, "
              f"zero={(entries == 0).sum()}, nonzero_integer_like={integer_like.sum()}")
        normal.append(np.count_nonzero(finite & ~small & (entries != 0)) / max(len(entries), 1))
        subnormal.append(small.sum() / max(len(entries), 1))
        zero.append(np.count_nonzero(entries == 0) / max(len(entries), 1))
    x = np.arange(len(normal))
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.bar(x, normal, label="Normal nonzero", color="C0")
    ax.bar(x, subnormal, bottom=normal, label="Subnormal", color="C3")
    ax.bar(x, zero, bottom=np.asarray(normal) + subnormal, label="Zero", color="0.8")
    ax.set_xticks(x, FIELD_LABELS)
    ax.set_ylim(0, 1.7)
    ax.legend(loc="upper right", frameon=False)
    finish(fig, ax, output / "storage_status.png", sample, "Saved field", "Fraction")
    for name in ("gas_sigma", "liquid_expected"):
        entries = np.asarray(values["stored"][name], dtype=np.float32)
        small = np.isfinite(entries) & (entries != 0) & (abs(entries) < np.finfo(np.float32).tiny)
        if not small.any():
            continue
        bits = entries[small].view(np.uint32)
        codes, counts = np.unique(bits, return_counts=True)
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        ax.bar(np.arange(len(codes)), counts, yerr=np.sqrt(counts), color="C3", capsize=2)
        shown = np.linspace(0, len(codes) - 1, min(len(codes), 10), dtype=int)
        ax.set_xticks(shown, codes[shown])
        ax.set_ylim(0, counts.max() * 1.3)
        finish(fig, ax, output / f"{name}_raw_bits.png", sample, "Subnormal float32 bits", "Saved entries")
    print(f"{sample} {source}: {values['counts']}")
    for name in ("gas_angle", "liquid_angle", "gas_photons", "liquid_photons"):
        entries = np.asarray(values["stored"][name])
        entries = entries[np.isfinite(entries)]
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        histogram(ax, entries, np.histogram_bin_edges(entries, bins=50), "Stored field", "C0", 1)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.4)
        ax.text(0.02, 0.86, "Saved RICH fields", transform=ax.transAxes, fontsize=18)
        xlabel = r"Saved $\theta_C$ [rad]" if name.endswith("angle") else r"Saved $N_{\mathrm{ph}}$"
        finish(fig, ax, output / f"{name}.png", sample, xlabel, "Saved entries / bin")


def source_comparison(files, output, sample):
    """Pair by run/event and original PA identity, never by VECP row position."""
    equal, compared = Counter(), Counter()
    gen_equal = Counter()
    events = Counter()
    pairs = list(zip(FIELDS["nanoaod_raw_sdst"].values(), FIELDS["nanoaod"].values(), strict=True))
    gen_pairs = [("GenPart_pdgId", "GenPart_pdgId"), ("GenPart_status", "GenPart_status")]
    gen_pairs += [(f"GenPart_fourMomentum.fCoordinates.f{axis}", f"GenPart_vector.fCoordinates.f{axis}")
                  for axis in "XYZT"]
    for path in files:
        reference = path.with_name("nanoaod.root")
        if not reference.exists():
            continue
        with uproot.open(path) as raw_file, uproot.open(reference) as ref_file:
            if ("Part_paIdx" not in ref_file["Events"] or
                    not raw_file["Events"]["HaidRaw_paIdx"].description.startswith("Zero-based PA index")):
                print(f"{sample}: RICH source comparison skipped; old files lack a common PA identity.")
                continue
            raw_names = ["Event_runNumber", "Event_eventNumber", "HaidRaw_paIdx",
                         *FIELDS["nanoaod_raw_sdst"].values(), *(a for a, _ in gen_pairs)]
            ref_names = ["Event_runNumber", "Event_evtNumber", "Part_paIdx",
                         *FIELDS["nanoaod"].values(), *(b for _, b in gen_pairs)]
            # Read dotted four-vector leaves directly; arrays() can omit them.
            raw = {name: raw_file["Events"][name].array(library="ak") for name in raw_names}
            ref = {name: ref_file["Events"][name].array(library="ak") for name in ref_names}
            raw_keys = list(zip(ak.to_list(raw["Event_runNumber"]), ak.to_list(raw["Event_eventNumber"]), strict=True))
            ref_keys = list(zip(ak.to_list(ref["Event_runNumber"]), ak.to_list(ref["Event_evtNumber"]), strict=True))
            repeated = {key for key, count in Counter(ref_keys).items() if count != 1}
            lookup = {key: row for row, key in enumerate(ref_keys) if key not in repeated}
            events["raw"] += len(raw_keys)
            events["reference"] += len(ref_keys)
            for row, key in enumerate(raw_keys):
                # Run/setup records can retain a previous physics event number.
                if len(raw["GenPart_pdgId"][row]) == 0:
                    continue
                if key not in lookup:
                    events["unresolved_event"] += 1
                    continue
                other = lookup[key]
                events["aligned"] += 1
                for a, b in gen_pairs:
                    gen_equal[a] += int(np.array_equal(np.asarray(raw[a][row]), np.asarray(ref[b][other]), equal_nan=True))
                ids = np.asarray(raw["HaidRaw_paIdx"][row], dtype=int)
                ids_count = Counter(ids)
                ref_ids = np.asarray(ref["Part_paIdx"][other], dtype=int)
                ref_counts = Counter(ref_ids)
                pa_lookup = {part: index for index, part in enumerate(ref_ids)
                             if part >= 0 and ref_counts[part] == 1}
                r = np.asarray([index for index, part in enumerate(ids)
                                if part in pa_lookup and ids_count[part] == 1], dtype=int)
                q = np.asarray([pa_lookup[part] for part in ids[r]], dtype=int)
                for name, reference_name in pairs:
                    a, b = np.asarray(raw[name][row])[r], np.asarray(ref[reference_name][other])[q]
                    compared[name] += len(a)
                    equal[name] += int(np.isclose(a, b, rtol=0, atol=0, equal_nan=True).sum())
    if not compared:
        return
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    labels = list(FIELDS["nanoaod_raw_sdst"])
    fractions = [equal[FIELDS["nanoaod_raw_sdst"][key]] / max(compared[FIELDS["nanoaod_raw_sdst"][key]], 1)
                 for key in labels]
    ax.bar(np.arange(len(labels)), fractions, color="C0")
    ax.set_xticks(np.arange(len(labels)), FIELD_LABELS)
    ax.set_ylim(0, 1.3)
    finish(fig, ax, output / "source_comparison.png", sample, "Saved field", "Equal fraction")
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    numerator = np.asarray([gen_equal[a] for a, _ in gen_pairs])
    denominator = np.full(len(numerator), events["aligned"])
    lower, upper = clopper_pearson(numerator, denominator)
    ratio = numerator / denominator
    ax.bar(np.arange(len(numerator)), ratio, color="C0",
           yerr=np.maximum([ratio - lower, upper - ratio], 0), capsize=2)
    ax.set_xticks(np.arange(len(numerator)), ["PDG", "status", r"$p_x$", r"$p_y$", r"$p_z$", "E"])
    ax.set_ylim(0, 1.3)
    finish(fig, ax, output / "gen_alignment.png", sample, "Gen field", "Equal event fraction")
    print(f"{sample}: aligned source comparison {dict(events)}; Gen equal={dict(gen_equal)}, "
          f"RICH equal={dict(equal)}, compared={dict(compared)}")



def run(values, files, source, base, sample):
    output = plot_path(base, "03_rich", "00_data_validation", sample, source, SOURCES[source])
    storage(values, files, source, output / "branch_checks", sample)
    if source == "nanoaod_raw_sdst":
        source_comparison(files, base / "03_rich" / "00_data_validation" / sample / "source_comparison", sample)


if __name__ == "__main__":
    from .plot import main
    main("00_data_validation")
