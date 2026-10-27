"""Two-track control selections with channel-specific detector requirements."""
import numpy as np

from sample_preparation.run_quality import TRACKING, QUALITY_MIN, QUALITY_MAX


DETECTORS = TRACKING + ("mub_b", "mub_d", "muf_a", "muf_c",
                        "hpc_0", "hpc_1", "emf_a", "emf_c")
CHANNELS = {
    "mumu": dict(mass=0.1056583755, label=r"$Z\to\mu^+\mu^-$",
                 mass_label=r"$m_{\mu\mu}$", pidkey="muid"),
    "ee": dict(mass=0.00051099895, label=r"$Z\to e^+e^-$",
               mass_label=r"$m_{ee}$", pidkey="elid"),
}
STAGES = {
    "mumu": ("two_os", "back_to_back", "fiducial", "standard_id", "detector_quality", "mass_window"),
    "ee": ("two_os", "back_to_back", "fiducial", "standard_id", "detector_quality", "mass_window"),
}


def standard_id(tags, channel):
    pid = tags[CHANNELS[channel]["pidkey"]]
    return ((pid >= 0) & ((pid & 4) != 0)) if channel == "mumu" else pid >= 4


def select_pairs(momentum, mass, tags, channel, flags=None):
    """Start from valid two-OS Part pairs; quality applies to Data only.

    mass uses the channel's muon or electron mass hypothesis.
    Both channels include their angular boundaries. Muon intervals follow the
    published muon fiducial; electron intervals are a conservative working choice.
    A is negative z; C is positive z. HPC 0/1 are both required, without
    assuming that their numbering identifies a z side.
    """
    theta = np.degrees(np.arctan2(np.hypot(momentum[:, :, 0], momentum[:, :, 1]),
                                  momentum[:, :, 2]))
    if channel == "mumu":
        barrel = ((53 <= theta) & (theta <= 88.5)) | ((91.5 <= theta) & (theta <= 127))
        forward_c = (20 <= theta) & (theta <= 42)
        forward_a = (138 <= theta) & (theta <= 160)
    else:
        barrel = ((44 <= theta) & (theta <= 88)) | ((92 <= theta) & (theta <= 136))
        forward_c = (20 <= theta) & (theta <= 35)
        forward_a = (145 <= theta) & (theta <= 160)
    fiducial = (barrel | forward_a | forward_c).all(axis=1)
    category = np.full(len(mass), "Other", dtype="U5")
    category[barrel.all(axis=1)] = "BB"
    category[(forward_a | forward_c).all(axis=1)] = "FF"
    category[fiducial & (category == "Other")] = "BF"

    tracking_quality = np.ones(len(mass), dtype=bool)
    detector_quality = np.ones(len(mass), dtype=bool)
    if flags is not None:
        good = (flags >= QUALITY_MIN) & (flags <= QUALITY_MAX)
        tracking_quality = good[:, :4].all(axis=1)
        if channel == "mumu":
            barrel_good = good[:, 4:6].all(axis=1)
            a_good, c_good = good[:, 6], good[:, 7]
        else:
            barrel_good = good[:, 8:10].all(axis=1)
            a_good, c_good = good[:, 10], good[:, 11]
        detector_quality = ((~barrel.any(axis=1) | barrel_good)
                            & (~forward_a.any(axis=1) | a_good)
                            & (~forward_c.any(axis=1) | c_good))

    identified = standard_id(tags, channel)
    phi = np.arctan2(momentum[:, :, 1], momentum[:, :, 0])
    difference = phi[:, 0] - phi[:, 1]
    delta_phi = np.abs(np.arctan2(np.sin(difference), np.cos(difference)))
    masks = {"two_os": tracking_quality}
    masks["back_to_back"] = masks["two_os"] & (np.pi - delta_phi < 0.35)
    masks["fiducial"] = masks["back_to_back"] & fiducial
    masks["standard_id"] = masks["fiducial"] & identified.all(axis=1)
    masks["detector_quality"] = masks["standard_id"] & detector_quality
    masks["mass_window"] = masks["detector_quality"] & (mass > 70) & (mass < 110)
    return category, masks
