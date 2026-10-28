"""Raw Transformer inputs, grouped for the final-selection Data/MC plots.

Each row is (branch, common bin edges, axis label, interpretation note).
The 57 fields in docs/branches/classification.md become 65 distributions when
the four XYZ vectors are split into their stored components. Bins show the
selected samples' central ranges; the reader records every finite tail in its
flow table. No negative values or stored sentinels are removed.
"""
import numpy as np


VARIABLES = {
    "01_track": [
        ("TracRaw_theta", np.linspace(0, np.pi, 65), r"Track $\theta$ [rad]", ""),
        ("TracRaw_phi", np.linspace(0, 2 * np.pi, 65), r"Track $\phi$ [rad]", ""),
        ("TracRaw_invR", np.linspace(-0.001, 0.001, 101), r"Track inverse radius [cm$^{-1}$]", ""),
        ("TracRaw_charge", np.arange(-1.5, 2.5), "Track charge", ""),
        ("Event_bFieldGevCm", np.linspace(0.003689, 0.003693, 41), "Magnetic-field conversion [GeV/cm]", ""),
        ("TracRaw_impactRPhi", np.linspace(-1, 1, 101), r"Track impact parameter $R\phi$ [cm]", ""),
        ("TracRaw_impactZ", np.linspace(-10, 10, 101), r"Track impact parameter $z$ [cm]", ""),
        ("TracRaw_trackLength", np.linspace(0, 400, 81), "Track length [cm]", "Absolute stored length"),
        ("TracRaw_chi2VD", np.linspace(0, 120, 61), r"Track VD $\chi^2$", ""),
        ("TracRaw_ndfVD", np.arange(-0.5, 76.5), "Track VD degrees of freedom", ""),
        ("TracRaw_chi2NoVD", np.linspace(-20, 20, 101), r"Track non-VD $\chi^2$", "Negative/extreme source values retained; tails in flow table"),
        ("TracRaw_firstPointR", np.linspace(0, 200, 81), "First track point radius [cm]", ""),
        ("TracRaw_firstPointZ", np.linspace(-250, 250, 101), "First track point z [cm]", ""),
    ],
    "02_tpc": [
        ("MtpcRaw_dEdx80Max", np.linspace(0, 4, 81), "TPC dE/dx 80 maximum [stored units]", ""),
        ("MtpcRaw_dEdx80Sigma", np.linspace(0, 2, 81), "TPC dE/dx 80 width [stored units]", "Distribution width, not measurement error"),
        ("MtpcRaw_dEdx80Integrated", np.linspace(0, 4, 81), "TPC integrated dE/dx 80 [stored units]", ""),
        ("MtpcRaw_dEdx65Max", np.linspace(0, 4, 81), "TPC dE/dx 65 maximum [stored units]", ""),
        ("MtpcRaw_dEdx65Sigma", np.linspace(0, 2, 81), "TPC dE/dx 65 width [stored units]", "Distribution width, not measurement error"),
        ("MtpcRaw_zFitChi2", np.linspace(-10, 220, 93), r"TPC z-fit $\chi^2$", "Negative source values retained"),
    ],
    "03_em": [
        ("EmShower_energy", np.linspace(0, 70, 71), "EM shower energy [GeV]", ""),
        *[(f"EmShower_position.fCoordinates.f{axis}", np.linspace(-350, 350, 101),
           f"EM shower position {axis.lower()} [cm]", "") for axis in "XYZ"],
        ("EmShower_nLayers", np.arange(-0.5, 11.5), "HPC shower layers", "HPC only (detector 9); FEMC has no layer information"),
        ("EmShower_detector", np.arange(8.5, 27.5), "EM detector code", "HPC = 9; FEMC = 26"),
    ],
    "04_hadronic": [
        ("HadShower_energy", np.linspace(0, 60, 61), "Hadronic shower energy [GeV]", ""),
        *[(f"HadShower_direction.fCoordinates.f{axis}", np.linspace(-500, 500, 101),
           f"Hadronic shower position {axis.lower()} [cm]", "Stored direction is an XYZ position") for axis in "XYZ"],
        ("HadShower_nHits", np.arange(-0.5, 5.5), "Hadronic shower layer rows", "Counts layer rows, not total towers"),
    ],
    "05_stic": [
        ("Stic_energyFromMain", np.linspace(0, 60, 61), "STIC MAIN energy [GeV]", ""),
        *[(f"Stic_directionFromMain.fCoordinates.f{axis}",
           np.linspace(-60, 60, 121) if axis == "Z" else np.linspace(-10, 10, 81),
           f"STIC MAIN momentum {axis.lower()} [GeV]", "Stored direction is momentum") for axis in "XYZ"],
        ("Stic_numHitTowers", np.arange(-0.5, 25.5), "STIC hit towers", ""),
    ],
    "06_muon_id": [
        ("MuidRaw_tag", np.arange(-0.5, 32.5), "Stored raw muon tag", "Storage diagnostic: module number 24 replaces the muon tag"),
        ("MuidRaw_looseChi2", np.arange(-0.5, 128.5), r"Stored raw muon loose $\chi^2$", "Storage diagnostic: block bit pattern replaces chi-squared"),
        ("MuidRaw_hitPattern", np.arange(-0.5, 544.5, 8), "Stored raw muon hit pattern", "Storage diagnostic: wrong pattern before inefficiency correction"),
    ],
    "07_electron_id": [
        ("ElidRaw_tag", np.arange(-0.5, 32.5), "Stored raw electron tag", "Storage diagnostic: baseline files contain ELID header value 25"),
        ("ElidRaw_gammaConvTag", np.arange(9035.5, 9045.5), "Stored raw conversion tag", "Storage diagnostic: baseline files contain unrelated value 9040"),
    ],
    "08_hadron_id": [
        *[(f"HaidRaw_{name}", np.arange(-1.5, 4.5), label, "Stored PID codes, including -1")
          for name, label in (("kaonDedx", "Kaon dE/dx tag"),
                              ("protonDedx", "Proton dE/dx tag"),
                              ("kaonRich", "Kaon RICH tag"),
                              ("protonRich", "Proton RICH tag"),
                              ("pionRich", "Pion RICH tag"))],
        ("HaidRaw_kaonCombined", np.arange(-1.25, 4.76, 0.5), "Combined kaon PID", "Stored values, including -1 and half-integer codes"),
        ("HaidRaw_protonCombined", np.arange(-1.25, 4.76, 0.5), "Combined proton PID", "Stored values, including -1 and half-integer codes"),
        ("HaidRaw_richQuality", np.arange(-128.5, 128.5), "Stored RICH quality word", "Storage diagnostic: baseline status word was narrowed to int8"),
        ("HaidRaw_thetaGas", np.linspace(0, 0.4, 81), "Stored gas RICH angle [rad]", "Storage diagnostic: baseline gas bank has incorrect layout"),
        ("HaidRaw_thetaLiq", np.linspace(0, 0.8, 81), "Stored liquid RICH angle [rad]", "Storage diagnostic: baseline liquid bank has incorrect layout"),
    ],
    "09_vertex": [
        *[(f"Vtx_position.fCoordinates.f{axis}",
           np.linspace(-10, 10, 101) if axis == "Z" else np.linspace(-1, 1, 101),
           f"Raw vertex position {axis.lower()} [cm]", "Raw PV bank; tails in flow table") for axis in "XYZ"],
        ("Vtx_chi2", np.linspace(0, 100, 101), r"Raw vertex $\chi^2$", "Raw PV bank; tails in flow table"),
        ("Vtx_nOutgoing", np.arange(-0.5, 16.5), "Raw vertex outgoing count", "Raw PV bank"),
    ],
    "10_multiplicity": [
        ("nTracRaw", np.arange(-0.5, 11.5), "Raw tracks per event", ""),
        ("nEmShower", np.arange(-0.5, 16.5), "EM showers per event", ""),
        ("nHadShower", np.arange(-0.5, 11.5), "Hadronic showers per event", ""),
        ("nMtpcRaw", np.arange(-0.5, 11.5), "TPC module rows per event", ""),
        ("nVdAssocHit", np.arange(-0.5, 41.5), "Associated VD hits per event", ""),
        ("nVtx", np.arange(-0.5, 11.5), "Raw vertices per event", ""),
        ("nStic", np.arange(-0.5, 11.5), "STIC rows per event", ""),
    ],
    "11_links": [
        ("EmLayer_emShowerIdx", np.arange(-0.5, 16.5), "EM layer to shower row index", "Zero-based row index, not a physics observable"),
        ("HadHit_hadShowerIdx", np.arange(-0.5, 11.5), "Hadronic hit to shower row index", "Zero-based row index, not a physics observable"),
        ("VdAssocHit_tracRawIdx", np.arange(-0.5, 11.5), "VD hit to raw track row index", "Zero-based row index, not a physics observable"),
    ],
}
