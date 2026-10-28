"""Native high-level inputs in the final two-lepton selections.

Each row is (branch, common bin edges, axis label, interpretation note).
The 48 distributions are the native fields in docs/branches/classification.md,
with the stored four-momentum split into its four components. Particle and PID
rows refer to the two selected Part legs; event scalars enter once, and vertex
rows include every native vertex in a selected event. Fixed bins are shared by
Data and MC in both channels. Finite negative values, sentinels and tails are
retained by the reader and recorded in the flow table.

The four external features below have no established local definition. They
are deliberately absent from the histogram specifications.
"""
import numpy as np


UNDEFINED_FEATURES = ("n_extra", "photon_E_sum", "max_Muid_tag", "n_muid")


HIGH_LEVEL_VARIABLES = {
    "01_particles": [
        *[(f"Part_fourMomentum.fCoordinates.f{axis}", np.linspace(-70, 70, 141),
           rf"Selected particle $p_{axis.lower()}$ [GeV]", "Stored native momentum component")
          for axis in "XYZ"],
        ("Part_fourMomentum.fCoordinates.fT", np.linspace(0, 70, 141),
         "Selected particle energy [GeV]", "Stored native energy; finite tails retained"),
        ("Part_charge", np.arange(-1.5, 2.5), "Selected particle charge",
         "Stored int8 preserves -1/0/+1, but cannot preserve the source unknown code 999"),
    ],
    "02_calorimeter": [
        ("Part_hpcTotalShowerEnergy", np.linspace(0, 70, 141),
         "Associated HPC energy [GeV]", "Sum of associated HPC shower energies"),
        ("Part_hpcNumLayers", np.arange(-0.5, 10.5), "Associated HPC layers",
         "Layer count of the strongest associated HPC shower"),
        ("Part_hacTotalShowerEnergy", np.linspace(0, 70, 141),
         "Associated HAC energy [GeV]", "Sum of associated HAC shower energies"),
        ("Part_hacNumTowers", np.arange(-0.5, 61.5), "Associated HAC towers", ""),
        ("Part_sticShowerEnergy", np.linspace(0, 70, 141),
         "Associated STIC energy [GeV]",
         "Original Florian native STIC COMMON stride was wrong; stored values shown without correction"),
        ("Part_sticNumTowers", np.arange(-0.5, 25.5), "Associated STIC towers",
         "Original Florian native STIC COMMON stride was wrong; extreme stored values enter flow table"),
    ],
    "03_tracking": [
        ("Trac_length", np.linspace(-400, 400, 161), "Signed native track length [cm]",
         "Sign encodes primary-vertex fit participation; stored sign retained"),
        ("Trac_impParRPhi", np.linspace(-1, 1, 101),
         r"Native track impact parameter $R\phi$ [cm]",
         "Origin-relative scalar; Data/MC beam references differ"),
        ("Trac_impParZ", np.linspace(-10, 10, 101),
         r"Native track impact parameter $z$ [cm]", "Origin-relative scalar"),
        ("Trac_chi2VD", np.linspace(0, 160, 81), r"Native track VD $\chi^2$", ""),
        ("Trac_ndfVD", np.arange(-0.5, 106.5), "Native track VD degrees of freedom", ""),
        ("Trac_chi2NoVD", np.linspace(-20, 20, 101), r"Stored native non-VD $\chi^2$",
         "Negative/extreme source values retained; conventional chi-squared meaning is not established"),
        ("Trac_firstPointR", np.linspace(0, 20, 101), "Native first track point radius [cm]",
         "Central range shown; finite tails retained in flow table"),
        ("Trac_firstPointZ", np.linspace(-50, 50, 101), "Native first track point z [cm]",
         "Central range shown; finite tails retained in flow table"),
    ],
    "04_pid": [
        ("Dedx_value", np.linspace(0, 4, 81), "Native corrected dE/dx [MIP units]",
         "Native corrected response; negative stored sentinels retained in flow table"),
        ("Muid_tag", np.arange(-0.5, 32.5), "Native muon tag",
         "Bitmask; dimuon sample is conditioned on both legs passing native muon ID"),
        ("Muid_looseChi2", np.linspace(0, 20, 101), r"Native muon loose $\chi^2$",
         "Conditional selected-leg shape; finite tails retained in flow table"),
        ("Elid_tag", np.arange(-0.5, 6.5), "Native electron tag",
         "Dielectron sample is conditioned on both legs passing native electron ID"),
        ("Elid_gammaConversion", np.arange(-0.5, 8.5), "Native electron conversion tag", ""),
        *[(f"Haid_{name}", np.arange(-1.5, 4.5), label,
           "Stored hadron-hypothesis PID codes, including -1; lepton controls do not measure hadron-ID efficiency")
          for name, label in (("kaonDedx", "Native kaon dE/dx tag"),
                              ("protonDedx", "Native proton dE/dx tag"),
                              ("kaonRich", "Native kaon RICH tag"),
                              ("protonRich", "Native proton RICH tag"),
                              ("pionRich", "Native pion RICH tag"))],
        ("Haid_richQuality", np.arange(-1.5, 513.5), "Native RICH quality word",
         "Stored int32 word; large codes enter flow table; old Rich layout defects do not invalidate every Haid tag"),
        ("Rich_theg", np.linspace(0, 0.4, 81), "Stored native gas RICH angle [rad]",
         "Original Florian native gas RICH COMMON layout was wrong; stored values shown without correction"),
        ("Rich_thel", np.linspace(0, 0.8, 81), "Stored native liquid RICH angle [rad]",
         "Original Florian native liquid RICH COMMON layout was wrong; stored values shown without correction"),
        *[(f"Haidc_{name}Tag", np.arange(-1.5, 4.5), f"Native combined {name} tag",
           "Stored KHAIDC PID codes, including -1; lepton controls do not measure hadron-ID efficiency")
          for name in ("pion", "kaon", "proton", "heavy")],
    ],
    "05_event": [
        ("nPart", np.arange(-0.5, 21.5), "Native particles per event", "Full native Part collection count"),
        ("nPhoton", np.arange(-0.5, 16.5), "Native photon candidates per event",
         "Neutral EM candidates; excludes PhotonConv; original STIC layout defect can affect candidate membership"),
        ("nV0", np.arange(-0.5, 6.5), "Native V0 candidates per event", ""),
        ("Event_chargedMult", np.arange(-0.5, 11.5), "Native charged multiplicity",
         "Original PA charged multiplicity"),
        ("Event_neutralMult", np.arange(-0.5, 21.5), "Native neutral multiplicity",
         "Original PA neutral multiplicity"),
        ("Event_totalChargedEnergy", np.linspace(0, 140, 141),
         "Native selected charged energy [GeV]", "Charged energy passing the Team4 selection"),
        ("Event_totalEMEnergy", np.linspace(0, 100, 101),
         "Native neutral EM energy [GeV]", "Neutral EM energy; not the sum over charged Part showers"),
        ("Event_totalHadronicEnergy", np.linspace(0, 70, 71),
         "Native neutral hadronic energy [GeV]", "Neutral hadronic energy"),
        ("Event_hadronTagT4", np.arange(-0.5, 2.5), "Native Team4 hadron tag", "Stored boolean"),
    ],
    "06_vertex": [
        ("nVtx", np.arange(-0.5, 11.5), "Native vertices per event", "Final SKELANA vertex count"),
        ("Vtx_chi2", np.linspace(-1, 100, 102), r"Stored native vertex $\chi^2$",
         "All native vertices in selected events; negative values retained; final bank differs from Raw PV"),
        ("Vtx_nOutgoing", np.arange(-0.5, 16.5), "Native vertex outgoing count",
         "All native vertices in selected events; final bank differs from Raw PV"),
    ],
}
