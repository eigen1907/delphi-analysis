# Photon validation

Run `bash scripts/03_sample_validation/02_gamma/plot.sh` from the repository root.
Defaults: `output/dataset/20261015_jongwon`, `train`, all stages. Use `--input`,
`--split` and `--stage` to change them. The implementation is in
[data.py](../src/sample_validation/gamma/data.py),
[angular_matching.py](../src/sample_validation/gamma/angular_matching.py) and
[detector_plots.py](../src/sample_validation/gamma/detector_plots.py).

## Gen selection and acceptance

Start with stable photons: `GenPart_pdgId == 22` and `GenPart_status == 1`.
Ancestry separates beam ISR, other ISR, FSR and the remaining `decayed` category.
Beam ISR has `px == py == 0` and an incoming-electron parent.

| Matching population | Selection |
| --- | --- |
| `gen_gamma` | All photons except beam ISR, within EM acceptance |
| `gen_isr` | Non-beam ISR within EM acceptance |
| `gen_no_isr` | FSR or `decayed`, within EM acceptance |

The fiducial regions use the **Gen photon polar angle in degrees**. All bounds
are open; boundaries and gaps are excluded.

| Detector | Accepted θ |
| --- | --- |
| HPC | `44 < θ < 88` or `92 < θ < 136` |
| FEMC | `11 < θ < 36` or `144 < θ < 169` |
| STIC | `2 < θ < 10` or `170 < θ < 178` |

EM acceptance is their union. **All ranges (0–180°)** combines accepted photons;
it is not a 4π denominator. Reco candidates receive no preliminary θ cut.
Gen four-vectors must be finite, with positive energy and nonzero momentum.

## Reco collections and matching

| Collection | Meaning |
| --- | --- |
| Native `Photon` | Stored neutral EM candidates; not a confirmed photon-PID sample |
| Native `Photon_PhotonConv` | Photon plus conversions, deduplicated through stored parent/daughter relations |
| Native `Part_pdgId_21` | DELPHI mass code 21; may include unidentified neutral particles |
| Raw `EmShower`, `Stic`, `EmShower_Stic` | Clusters with finite positive energy and finite nonzero direction; charged-PA clusters are retained |

Raw EM directions use `cluster position − simulated PV`; STIC uses the stored
MAIN momentum. Raw collections have no Reco→Gen Sim links, so only angular
matching is available. Their results measure **angular EM coverage**, not photon
identification efficiency.

- **Angular:** independent 3D opening-angle cones around each Gen photon. One
  Reco may match several Gen photons; this is not one-to-one matching.
- **Truth link:** use a unique Gen origin from the stored Native relations in
  both directions. Stop Sim ancestry at the first Gen anchor. Missing or
  conflicting origins are not replaced by angular matches.
- **Default cone: 1°**, identical for both campaigns and all Native/Raw
  collections. The `0.25, 0.5, 1, 2, 3, 5, 10, 20°` scan is supplementary.

The 1° cut follows §3.4 of the
[DELPHI μμγ paper](https://pure.uva.nl/ws/files/2867610/1024_11797y.pdf), where it
selects an energy-response sample. Applying it to STIC is an extension; that
paper does not include STIC. Truth efficiency also depends on link completeness:
missing/conflicting links and overlapping photons can lower it independently
of reconstruction efficiency.

## Detector efficiency

For each energy band, population and Gen region, let `G` be the selected Gen
photons. Both plots keep unmatched photons in the denominator and omit empty
denominators.

| Plot | Denominator | Numerator |
| --- | --- | --- |
| `eff_det_regional` | All photons in `G` | Gen with at least one matched Reco |
| `eff_det_flag` | The same `G` | Gen with a matched Reco carrying the required detector tag |

Detector bins use Gen θ. Regional efficiency does not require a Reco detector
tag; the All flag accepts any HPC/FEMC/STIC tag. Native uses the stored
`Photon_isHPC/isFEMC/isSTIC`; conversions receive no invented detector tag. Raw
uses EM detector ID 9 for HPC, 26 for FEMC, and the STIC collection. A conversion
alone can therefore pass regional efficiency but fail flag efficiency.

## Energy bands and plots

Bands use **Gen energy in GeV**, with the same selection in numerator and
denominator. They do not restrict the Reco pool.

| Band | Gen energy | Filename suffix |
| --- | --- | --- |
| Default | `0 ≤ E`, no upper bound | None |
| Low | `0 ≤ E < 1` | `_energy0-1` |
| Mid | `1 ≤ E < 10` | `_energy1-10` |
| High | `10 ≤ E` | `_energy10-` |

The bands apply to `eff_det_flag`, `eff_det_regional`, `eff_vs_cos_theta` and
`event_energy_ratio`. Energy sums retain unmatched events as zero; ambiguous
overlapping Reco energies remain NaN. Angular event sums count a Reco only once,
even if it matches several Gen cones.

`01_sample_distribution` includes beam ISR and photons outside acceptance, so
its Gen counts differ from matching denominators. Energy plots are:

| File | Axes |
| --- | --- |
| `energy.png` | Linear x and y |
| `energy_log.png` | Linear x, logarithmic y |
| `energy_loglog.png` | Logarithmic bins and y; symlog x retains zero energy |

In combined Reco distributions, PhotonConv is stacked **below** Photon after
deduplication. Both use the same color and opposite diagonal hatches; their
total equals the combined candidate count.

Matching output:
`output/plots/<campaign>/03_sample_validation/<split>/02_gamma/<stage>/matching_result/<sample>/<source>/<collection>/`.
