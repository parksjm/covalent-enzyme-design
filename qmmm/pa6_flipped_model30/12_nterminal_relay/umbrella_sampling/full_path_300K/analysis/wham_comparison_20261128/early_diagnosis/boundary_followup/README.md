# Conditional-boundary sensitivity and dominant-frame receipt

Separate followup; the completed parent report and figures are unchanged.
Same immutable 82,000-frame pool, MBAR weights, physical q_act bins and
reference cell. No statistical refit, new trajectories or simulation changes.

| q_acyl lower bound (upper remains −0.8 Å) | Selected frames | Monotonic q_act −1.05→+0.25 Å? | PMF at q_act=−0.85 | ESS there | Dominant hard-bin frame's fraction there |
|---:|---:|:---:|---:|---:|---:|
| −2.00 Å | 34,978 | Yes | 0.91992 | 97.64 | Excluded |
| −2.05 Å | 34,985 | Yes | 0.89382 | 89.15 | 4.283% |
| −2.10 Å | 34,987 | Yes | 0.89382 | 89.15 | 4.283% |

The boundary concern is real: the frame that dominates the hard s≈0.2625
bin lies only 0.001763 Å outside the original q_acyl≥−2 restriction.
Both wider restrictions include it. Across the stated q_act interval the
largest absolute PMF change is 0.41378 kcal/mol; no downward steps appear.
The endpoint at q_act=+0.25 remains 13.16029 kcal/mol. Thus the different
appearance of this physical marginal is not produced solely by censoring
that frame at q_acyl=−2. Nevertheless, this is still a different coordinate
and conditional ensemble; it cannot establish that every path feature is
an artifact, or certify any physical basin or barrier.

## Native origin and bias of the high-weight frame

- State `activation_grid33/activation_support_3_03`; pooled zero-based index
  78,151. Stage `04_production_5ps`, raw CSV zero-based row51 / native frame52.
- Exact retained q matches the native raw-feature row to 1e−14 Å.
  Native trajectory time12.6000003815 ps corresponds to production age7.60 ps
  after the 5 ps initialization offset, 2.60 ps past the once/history discard.
- Restrained coordinates: q_acyl and q_act; centers(−1.75,−0.60) Å and
  native Amber force parameters(25,25) kcal mol⁻¹ Å⁻², outer extent10 Å.
  Displacements(−0.25176264,−0.21377617) Å remain in the harmonic region.
- Existing implementation uses **U=ΣkΔq²**, without a factor of 1/2.
  Contributions1.58461073+1.14250630 give **2.72711702 kcal/mol =4.57446 kT**.
  Independently recorded native bias2.73250 differs by0.005383 kcal/mol,
  within the row's0.012394 kcal/mol serialization/quantization bound.
- Native process exited0. Stage acceptance and analysis acceptance are true;
  every integration step SCC-converged, with no recovered or unresolved trials.
  This rules out a known SCC failure for this accepted stage, not statistical
  nonstationarity or chemical shortcomings.

The source window does not restrain q_del. This frame's q_del=−2.382732 Å is
the second-lowest value among its100 retained source frames, whose range is
[−2.398433,−1.439894] Å. Its total source bias is moderate rather than an
enormous restraint excursion. A neighboring existing grid state3_02 has an
even lower bias1.58936 kcal/mol at this configuration. The fitted MBAR
denominator receives56.25%,25.33% and16.28% from grid3_03,3_02 and3_04,
respectively. These are algebraic denominator responsibilities, **not observed
independent visits or temporal overlap evidence**.

This distinguishes two problems: the configuration lies in reach of several
biases, while its physical probability contribution to the selected hard
histogram bin is represented by effectively one retained frame. Presence of
nearby umbrella centers does not certify coverage of the relevant joint
q_del/transverse tail. Its actual importance-tail persistence must be assessed
with additional accepted sampling; this followup launches none.

## Distance from the reference polyline

| Hard s bin center | Weighted mean transverse distance (Å) | Largest-weight frame's distance (Å) |
|---:|---:|---:|
| 0.1125 | 0.6182 | 0.7238 |
| 0.1625 | 1.0514 | 1.1889 |
| 0.2375 | 1.4864 | 1.5569 |
| 0.2625 | 1.6607 | 1.6702 |

Distances use the original equal-Angstrom three-CV metric. These weighted
histogram bins receive substantial probability away from the reference
polyline, making path progress an imperfect label for reference-like local
chemistry. This observation does not by itself invalidate such configurations.

`sensitivity_summary.json`, `physical_profiles.json` and its CSV,
`dominant_frame_receipt.json`, and `transverse_support.json` retain numerical
evidence. Native acceptance and analysis receipts are embedded in the frame
receipt. `input_sha256.json`, source snapshots and `job/` retain provenance.
Ruff passed; original physical profile reproduction and exact native-row
matching assertions passed; one-thread nice19 job exited0. No coverage
measurement or sampling-convergence claim.
