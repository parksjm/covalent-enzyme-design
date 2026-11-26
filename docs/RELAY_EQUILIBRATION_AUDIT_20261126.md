# Audit of the unusually late relay-US truncation times

## Authorized replacement analysis

The PI subsequently authorized a new primary postprocessing protocol: one
fixed 5 ps additional production burn-in per continuous window, all later
frames, and 200 full MBAR refits using within-window contiguous 5 ps chunk
resampling. Native heating (1 ps) and equilibration (5 ps) were already excluded.
The extra discard is a transparent provisional choice informed by this audit,
not a claim of universal equilibration. Fixed 0/2.5/10 ps discards, temporal
halves and 2.5 ps bootstrap chunks quantify analysis sensitivity. Correlation
estimates remain diagnostics and no longer determine either discard or frame
selection. Insufficiently equilibrated windows remain scientific limitations;
they are not chemically censored or fixed by choosing a final quiet tail.

Keep all earlier analyses as provenance. New source-bound evidence is in
`full_path_300K/analysis/fixed_burnin_20261126/`, with figures in the existing
`full_path_300K/plots/` directory. The scripts reuse the same native acceptance
receipts, Hamiltonian, biases, projection, bins and E–S reference. See Methods
for R55–R60 and the distinction between stationary sampling uncertainty and
slow-state/preparation/model error. Native US and its 60 ps extension are
unchanged.

### Applicability to short QM/MM umbrella trajectories

The truncation and correlation statistics are not specific to classical force
fields or long trajectories. Their useful resolution depends on trajectory
length relative to the slow modes of the observable. A fast proton vibration
does not ensure that side-chain conformations, hydration or polymer orientation
have equilibrated. On short QM/MM records, uncertain correlation estimates and
too few transitions make automated truncation less reliable, not inapplicable.

For the current campaign, evaluate the fixed-discard PMF alongside the burn-in
sensitivities, temporal halves, block-size sensitivity, full-CV overlap,
proton ownership and slower environmental/geometry diagnostics. Extend the
existing windows to the authorized 60 ps. Persistent distribution gaps can
require bridge windows; temporal drift or hidden-state trapping instead calls
for equilibration, independent seeds or enhanced sampling. Neither a low
barrier nor small bootstrap errors is an acceptance criterion by itself.

For the final matched direct/relay comparison, prespecify the same protocol and
require a stable barrier region within a practical declared tolerance, stable
basin references and adequate overlap, together with independent-initialization
checks. A prospective target of about 1 kcal/mol statistical precision or
change between extended profiles is useful for planning, not a theorem that
the underlying PMF is equilibrated. Do not adapt cuts to match experimental
rates. Optimize sampling before interpreting a DFTB3/MM barrier quantitatively;
Hamiltonian uncertainty remains separate.

Recent primary examples are R59 (Jäckering et al., JCIM 2024, adaptive string,
explicit polymer-orientation coordinate, umbrella integration and confidence
criterion) and R60 (Zheng et al., Commun. Chem. 2026, 30 ps/window with last 15 ps
analyzed, convergence checks and committor analysis). R58 (Lu et al., 2016)
reviews the specific challenges of QM/MM free-energy sampling. Full verified
citations and links are in the manuscript bibliography. Their times and
acceptance choices are precedents rather than direct Nyl12 validation.

## Original audit

26 November 2026. Native simulations and analysis selections were not changed.
This audit examines the eight latest cutoffs in the complete 20 ps/window data.

## What the selector actually does

PyMBAR 4.0.3 `timeseries.detect_equilibration` tests candidate discard indices
and maximizes a correlation-based effective sample count for one observable.
In the installed implementation, its objective is `(N - t0 + 1) / g(t0)`;
`g` estimates statistical inefficiency from the suffix autocorrelation.
The default uses the approximate `fast=True` lag summation. The one-sample
offset is an implementation detail and does not explain the large discards.

Our `select_uncorrelated` function applies this procedure to three reaction
coordinates and ten distance/coordination diagnostics separately. It then
takes the **latest of the 13 optimized cutoffs**, recomputes the largest
statistical inefficiency over that common suffix using the full lag summation,
and subsamples at `ceil(g)`. This conservative aggregation is a project choice;
it is not a validated multivariate equilibration detector or a joint
optimization of the final common sample count.

The diagnostic objective shown in the plots concerns temporal correlation for
one observable. It is distinct from the bin-local MBAR weight ESS.

## Measured results

| Activation window | Applied cutoff, ps | Driving observable | Latest cutoff using full lag summation, ps | Biased-CV-only cutoff, ps |
|---|---:|---|---:|---:|
| 023 | 18.20 | q_del | 18.15 | 5.60 |
| 020 | 17.80 | HG1–leaving-N distance | 15.30 | 5.50 |
| 022 | 17.65 | Proton coordination progress | 15.90 | 6.20 |
| 006 | 17.40 | Proton coordination progress | 17.35 | 5.25 |
| 021 | 17.10 | Leaving-N proton coordination | 17.10 | 0.85 |
| 002 | 17.00 | Proton coordination progress | 6.85 | 1.70 |
| 024 | 15.80 | q_acyl | 17.20 | 15.80 |
| 011 | 15.25 | Proton coordination progress | 13.25 | 6.05 |

Biased-CV-only values are a diagnostic sensitivity check, not an adopted policy:
orthogonal coordinates remain relevant to whether the sampling is representative.
Both search variants use the same underlying time series. A difference between
them measures estimator sensitivity, not independent sampling uncertainty.

- **Window 023:** q_del changes substantially around 5–6 ps, then exhibits an
  excursion around 14.5–18 ps before returning toward its earlier distribution.
  The 18.2 ps choice isolates the final quieter suffix. The full-lag selector
  also chooses a late suffix, so approximate lag summation alone is not the
  cause. The trajectory supports an initial transient and slower fluctuations,
  but does not establish 18.2 ps as a unique physical equilibration time.
- **Window 002:** the four 5 ps means of the driving coordination metric are
  −0.0434, −0.0477, −0.0441 and −0.0461. There is no evident monotonic initial
  relaxation lasting 17 ps. The objective has competing local extrema; changing
  lag summation moves the cutoff from 17 to 6.85 ps. A 17 ps discard is not
  convincingly justified by these data.
- **Window 020:** the leaving-N–HG1 distance makes a temporary excursion around
  12–15 ps and returns. Keeping only the last 2.2 ps is not strong evidence of
  removing a sustained initial transient.
- **Windows 022, 006, 021 and 011:** small coordination changes and excursions
  drive the late choices. These are useful slow-mode diagnostics, but the
  short suffixes cannot certify adequate equilibrium sampling.
- **Window 024:** the biased acylation coordinate itself shifts around 15–16 ps.
  This is a stronger warning of nonstationarity; a universal early discard
  cannot be justified merely by rejecting the other extreme cutoffs.

Across all eight windows, HG1's nearest bonded heavy-atom owner remains
unchanged in all 400 frames. It remains on OG1 in windows 002/006/011 and on
the Thr N-terminal N in windows 020/021/022/023/024. The two original N-terminal
N–H distances remain below 1.18 Å. Thus the observed cutoff-driving changes
do not correspond to these protons switching bonded owners.

## Why the maximum-of-13 rule needs revision

After the final common-g thinning, a diagnostic 5 ps initial discard would
retain more configurations than the automatic cutoff in seven of these eight
windows. For window 023, the counts are 18 versus 6; for window 002, 18 versus 9
(31 if no production burn-in were removed). Therefore choosing the latest
individual optimum does not even guarantee maximizing the intended shared
effective sample count. Larger counts do not prove equilibration, and these
alternatives were **not used to replace the PMF**.

Treat these extreme cutoffs as uncertain analysis choices, not accepted
physical equilibration times. Before accepting a barrier, compare prespecified
initial discards, subsequent time-block stability, orthogonal populations and
the cutoff's stability as the trajectories lengthen. Do not choose a discard
because it produces a smoother or lower PMF. Keep the ongoing 60 ps extension.

The methodological concern is consistent with the distinction between
truncation selection and equilibration established by Chodera's original
method and emphasized by Clark, Cole and Michel's later benchmark. The latter
shows that noisy correlation estimates can produce spurious late optima and
that truncation heuristics do not reliably diagnose insufficient sampling.
No new detection package was installed or substituted in this audit.

## Artifacts and reproduction

Campaign root:
`qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/`.

- Data/source: `analysis/continuous_window_20261126/cutoff_audit_20261126/`.
  `extract_cutoff_audit.py`, `run_extraction.sh`, JSON and NPZ retain the
  hash-validated inputs, all 13 observables and both detector objective curves.
  The diagnostic objective argmax was checked against the installed PyMBAR API.
- Fixed-cutoff sensitivity: `assess_shared_cutoffs.py` and
  `shared_cutoff_sensitivity.json` in that directory.
- Figures: `plots/equilibration_worst_windows_20ps_20261126.png` and
  `plots/equilibration_activation_full_XXX_20ps_20261126.png`, with PDFs.
  Grey curves are raw frames; blue curves are 0.5 ps display averages; red
  dashed lines are applied cutoffs. Display averages never enter the detector.
- Forecast: `analysis/continuous_window_20261126/eta_25ps_20261126.json`.
  At 13:41 EST all 57 workers remained in stage07, steps 8899–9049. Recent rates
  predict complete 25 ps coverage around 14:08 EST; allow validation/analysis time.

## References

The manuscript Methods and Results now cite these as R55–R56, with
Grossfield et al.'s sampling-quality guide as R57. Complete entries are in
`manuscript/REFERENCES.md` and the companion BibTeX file. The Methods explicitly
separate the published single-observable criterion from our study-specific
latest-of-thirteen/common-maximum-g aggregation.

1. Chodera JD. A Simple Method for Automated Equilibration Detection in
   Molecular Simulations. *J Chem Theory Comput.* **2016**, 12, 1799–1805.
   [DOI:10.1021/acs.jctc.5b00784](https://doi.org/10.1021/acs.jctc.5b00784).
2. Clark F, Cole DJ, Michel J. Robust Automated Truncation Point Selection for
   Molecular Simulations. *J Chem Theory Comput.* **2025**, 21, 88–101
   (online 23 December 2024).
   [DOI:10.1021/acs.jctc.4c01359](https://doi.org/10.1021/acs.jctc.4c01359).
