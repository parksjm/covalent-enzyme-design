# Updated reaction-distance diagnostics — 28 November 2026

MBAR remains the primary estimator. This report uses exactly the frozen
100-history / 82,000-frame pool in
`../review_20261128_01/checkpoint_001_main70ps_supplement43/`:
57 main histories at 70 ps, 33 corrective histories at 10 ps, four targeted
histories at 25 ps, and six initialization checks at 30 ps. Discard the first
5 ps once per continuous history. The exact retained order, all three CV
identities, and 922 source hashes were verified.

## Figures and interpretation

- `five_distances_vs_s.png/pdf`: complete path and all five requested distances.
- `five_distances_vs_s_early.png/pdf`: early-region enlargement.
- `five_distances_main57_70ps_independent_traces.png/pdf`: native trajectories;
  no lines join separate umbrella histories.
- `early_grid_key_windows_5to10ps.png/pdf`: relevant corrective histories.

Gray points show raw, biased observations. Red curves are MBAR-weighted
medians; orange envelopes show weighted empirical 10–90% distributions,
**not confidence intervals**. `distance_statistics.json`,
`bin_diagnostics.json`, and `early_bins_summary.csv` retain the numerical
results. `retained_rows.csv` preserves all reconstructed rows locally.

The observation contributing 96.84% of the weight at s=0.2625 has distances
OG1–HG1 0.9521, HG1–N(leaving) 4.1486, N(alpha)–HG1 1.7658,
OG1–C 3.3707, and C–N 1.3690 Å. These are reactant-like distances, not a
completed proton transfer or addition. Its displacement from the reference
path is 1.67 Å in the declared three-coordinate metric. The bin weight ESS
is 1.066; this Kish ESS measures weight concentration, not independently
sampled configurations after time correlation.

At s=0.2875, 90.33% of raw observations place HG1 nearest N(alpha), whereas
that population contributes only 0.138% of the reweighted probability.
Numerous biased frames therefore do not ensure reliable equilibrium support.
The early coordinate also includes nucleophile retreat without proton transfer.
These observations diagnose projection/coverage limitations; they do not
identify the physical origin of every PMF feature.

The existing 33-window grid ends at q_acyl=-1.75 Å, while the dominant frame
has q_acyl=-2.0018 Å. The approved boundary extension targets q_acyl=-2.0 Å
with q_act=-1.0, -0.8, -0.6 Å and two geometric starts per center. Matching
proton-state, five-distance, transverse-distance, weight-support and
initialization distributions are required before claiming improvement.

## Separate omission sensitivity

`../frame_omission_20261128/attempt_02/results/` contains a curiosity-only
full MBAR refit with this one accepted observation omitted. The affected
point rises from 2.89979 to 4.95910 kcal/mol; the profile maximum stays
26.2800 kcal/mol. Earlier features remain. It is explicitly **not for
publication** and supplies no scientific basis for excluding the frame.
All source observations and the primary analysis remain unchanged.

## Reproducibility and checks

`job_attempt02/command.json` and `run.sh` retain the analysis environment
and execution. `source_snapshot/` preserves the original source; the final
figure-only correction is in `figure_correction_source_snapshot/` with its
`local_plot_job/` receipt. Use PYTHONPATH pointing to the pinned project
source, OMP/BLAS threads=1, and the recorded existing uv environment.
The failed first import attempt is retained in `job/`.

Ruff and the weighted-quantile regression passed. New analysis-code coverage
was not measured. The separate omission core passed three tests with 93.94%
coverage. Job completion and these consistency checks do not establish
sampling convergence, a TI1 minimum, or a transition state.
