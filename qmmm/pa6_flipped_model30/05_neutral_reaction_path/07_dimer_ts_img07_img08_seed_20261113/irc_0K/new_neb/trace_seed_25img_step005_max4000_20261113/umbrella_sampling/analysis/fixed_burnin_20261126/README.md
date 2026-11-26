# Matched direct-route umbrella-sampling reanalysis — 2026-11-26

This is derived analysis of the existing direct-NEB-seeded campaign, not a new
simulation. All native trajectories, restarts, biases and failed attempts are
preserved. The relay US and its queued 60 ps continuation remain running.

## Data and statistical protocol

All 63 windows have six accepted 5 ps production segments (30 ps/window).
Heating (1 ps) and dedicated equilibration (5 ps) are separately excluded.
Concatenate continuous production, discard its first 5 ps once, then retain
all frames. Analyze both a matched 25 ps horizon (400 retained frames/window)
and the full 30 ps (500/window). Use the shared `fixed_us_statistics` module,
200 full MBAR refits with within-window contiguous 5 ps resampling, 2.5 ps
block sensitivity, 0/2.5/10 ps discard sensitivities and early/late halves.
Intervals are pointwise conditional precision, not barrier confidence intervals
or evidence of equilibration.

## Coordinate compatibility

The historical 47.05 kcal/mol direct result is a marginal along q_acyl. It must
not be compared directly with the relay 3-CV path histogram. For this revision,
recover r(Nalpha,HG1) from five-atom reads of every saved native frame; verify
all four previously measured distances and original times, source hashes,
process/acceptance records, q_H=q_act+q_del, and all archived cross-bias energies.
All 378 augmentations passed without censoring frames. Raw trajectories remain
on their owning hosts (96 segments Violet, 282 Azalea).

Both routes now use the same (q_acyl,q_act,q_del) Angstrom metric and extended-
tangent normalized-arclength projection. The direct physical reference is its
original 24 NEB images, verified against prepared seed hashes; the relay retains
its original reference. The numerical s values identify positions on each
route's own reference, not identical configurations. Preserve the relay 0.025
bin grid and s=0.0125 reference. One extra terminal bin is needed to retain all
direct fluctuations (raw range -0.055797 to 1.270101); no clipping is applied.
The grid extension changes no common interior bins or reference definition.

Retain actual native q_acyl/q_H rk=100 biases during reweighting; do not replace
them with relay activation/delivery rk=25 biases. The direct-seeded trajectories
can visit the N-terminal proton channel and are not chemically filtered into
a pure direct mechanism. Compare proton occupancy, overlap, projection ambiguity,
transverse excursions and temporal stability before mechanistic conclusions.

## Reproduction and files

- `coordinate_source/` and host plans: immutable sparse extraction code/inputs.
- `coordinate_job_{violet,azalea}/`: retained commands, hashes, process logs.
- `augmented/window_*/03_production_block_*/`: measured CSVs and source receipts.
- `analysis_source/` and its manifest: immutable analysis implementation.
- `fixed_{25,30}ps_job/`: tmux commands, environments, logs and status.
- `direct_{25,30}ps_analysis.json`: numerical PMFs, intervals and diagnostics.
- `physical_model_audit.json`: bounded physical-model identity comparison.
- Figures are in the existing campaign `plots/` directory.

No convergence or globally accepted activation barrier is asserted. The existing
relay PMF still has temporal drift and a weakly sampled acyl basin. Route-specific
projected densities are diagnostic comparisons, not a direct kinetic rate ratio.

## Completed results

| Campaign | Production/window | Retained frames | Sampled maximum, kcal/mol | Pointwise 95% interval at maximum |
|---|---:|---:|---:|---|
| Direct-NEB-seeded | 25 ps | 25,200 | 50.920 | Unresolved: 194/200 draws support the bin |
| Direct-NEB-seeded | 30 ps | 31,500 | 50.947 | Unresolved: 195/200 draws support the bin |
| N-terminal-relay-seeded | 25 ps | 22,800 | 24.807 | 23.938–25.553 |

These are sampled histogram maxima on route-specific coordinates, not accepted
activation barriers. The direct maximum at s=0.4125 has only 15 frames at both
horizons. The direct 25 ps peak changes from 49.229 to 52.670 kcal/mol across the
0/10 ps discard sensitivity. Early/late peaks are 50.333/48.540 kcal/mol. The relay
also drifts temporally (23.008/26.318 kcal/mol). All bootstrap MBAR networks are
connected, which does not establish sampling convergence.

HG1 is nearest to Nalpha in 2,763/25,200 direct retained frames (10.96%) and
3,463/31,500 (10.99%) at 30 ps. These are biased-window nearest-owner counts,
not equilibrium state populations or a strict covalent-bond assignment. They
confirm mixed channels in the direct-seeded data. Do not compute a pure
mechanism A/B rate ratio from this comparison.

The physical-model audit supports the same intended topology/QM partition and
frozen environment (maximum coordinate difference 6.32 × 10⁻⁶ Å), with different
SCC thresholds and different native bias schemes. Historical direct 3ob file
hashes were not pinned, so exact direct-era parameter identity is unverified.
Both analysis jobs completed with exit 0. Numerical summaries and source hashes
are in `comparison_summary.json`; use the matched 25/25 ps figure for comparison.
The full 30 ps direct plot is `plots/pmf_latest.png`. The historical
`analysis/pmf_latest/` remains the archived q_acyl marginal and is labeled as
historical by the pull script.
