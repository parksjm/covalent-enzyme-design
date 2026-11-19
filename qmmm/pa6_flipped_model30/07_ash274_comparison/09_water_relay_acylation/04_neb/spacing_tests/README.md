# Parallel NEB spacing and initialization comparison

All trials run on violet, under this existing NEB campaign. The original NEB
continues as historical context; it is not a matched restart control.

| Trial | Regular images | Seed controls | Optimizer | Purpose |
|---|---:|---|---|---|
| 01_uniform_24_lbfgs | 24 | Uniform QM arc length | L-BFGS | Matched numerical baseline |
| 02_uniform_48_lbfgs | 48 | Same uniform controls as 01 | L-BFGS | Image-density comparison |
| 03_adaptive_24_lbfgs | 24 | Double sampling weight in source images 9–10 | L-BFGS | Seed-detail comparison |
| 04_uniform_24_cg | 24 | Same uniform controls as 01 | CG | Short curvature-history diagnostic |

Inputs are atom-mapped, immutable reconstructions of original NEB cycle 268,
update_000075. Every trial retains DFTB3/3ob, electronic temperature 0 K,
QM charge −1, 102 real QM atoms/54 H/eight link atoms, ASP272/ASH274, neutral
Thr233 N-terminus and the same 4,662 movable atoms in the 12 Å region.
Frozen endpoints, exterior and force tolerances are unchanged. There are no
additional reaction-coordinate or positional restraints. Common `maxstep=0.02`
is in bohr per Cartesian component; it does not cap native kink repairs.
Native initialization still uses piecewise Cartesian interpolation. These tests
change its control-point distribution or image density, not its formula.

L-BFGS remains the primary optimizer. CG is not expected to be more efficient;
it tests whether the quasi-Newton history contributes to instability after
image resets. Its `maxcycle=250` produces the adapter's 1,750 energy/gradient
budget, roughly 79 full 24-image cycles. Budget exhaustion is not convergence.
The L-BFGS trials retain `maxcycle=16000`/112,000 evaluations.

## Explicit seed review

The 48-image seed contains one shared transferring proton (H20318, image 15),
1.521 Å from leaving N and 1.590 Å from OG1. Its original strict retention gate
remains false. `inputs/uniform/shared_proton_review.json` explicitly permits
that exact geometry as an unconverged path seed, pinning its manifest SHA.
All other 47 images pass the strict gate; all spectator-H, ASH274 and
compression checks pass. This does not certify a stationary state or impose a
restraint. Inter-image midpoint diagnostics are recorded separately from the
actual starting-image gates.

## Inspection and decisions

- Engine input: `<trial>/calculation/minimize.mdin`; full output stays on violet.
- Portable movie: `<trial>/progress/latest/vmd/view_path.vmd`
  when the first coherent export is available. Each trial has its own watcher.
- Energy/force comparison: `plots/comparison.png`, JSON and CSV; updated by
  `nyl12_neb_spacing_comparison_20261119`.
- Each trial has `launch.sh`, `job/`, `preflight.json` and retained source copies.
- The laptop pull includes current plots/configs/scripts and latest movies;
  it excludes raw calculation files and accumulated history. Refresh the
  laptop's pull script after these filter changes.

Compare force reduction, reset count/timing, image-spacing/chemical-coordinate
gaps and elapsed wall time at 50, 100 and 250 cycles where budgets allow.
Also use native `Calculating the energy` counts for work comparisons because
48-image cycles cost about twice as much as 24-image cycles. Lower profile
maxima alone are not a success criterion. Native convergence, chemical
connectivity and eventual TS curvature verification remain separate gates.

The first profiles are initialization diagnostics, not barriers: at early
cycles they contain maxima of approximately 65–82 kcal/mol. No improvement in
convergence or kinetic plausibility has yet been demonstrated.

Custom-code checks: 50 focused tests passed in the existing analysis environment;
Ruff, shell syntax and whitespace checks passed. The initial combined `uv`
collection lacked matplotlib and stopped; using the established analysis Python
with existing test tools corrected the environment without installing anything.
Measured new-helper/monitor/launcher coverage is 85–92%; the thin preparation
recipe measured 90% during actual accepted-input preparation.
