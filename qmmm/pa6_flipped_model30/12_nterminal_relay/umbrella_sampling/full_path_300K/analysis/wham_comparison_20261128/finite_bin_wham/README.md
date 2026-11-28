# Finite-bin WHAM sensitivity diagnostic, 28 November 2026

This derived analysis compares the frozen 82,000-frame / 100-history pooled
MBAR estimate with finite-bin WHAM, using the existing PyMBAR 4.0.3 numerical
solver. It is an estimator-discretization diagnostic; it is not an independent
software cross-check, new sampling, or evidence of a converged activation barrier.

`provenance.json` records the hypothesis, exact input/source hashes, installed
versions and numerical acceptance gates. All three input artifacts and relevant
native Python modules are retained in `frozen_inputs/` and `source_snapshot/`.
The source snapshot is the executed implementation. No underlying molecular
Hamiltonian, protonation, atom mapping, restraint or sampling job was changed.

The 3D physical-CV grid uses q_acyl, q_act and q_del, with widths 0.1, 0.05,
0.025 and 0.0125 Å and edges anchored at zero. The finest grid is also shifted
by half a cell along every axis. Exact native Amber bias functions, including
linear tails, are evaluated at occupied cell centers. All 100 original retained
counts are supplied to the robust solver, warm-started from the cached MBAR
free energies. Original-frame MBAR weights are reused after explicit parity
checks; no baseline refit is performed.

Finite-bin WHAM equations are solved by expanding each center-bias column for
each original frame assigned to that cell. Each frame within a cell receives
its cell probability divided by its original occupancy. Primary plots retain
the original frames' hard-s values and the frozen final bins/reference. The
separately labeled center-projection plot bins cell masses at projected cell
centers; this adds projection discretization and is not the primary comparison.
Empty bins remain unsupported and plotted curves are unsmoothed.

`attempt_01/` retains the persistent job command, environment overrides,
combined log and process `status.json`; `results/scientific_acceptance.json`
records separate numerical/scientific scope. The launch uses nice 19, one
numerical CPU thread and no GPU in `relay_finite_bin_wham_20261128`.

Qualification: Ruff formatting/lint passed; 10 analytic pytest tests passed,
covering one-state correction, known multistate masses, state-gauge invariance,
identical-bias splitting, within-cell equality, projection/mass equivalence,
fixed reference and invalid-input rejection. Core unit-test coverage is 98%;
runner/figure code coverage was not measured. The original pytest-cov command
also wrote the repository's ignored `.coverage`; this bookkeeping side effect
was disclosed, and its data were copied to `qualification.coverage`.

Interpretation belongs with the lead campaign review. In particular, a smoother
coarse-grid curve does not establish greater accuracy or sampling convergence.
The 0.1 kcal/mol finest-grid/shift threshold was declared before the run as a
numerical diagnostic only; exceeding it is reported without protocol tuning.

## Completed result

The tmux attempt completed with exit 0 in 2 minutes 46 seconds. All frozen
source/input hashes, baseline parity, normalization and fixed-point checks
passed. All five grids retain the projected maximum at s = 0.6875. MBAR gives
26.28003 kcal/mol; the 0.0125 Å WHAM result is 26.28371, or 26.27029 with the
half-cell shift. These are projected maxima, not accepted activation barriers.

For 0 ≤ s ≤ 0.3, the finest grid has RMSD 0.02918 kcal/mol and maximum absolute
deviation 0.09941 from MBAR. The half-shift grid gives 0.01352 and 0.02378,
respectively. The two finest unshifted grids differ by at most 0.02626.
The finest origin change reaches 0.12319 at s = 0.2625, exceeding the declared
0.1 diagnostic. This local residual sensitivity is retained explicitly;
there is no further grid tuning. Early jagged structure persists across methods.

See `attempt_01/results/comparison.json` (authoritative), `comparison.csv`
(table), and `comparison_original_frames.png`/`.pdf` in that result directory.
The `comparison_cell_centers` figures show the separate projection diagnostic.
`completion_review.json` records final integrity, visual inspection, the
unresolved origin-tolerance check, and result-file SHA256 values.
