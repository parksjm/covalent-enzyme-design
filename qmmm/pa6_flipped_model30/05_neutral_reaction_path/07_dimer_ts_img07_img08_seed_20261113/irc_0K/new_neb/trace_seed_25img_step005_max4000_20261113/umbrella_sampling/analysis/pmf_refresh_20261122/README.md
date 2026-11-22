# NEB-seeded PMF refresh — 22 November 2026

All 63 windows contribute six accepted 5 ps production blocks: 30 ps/window,
1,890 ps total. Immutable source and scalar trajectory-data snapshots are
recorded in input_source_manifest.json. The original two-CV Amber biases enter
MBAR; all proton branches are retained without chemical censoring.

The q_acyl marginal uses the same histogram grid and fixed E–S reference
bin as the previous refresh. It is a density-reference difference, not an
integrated E–S basin free energy. The apparent maximum is 47.0475 kcal/mol
(previously 46.9414). Fifty within-window 5 ps block-bootstrap refits give a
pointwise interval; they do not capture all slow-state or model uncertainty.

Overlap is connected, but unbiased-target ESS is only 7.29 and block drift
remains in windows 028, 030 and 032. Sampling convergence and a validated
activation barrier have not been established. Proton ownership counts are raw
biased-frame diagnostics, not equilibrium populations.

The independently verified compact_export is published through the stable
analysis/pmf_latest pointer and included by the no-argument laptop pull script.
All prior simulation data and analysis exports are retained.
