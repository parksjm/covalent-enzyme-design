# Overnight coordinate-driving pilot: incomplete construction

Run: `output_nyl12_pa6_dftb_pilot_20261024/coordinate_driving_forward_reverse_attempt3/`.
All 18 forward windows completed engine minimization and stripped-coordinate
export. The controller then stopped at its forward TI1 geometry gate; no
reverse window ran. No calculation from this attempt remains active.

The final achieved OG1–C distance is 2.60843 Å versus a 1.55 Å target.
The proton coordinate OG1–HG1 minus N–HG1 is −0.73086 Å versus +0.50 Å.
OG1–HG1 is 1.11059 Å, N–HG1 is 1.84145 Å, and substrate C–O/C–N
are 1.26549/1.36715 Å. These coordinates are not the requested TI1.

The final window stopped after seven minimization cycles with global RMS
0.04279 kcal/mol/Å, below the inherited 0.05 criterion, while its maximum
force remained 6.8292 kcal/mol/Å on QM atom O49. A global RMS over this
large system is insufficient to establish local reactive-region convergence.
The chosen finite-strength restraints (25 kcal/mol/Å² in Amber's kΔ²
convention) also allow target errors. The controller correctly prevented a
reverse run from this endpoint, but per-window acceptance was too permissive.
This does not establish that TI1 is unstable: TI1 was never reached.

Energy analysis must sum full printed components, including CMAP and DFTBESCF,
without double-counting EAMBER. RESTRAINT includes positional and NMR driving
terms. Explicitly separate driving bias from the energy at the sampled
structures; removing bias does not turn those structures into a minimum-energy
path, and the original MM positional restraints remain. No free energies or
validated barriers may be reported from this pilot.

Final energy relative to window 00 is +91.7884 kcal/mol including driving
bias. The explicit driving penalty is 65.88241 kcal/mol; after removing
that term the relative value is +25.90599 kcal/mol, still including the
MM positional restraints. These are diagnostic numbers, not a barrier.
An independent coordinate-based audit verified, for all 18 endpoints, that
printed RESTRAINT equals the MM positional term plus the two driving terms
within 0.00005 kcal/mol (printed precision). Thus the poor coordinate tracking
is not explained by a missing driving term in the reported energy.

Derived plots/data and runnable scripts are retained in
`output_nyl12_pa6_dftb_pilot_20261024/overnight_review_20261025/`.
Original engine logs, structures and successful window movies are preserved.
Next: tighten and verify reactive-region convergence and target tracking in a
few windows before another full path. Do not start a reverse path from the
failed forward endpoint, infer a TS from this scan, or expand the project scope.
