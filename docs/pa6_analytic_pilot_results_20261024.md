# PA6 analytic-gradient ES pilot — 2026-10-24

The second bounded Sander DFTB3/MM minimization completed 25 cycles with process
exit 0. It **did not converge**: final printed RMS was 2.4646 kcal/mol/A against
the declared 0.05 target; GMAX was 38.646 kcal/mol/A. Geometry acceptance remains
false. Finite coordinates and the fixed periodic box were preserved. These are
diagnostic ES coordinates, not an accepted minimum, TI/TS or mechanism result.

Retained attempt:
`output_nyl12_pa6_dftb_pilot_20261024/analytic_relaxation_attempt2/`.
The original inputs, executed sources, analytic-gradient outputs, 25 saved
coordinate states and endpoint review remain there. This pilot used analytic
Sander gradients; the user stopped further finite-difference calculations.
A separately authorized continuation with a 500-cycle budget completed in
tmux `nyl12_es_continuation_20261024`, retaining the original positional
reference. It exited zero after 107 cycles with RMS gradient 0.042903,
below the 0.05 target, and saved 21 optimizer snapshots. The periodic box
was unchanged and coordinates remained finite. Its artifacts are in
`analytic_relaxation_continuation_attempt1/`; the original attempt2 exports
remain the earlier unconverged pilot. Final endpoint review retains Thr OH/NH2:
OG1–HG1 = 0.98493 Å, N–H = 1.02118/1.01706 Å. Substrate C–O/C–N
distances are 1.25974/1.35877 Å. Maximum QM/MM heavy-atom shifts from the
original reference are 0.28286/0.32613 Å. This is a usable restrained ES
reference; no TI bond formation or mechanism validation is established.

| Measurement | Initial | Endpoint |
|---|---:|---:|
| Thr929 OG1–C54, A | 3.1887 | 3.1180 |
| Phe113 N–O49, A | 3.0149 | 2.9607 |
| Asn186 ND2–O49, A | 2.6560 | 2.8027 |
| C54–O49, A | 1.2243 | 1.2540 |
| C54–N59, A | 1.3432 | 1.3613 |
| OG1–HG1, A | 0.9600 | 0.9960 |

The two N-terminal N–H distances remain approximately 1.015 A. Maximum QM
heavy-atom displacement was 0.1932 A. These measurements describe the supplied
ES model; they do not establish a catalytic proton-transfer mechanism or
protonation populations.

Printed energy components, including CMAP, DFTBESCF and restraints, sum to
-332062.9866 initially and -379914.9706 kcal/mol at the endpoint. Their difference
(-47851.9840 kcal/mol) is whole-system minimization accounting, **not a reaction
energy or barrier**. The initial sum agrees with the API energy within the
printed precision. Use the component sums rather than the coarse NSTEP header
for this comparison. The unconverged endpoint does not justify energetic ranking.

The laptop export includes the local pocket movie, complete protein/substrate
endpoint viewer, QM subsystem XYZ files with atom mapping and endpoint review.
All original solvated coordinates remain on the workstation. Viewer bonds and
ES labels do not certify convergence or chemical-state validity. Further
scientific work requires review of the endpoint, restraint/embedding treatment
and convergence before extending the calculation.

Evidence: `process.json`, `geometry_acceptance.json` and
`inspection/endpoint_review.json` in the retained attempt. No accepted stationary
point or full-cycle chemical mechanism is claimed.
