### TS diagnosis: recovered SCC trials were misclassified

Installed Amber source distinguishes `dftb_maxiter=70` inner/Broyden trials
from the overall `itrmax=1000` budget. The minus/plus descent logs contain
5/9 failed inner trials, respectively, followed by successful convergence at
0 K. There are no unrecovered/fatal SCC events in those logs. The old generic
regex treated those recoveries as terminal failures. A shared parser now
accepts only the exact recoverable message with subsequent same-step native
convergence, retains the retry history, and rejects unrecovered/fatal errors,
NaN, inaccurate-force warnings and electronic-temperature changes. Independent
process, finite-energy, native-termination and force checks remain required.

A bounded five-evaluation diagnostic completed with all native gates passing:
three fresh TS evaluations have identical gradient hashes and energies; fresh
minus/plus endpoint maximum forces are 9.38e-6/9.51e-6 Eh/bohr. They classify
as E-S and acylenzyme, respectively, and lie 46.41/43.70 kcal/mol below the
candidate. HG1 is on OG1 in E-S and on the leaving substrate nitrogen in the
acylenzyme; Thr233 N remains NH2. These are mode-directed downhill minima,
not a formal IRC or a certified TS. The high candidate-to-E-S difference
still does not establish a biochemically viable route.

The candidate maximum force is 9.24e-5 Eh/bohr and meets its **original**
1e-4 criterion. The diagnostic's uniform 1e-5 flag is a tighter reference,
not evidence that the original TS force-convergence test failed. There is
no current reason to change TS optimizer or repeat a long optimization first.

The remaining blocker is Hessian precision: one strong negative mode exists,
but the conservative global screen leaves 4,535 soft signs unresolved.
Most nonreciprocity is localized around Thr930 HB; identical same-point forces
establish determinism, not differentiability or symmetric-Hessian accuracy.
The eigensolver residuals are already small; increasing its iteration budget
will not resolve inaccurate matrix entries.

**Next TS step (recommended, not yet launched):** at the unchanged candidate,
recompute the three dominant Thr930 HB Cartesian columns at
h=6.25e-5 and 3.125e-5 Angstrom: 12 native analytic-gradient evaluations.
Require a derivative plateau and improved mixed-derivative reciprocity before
using any replacement columns in a separate derived Hessian. Then reassess
low-mode signs and measured uncertainty. Do not lower the global threshold
just to obtain acceptance, infer an error bound from skew alone, shrink the
12 Angstrom movable region, change protonation, introduce smearing, or repeat
the entire Hessian without this bounded gate. Any later Hessian at a newly
optimized geometry must be recomputed for that geometry.

Evidence: `.../validation_loop_20261121/diagnostic_repeat_20261123/`
contains `results_with_amberhome/diagnosis.json`, `interpretation.json`,
`job_with_amberhome/`, frozen sources and all native force evidence. Original
Hessian/validation pointers and the old false-negative receipts remain intact.

# Nyl12 stationary-point validation: native force precision repair

The neutral Thr233-N / ASP272 / ASP274 direct PA6 acylation candidate remains
**unverified as a first-order saddle**. This work tests that existing candidate;
it does not introduce a new protonation state, substrate X-ray structure,
mechanism campaign, smaller movable region or assumed stable TI1.

The candidate's native zero-smearing DFTB3/MM dimer optimization converged at
cycle99: energy −718.134034102393 Eh, maximum midpoint gradient
9.2423e−5 Eh/bohr versus1e−4, RMS4.9186e−6. The QM charge is−2. All4,165
movable atoms (12,495 Cartesian coordinates) and the frozen exterior remain.
Native job completion, force convergence, electronic convergence, Hessian
inertia and chemical endpoint connectivity are recorded separately.

## Measured precision repair

DL-FIND's `write_xyz` writes coordinates **and force vectors** with seven decimal
places (`3f15.7`). An isolated library rebuild reuses all original objects and
recompiles only this writer with `3es25.16`. Installed Sander/library files and
running MD/US calculations are untouched. The small
[`output precision patch`](../patches/dlfind_xyz_output_precision.patch) and
[`source/hash identity`](../provenance/nyl12_stationary_validation_20261121.json)
are reviewable; full engine sources/objects/binaries remain workstation artifacts.
The retained source identifies DL-FIND as LGPLv3-or-later.

Fourteen actual native gradients established exact printed energy parity and
maximum gradient differences1.78551e−7 Eh/Å, within the old-format rounding
bound1.78557e−7. This repairs an identified serialization limit, but the original
electronic accuracy still failed derivative stability: chemical-coordinate
h=.002 versus.001 Å columns differed by.00149785 Eh/(Å² amu), while two remote
MM columns differed by6.93e−7 and2.05e−6. The failed acceptance is preserved.
It does not prove SCC error: central-difference truncation, cutoff behavior and
electronic-state sensitivity remain possible explanations.

Native SCC source shows original `scfconv=1e-10`, `tight_p_conv=0` implies
`density_conv=5e-7`. The bounded continuation uses `scfconv=1e-12`,
`tight_p_conv=1` (density tolerance1e-12), verbose native SCC rows, and chemical
column steps.008,.004,.002,.001 Å on the same Hamiltonian/geometry. It reports
absolute/relative differences and half-step scaling before any full matrix.
Selected-column stability is an engineering screen, not a global Hessian
error bound. The full Hessian's sign margins and symmetry diagnostics remain
separate scientific gates; no universal1e−4 curvature tolerance is asserted.

The first tighter-density attempt exited before any usable gradient. Its
failure log/status and immutable sources remain. Attempt02 uses short native
topology filenames, preserves process-failure logs before removing temporary
scratch, and retains the same scientific inputs. Exact SCC warnings
(`Convergence could not be achieved`) fail even when native Sander continues.

Attempt02 completed all eight native gradients:26 SCC iterations each, final
energy changes1.3–3.0e−13 and charge changes3.7–4.9e−13. Its chemical-column
error.00149784919005 is unchanged from the original.00149784918991, excluding
loose SCC as the explanation for that discrepancy. The initial finest difference ratio
was.23876, but finer measurements establish a plateau rather than continued
quadratic scaling. Errors h=.0005→.00025 and.00025→.000125 Å are2.36e−7 and
3.78e−7 Eh/(Å² amu), relative1.04e−6 and1.67e−6. The smallest step increases
noise. Select the intermediate.00025 Å; larger-step behavior was nonmonotonic.

An actual two-gradient native frame audit passed: both ±.00025 Å calculations
have exactly one coordinate and one force frame. After the exact source-defined
coordinate scale.999999999068864, maximum active-coordinate error is1.42e−14 Å
and signed perturbation-pair error9.49e−15 Å. Source ordering confirms the same
gradient call precedes both coordinate/force outputs, before any optimizer step.
Attempt05 enforces this coordinate/frame-count gate on every new gradient and
archives native coordinate frames. A wrong-coordinate regression is tested.

The intermediate attempt03 stopped on a validator false positive: an unbounded
`NaN` regex matched a random scratch filename containing `nan`. Word-boundary
matching and an explicit filename regression fixed it. Attempt04 completed the
fine diagnostic but its asymptotic-ratio rule held the matrix at the numerical
plateau. Attempt05 uses the measured stable three-step plateau and actual frame
identity; all earlier evidence remains. No failed attempt was overwritten.

## Persistent conditional controller

Canonical workstation tree:
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/ts_image08_attempt1/validation_loop_20261121/`.
Use its `current.json`; active immutable attempt is
`verified_frame_validation_attempt05/`, tmux
`nyl12_stationary_verified_frame_attempt05_20261121`.

1. Attempt05 reuses the complete immutable derivative diagnostic. Three fine
   steps establish a stable plateau; its intermediate h=.00025 Å is selected.
   This selected-column screen does not certify every Hessian derivative.
2. Only a passing screen runs native parity and three half-step column checks.
3. Only passing precision starts all12,495 Hessian columns, up to eight nice,
   single-thread native workers; completed columns have restartable hashes.
4. The complete mass-weighted Hessian is symmetrized with antisymmetric norm,
   all eigenvalues and lowest-mode residuals retained. Unresolved signs stop
   TS acceptance; a small QM-only Hessian cannot certify this full region.
5. Only numerical index one permits both signs of its actual negative mode
   followed by native endpoint minimization. Both starts must lower the energy;
   endpoints must converge forces and recover E–S/acylenzyme bond/proton geometry.
6. Complete endpoint Hessians must establish numerical minimum inertia before
   the overall scientific acceptance can pass. This is local downhill
   connectivity, not a formal IRC, PMF or experimental mechanism assignment.

Entry points are `scripts/run_stationary_accuracy_diagnostic.py` and
`scripts/run_stationary_validation.py`. Launches execute immutable snapshots,
not subsequently edited checkout files. Every accepted gradient retains exact
native force/mdout/mdin/log gzip, a lossless gradient array, and input/coordinate
hashes. Full-system derivative restart copies are explicitly temporary scratch;
coordinates reconstruct from immutable reference plus signed column displacement.
Original source inputs and historical attempts remain intact.

Before launch:48CPUs, load≈25,1.24TB free. Historical native gradients measured
4.4–7s; one full Hessian requires24,990 gradients (~31–49 CPUh native time),
plus setup/compression/eigendecomposition. Eight-worker wall time is several
hours if throughput holds; endpoint Hessians add two analogous stages.
Use actual `ts_hessian/progress.json` for measured ETA after launch. No TS,
minimum-inertia proof or full-loop completion is claimed in this launch note.

Focused analytic/negative tests cover mass scaling, whole-matrix inertia,
unresolved signs, exact SCC failure text, native failure archival, chemical
endpoint identity and conditional stage execution. Ruff and pytest results are
recorded in the append-only action journal; no production dependencies added.


## Completed Hessian outcome — November 21, 2026

All 12,495 columns completed in approximately 6.4 hours. The stored symmetric
matrix has one negative eigenvalue (-0.230729 in the declared mass-weighted
units) and 12,494 positive eigenvalues. The negative mode is approximately
1307i cm^-1, dominated by catalytic HG1 transfer and coupled to OG1 attack and
C-N elongation. This is a promising chemically relevant TS candidate.

Final verification remains incomplete. The global consistency gate used the
raw antisymmetry infinity norm as a heuristic threshold; 4,535 eigenvalues
fell within that threshold, not below zero. This is not proof of 4,535 unstable
or inaccurate modes, and antisymmetry alone is not a rigorous symmetric-error
bound. Most skew is concentrated in a few coordinates (95% of squared skew
norm in the top 20 rows), especially adjacent QM Thr930 HB. Tiny eigenpair
residuals validate matrix diagonalization, not the estimated physical Hessian.
The original failure record remains unchanged. The controller stopped before
downhill endpoint optimizations; no verified reaction connectivity is claimed.

The retained completion audit is `completed_hessian_audit_20261121.{json,md}`
in the validation campaign. Next use bounded directional/step checks and
candidate-mode diagnostic descents, preserving the Hamiltonian and 12-Angstrom
movable region. No complete Hessian rerun is justified by this audit alone.
