# TI1 parameter generation protocol — 2026-10-22

## User-selected reference and scope

The user authorized generation and validation now, then specified a truncated
ACE/Thr/NME anion with a retained tetrahedral modification and constrained
backbone/side-chain torsions. This is the first **bonded-parameter surrogate**.
It is distinct from the real N-terminal Thr NH3+/oxyanion TI1: do not transfer
its charges to that native state or declare the native force field ready.

Use the unmodified user-selected `cif/n12c.cif`, SHA-256
`e318fddc7142be1efabd6ebbd80ad6d45bf98f6d83f1ea7c1e60add9646d7b8c`.
The explicit selection is author chain A, Thr233. The mmCIF has one model and
multiple chains; other chains include sequence variants and are not pooled.
Measured psi = 134.110276428 degrees, chi1 = -58.964417878 degrees.
No preceding chain-A carbonyl exists, so native phi is undefined. ACE-cap phi
values -120 and -60 degrees are declared model choices for cap sensitivity.
Both R and S attacked-carbon configurations are explicit comparison states;
neither is silently declared the biologically preferred face.

## Local graph

Use ACE-THR(-O-C*(O−)(CH3)(NHCH3))-NME: Ac-Thr-NME plus an attacked
N-methylacetamide surrogate, with the Thr hydroxyl proton omitted to an
external reservoir. Formula C10H20N3O4−, 37 atoms, charge -1, 134 electrons,
singlet; Thr CA S and CB R. The C*–N bond remains intact. Methyl termini
represent the local PA6 amide environment; they are not a full L4 substrate.
Relay water is absent from this local bonded reference; the actual enzyme
QM subsystem will still include explicit relay water with TIP3P MM naming/LJ.

The native full-L4 graph is separately scripted: 103 atoms including explicit
water, NH3+ reacting Thr, substrate oxyanion, charge 0. It is reserved for
native-state charge/transfer validation. Real PA6 termini remain NH3+/COO−.
The current provisional protein topology remains barred from refinement.

## Calculation stages

1. Deterministic ETKDG seed construction, then set declared phi/psi/chi1.
   This constructs starting coordinates; it is not an energy optimization.
2. Gas-phase wB97XD/def2TZVP constrained optimization in Gaussian 16, with
   UltraFine integration, Tight optimization/SCF, CalcFC and XQC. No extra
   dispersion correction. Freeze only the three declared dihedrals; the
   covalent TI1 reaction-center bonds and all hydrogens remain free.
3. Check normal completion, ordered atoms, L-Thr and C* stereochemistry,
   covalent connectivity/proton inventory and frozen-torsion agreement before
   starting the fresh DFT frequency/Hessian calculation.
4. Obtain the exact Cartesian DFT Hessian in a separate checkpoint, followed
   by a separate HF/6-31G(d) MK ESP calculation/checkpoint. Surrogate RESP,
   when fitted, is for testing this surrogate; native charges require the
   actual NH3+/oxyanion coupled state and explicit cap-transfer constraints.

Constrained references are not free minima. Raw full-molecule frequencies
include constrained directions; their signs do not alone establish a minimum
on the mobile-coordinate manifold. Preserve gradients, constraints and raw
Hessians for conditional-curvature analysis before fitting/validation.

## Screens declared before fitting

Optimization screening requires the intended charge/multiplicity, all input
stereocenters/connectivity, no detected proton transfer and frozen torsions
within 0.10 degree. Reactive bond windows (angstroms) are conservative chemical
screens, not fitted equilibrium targets: OG1–C* 1.30–1.75, C*–O− 1.25–1.65,
C*–N 1.25–1.70, C*–C 1.35–1.75. Retain every failure; do not relax these gates
after observing results. Record O− contacts to cap/backbone H explicitly;
fixed dihedrals do not guarantee absence of artificial gas-phase H bonding.
The scripted checker explicitly rejects O−–H distances below 1.20 angstrom
as a conservative possible-proton-transfer screen, in addition to checking
all original bonds and newly short nonbonded contacts. This is a screening
criterion, not a fitted interaction distance.

Candidate bonded terms will be derived from QM data with documented Amber
harmonic conventions. Existing AmberTools modified-Seminario routines are
available; source/units must be checked before using them. No new dependency
is needed or authorized. Independently compare candidate MM structures and
conditional curvature/scans with held-out QM data, including the new ether
junction, oxyanion and intact amide. Prospective local geometry tolerances:
bond errors <=0.03 angstrom; angle RMS <=3 degrees and maximum <=5 degrees.
Compare cap orientations: local bond targets within 0.03 angstrom and force
constants within 20%; larger sensitivity requires a larger reference model.
Held-out scan energies should have RMSE <=1 kcal/mol within the tested local
basin. These are project acceptance criteria, not literature-certified bounds.

Native-state RESP must separately enforce exact total charge (residual <1e-6),
document standard peptide/cap constraints, use shared multi-conformer charges,
report fit/held-out ESP errors and retain the true reactive proton inventory.
Do not force arbitrary integer partial charges on the two covalently joined
residues. Complete transfer/geometry/torsion/curvature checks and exact-topology
digests before declaring parameters eligible for AF3-model MM relaxation.

All preparation, engine calculations and real-input checks use scripted tmux
jobs, source snapshots/hashes, explicit resources/environments, process logs
and separate scientific acceptance records. Current action journal remains
`output_agent_actions/20261021.jsonl`. AF3 fork changes remain the user's work.

## Conditional-curvature implementation and prospective thresholds

Before any completed TI1 Hessian was available, the following assessment was
implemented and checked against analytic constrained problems. With frozen
angular constraints g, use L = E - lambda.g, obtain lambda by solving
J.T lambda = grad(E), and project H(E) - sum(lambda_i H(g_i)) into the
constraint tangent after removing six independent rigid motions. The 37-atom
nonlinear reference therefore has 102 tested internal directions. A projection
of the raw energy Hessian alone omits the nonlinear constraint-force correction.
Eigenvalues use an orthonormal, unweighted Cartesian basis (Hartree/bohr^2);
they are not mass-weighted frequencies. No SciPy dependency was added.

Prospective numerical settings: central finite differences at 0.0002 and
0.0001 bohr; relative Frobenius changes in J and H(g) <=0.001; maximum
Lagrangian-gradient residual <=0.00003 Hartree/bohr; every reduced eigenvalue
>0.000001 Hartree/bohr^2 at both steps. Values within the sign tolerance are
unresolved, not accepted minima. Require the actual Gaussian Cartesian gradient;
never substitute zero. Require accepted optimization/frequency identity,
connectivity and torsions before this analysis. These are project thresholds,
not literature-certified error estimates. Passing establishes local stability
only under the declared frozen torsions, not a free minimum or validated FF.

The four interim geometries passed an independent derivative-only implementation
check: relative step changes below 8.7e-8 and 102 remaining directions. That
check uses no energy Hessian and cannot assess their energetic stability.
On 2026-10-22 the user visually inspected the interim optimization trajectories
and judged them good. This does not waive convergence, cap-artifact or parameter
transfer checks, including the recorded S-reference cap N-H/oxyanion contacts.

## Native-state atom mapping and proposed RESP constraints

Scripted R/S graph crosswalks identify 14 corresponding heavy atoms between
the capped bonded surrogate and the 103-atom native Thr-NH3+/PA6-L4/oxyanion/
NME/water graph. Surrogate methyls become native methylenes at two boundaries;
no hydrogen correspondence or capped-anion charge transfer is implied.
The native graphs have no accepted 3D geometry or ESP data yet.

A source-backed proposal protects the installed NTHR peptide C/O charges
(+0.6163/-0.5722), the full NME cap, and TIP3P water (-0.834/+0.417/+0.417).
These 11 fixed atoms sum to +0.0441 e; the remaining 92 atoms before equivalence
must sum to -0.0441 e for total charge zero. Fit the reactive ammonium and
substrate jointly without arbitrary integer residue charges. Proposed
ammonium/methyl/carboxylate equivalences are explicit in the retained proposal;
methylene hydrogen equivalence is not assumed. The retained Thr proton origin
HG1 is bookkeeping provenance for the native ammonium H3, not an O-H bond.

This is a proposed charge-fitting boundary, not an executed or accepted fit.
Fixed TIP3P charges while water participates in the QM cluster require fit and
held-out ESP-error validation. Water keeps TIP3P names/types and LJ parameters
in the MM topology; its electrostatics inside a QM/MM calculation follow the
QM/MM embedding policy rather than an extra duplicate MM charge contribution.

Actual installed ff19SB/TIP3P leaprc files were sourced with tleap, and NTHR,
NME and TP3 were independently exported and parsed with warnings treated as
errors. All 11 proposed constraints match those clean source units for both
native graphs. The first whole-library ParmEd read emitted warnings for an
unrelated CHCL3BOX legacy unit; preserve that log, but base source acceptance
on the subsequent clean target-unit extraction/verification. No library was
modified. Reports and exact helpers are under
output_nyl12_parameter_validation_20261022/{parallel_prerequisites_01,
native_library_extract_01,native_library_verification_01}/.

## Amber tools and force-field baseline — explicit clarification

Source audit on 2026-10-22 confirms that the legacy topology helper
`parameterize_poly_xama` uses original GAFF (`antechamber -at gaff`,
`parmchk2 -s gaff`, `leaprc.gaff`) with provisional AM1-BCC or Gasteiger
charges on the intact substrate. It hard-codes a neutral intact substrate.
It has not parameterized the new coupled TI1 state with Gaussian RESP and
must not be reused unchanged for the capped anion or native adduct.

The intended new candidate baseline is ff19SB for unmodified protein and
GAFF2 for ordinary substrate chemistry, with TIP3P water. Retain standard
ff19SB protein atom types/terms; do not antechamber-type the entire protein.
Use explicit GAFF2 selection consistently in antechamber, parmchk2 and LEaP;
record the exact installed parameter-file hashes. This is a candidate method
selection, not an executed migration or an accepted new force field.

Use Gaussian wB97XD/def2TZVP for reference geometry/Hessian and separate
HF/6-31G(d) MK ESP data for Amber RESP charges. Antechamber supports Gaussian
output and charge/type preparation; espgen/respgen/resp are the standard
AmberTools building blocks. The native multi-conformer, fixed-cap/backbone/
water constraints require retained explicit RESP control files and charge
mapping. A default one-command ligand RESP fit does not by itself implement
those constraints. Surrogate and native charge fits remain separate.

First assign and test the standard ff19SB/GAFF2 candidate, including review of
all TI1 atom types and parmchk2 estimates. Refit missing or demonstrably
inadequate reactive-region bond/angle/torsion terms against QM data, without
replacing unrelated validated protein terms. The optional modified-Seminario
candidate extraction and its manufactured convention tests are additional
checks, not a claim that antechamber normally derives all bonded parameters
from a Gaussian Hessian. Do not emit bespoke parameters merely because a
Hessian block is available. Acceptance still requires cap sensitivity,
independent local energy/geometry checks and native-state transfer validation.

Primary sources inspected: installed AmberTools26 source/CLI and
https://ambermd.org/antechamber/ac.html plus
https://ambermd.org/tutorials/basic/tutorial4b/index.php.
