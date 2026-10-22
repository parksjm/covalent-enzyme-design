# Scripted TI1 relay-water repair — 2026-10-21

## Chemistry and scope

The substrate names are PA6 L4 (`pa6_l4`) and MAXA|MAXA (`maxa_maxa`).
Historical `pa6_l2` and `poly_xama_maxa` input identifiers are provenance
aliases, not the current substrate names. L2 and MAXA remain the respective
hydrolysis products; renaming does not change their molecular graphs.

For these TI1 seeds, catalytic Thr 233 is PDB E1 / Amber residue 929. Its
N-terminal group is NH3+, OG1 has no proton and is covalently bonded to the
tetrahedral substrate carbon; the substrate oxygen is an oxyanion. The relay
must connect **N–H…Ow–H…OG1**. A pair of O…N/O distances alone is insufficient.
This protocol places neutral H2O and does not move a proton between atoms.

The original water positions fail the retained steric gate (maximum contact
shortfall 0.6866 Å PA6 / 0.6613 Å MAXA|MAXA; allowed 0.10 Å). The prior
5832-point fixed-pocket oxygen search found zero eligible sites. Excluding
only Asn 186 side-chain atoms and Thr 929 N hydrogens exposed 9 / 5 oxygen
sites. This motivates bounded local rotations, not manual coordinate editing.

## Reproducible search and acceptance

The versioned entry points are `scripts/launch_nyl12_relay_repair.py` and
`scripts/repair_nyl12_relay_water.py`. The launcher retains input copies and
source snapshots, hashes, software versions, explicit arguments/environment,
and a runnable tmux job script. Results use new attempt directories. Original
models and previous failed attempts remain intact. The live action journal
records preparation, launch, inspection, exports, and repository operations.

1. Validate the original TI1 atom inventory and connectivity. Search the same
   water-O sphere-intersection grid: distances to N and OG1 each 2.6–3.4 Å in
   0.1 Å increments, with 72 azimuths. Keep every site passing the optimistic
   fixed-environment screen, not merely the best-ranked site.
2. Rotate Asn 186 chi 1/chi 2 and Thr NH3 hydrogens in explicit 15° increments.
   These are delta rotations from the input structure. Attached hydrogens move
   with the Asn side chain; all backbone, substrate and Thr heavy atoms stay
   fixed. Require Asn ND2–oxyanion O 2.4–3.5 Å.
3. Check moved atoms against the full environment with only topology 1–2 and
   1–3 exclusions. Geometric contact limits (Å): H–H 1.40, H–C 2.10, H–N/O
   1.55 for polar H or 2.00 for nonpolar H, H–S 2.30, C–C 2.80, C–N/O 2.75,
   C–S 3.00, N/O–N/O 2.50, N/O–S 2.95, S–S 3.20. Maximum allowed shortfall
   is 0.10 Å. These are screening criteria, not a force field or energy model.
4. Require N–H…Ow angle ≥140° and H…Ow 1.4–2.5 Å. Construct rigid H2O
   (O–H 0.9572 Å, H–O–H 104.52°), aim H1 at OG1, and sample the H2 roll
   every 10°. Do not substitute an oxyanion-directed water for the OG1 bridge.
   Require all water atoms to pass the original full-environment steric gate.
5. Rank passing candidates geometrically by small Asn/heavy and NH3/H
   displacements, with a modest preference for 2.85 Å bridge distances.
   Retain sampling settings, stage failures, original and selected geometry,
   and the exact selected rotations. A finite-grid failure is not proof that
   no valid geometry exists in continuous coordinate space.
6. Independently verify fixed coordinates and non-water bond-length invariance,
   write an Amber restart, and rerun the unchanged TI1 preflight on that file.
   Export a PDB only after acceptance. Record d 1, d 2, d 3, a_BD and
   herschlag_phi before and after. Process exit status and scientific acceptance
   are separate records; completion can legitimately yield no accepted seed.

Accepted outputs, if found, are **geometric starting seeds**, not optimized
TI1 minima, TS structures, parameter validation or evidence of catalytic
activity. Fixed pre-existing distortions are not repaired by this local search.
No failed candidate is silently promoted. A remaining backbone obstruction
requires a separately documented bounded-repair protocol before refinement.

## Inspection export

`scripts/prepare_nyl12_inspection_download.py` generates canonical model names
and copies the versioned `scripts/nyl12_inspection_view.pml` into the bundle as
`protein/view_models.pml`. The viewer includes `set valence, 0`,
`show sticks, resi 113+186 and chain A`, and `util.cnc`. All exported files have
source paths and hashes. No coordinates are adjusted interactively in PyMOL.

## AF3 alternative: executed feature preflight, prediction pending

The user proposed a minimal alpha–beta/covalent-substrate/water test. Scripts
`launch_nyl12_af 3_water_test.py` and `check_af 3_water_features.py` launched the
PA6 L4 test in tmux with source/input snapshots. Explicit HOH is retained as a
non-polymer ligand despite the pipeline water-chain filter. Actual featurization
retained both N/OG1-to-water restraints. Inspect the attempt job status/log and
prediction outputs before claiming inference success. AF3 lacks explicit water
H orientations and does not validate terminal protonation.

## Authorized restrained local minimization — 2026-10-21, evening

The user approved using a roughly placed AF3 water oxygen as a starting point,
followed by explicit hydrogen construction and restrained local minimization.
Perfect AF3 water placement is not a prerequisite. Two water–Thr 233 N/OG1
distance restraints define the intended bridge; the 28 extra exclusion
restraints tried in AF3 attempt_02 are not adopted as the refinement protocol.
Attempt_01 meets its 2.6–3.2 A bridge bounds but has a 1.673 A water–C62 clash.
Attempt_02 loses the N bridge (4.847 A) and retains clashes. Both failed
predictions and independent audits are retained.

The first repair experiment uses the same alpha-A/beta-E PA6 L4 test system,
not a new full-octamer prediction. Extract the corresponding covalent topology
from the original product-gated octamer without changing its charges, atom
types or covalent connectivity. Map the AF3 heavy atoms explicitly and rebuild
attached hydrogen coordinates from local bonded frames; retain Thr 233 NH3+,
unprotonated covalently bonded OG1, intact substrate C–N, and neutral H2O.
Resolve any missing terminal heavy atoms explicitly and document their source.
A small scripted oxygen displacement may remove catastrophic overlap before
energy minimization; this geometric preparation is not an energy minimum.

Water remains the topology's **TIP3P WAT residue, atoms O/H1/H2**, preserving
its atom types and Lennard–Jones parameters for QM/MM nonbonded interactions.
Record and verify the actual source-topology types and parameters; naming alone
does not assign them. Water is included in the DFTB3 QM region, so its internal
geometry and electronic interactions are not governed by a classical TIP3P
water potential. Do not impose rigid-water SHAKE constraints on this QM water.
OBC denotes an implicit-solvent GB model, not an explicit-water atom template.
The existing protocol specifies GB8; do not silently relabel it as OBC.

Use a separate repair-stage input and acceptance record. Permit water and the
nearby pocket to relax, restrain the surrounding protein, and initially apply
two broad water-O to Thr N/OG1 distance bounds. Retain source/reference
coordinates and document the flexible atom selection, force constants, QM
boundary, convergence and final geometry. Reduce/release bridge restraints in
a separately recorded stage before treating persistence as supporting evidence.
Existing strict TI1 acceptance checks remain intact; allowing a crowded seed
into this explicitly authorized repair is not accepting it as a valid model.

No fitted force-field parameter set, accepted TI1 minimum, or octamer-compatible
repair is claimed by this decision. Launches and outcomes are recorded in the
active append-only action journal and the dated job directories.

### Superseding user decision: original oxygen and TIP3P preconditioning

The user rejected the geometric 1.12 A oxygen displacement. The active protocol
retains the exact AF3 oxygen when adding hydrogens. Shifted attempt 01 was stopped
and retained. Direct unshifted DFTB3 attempt 02 produced nonconverged SCC forces
and was stopped. The user approved a classical pre-minimization instead.

Use existing Amber/TIP3P parameters with only WAT O/H1/H2 mobile and every enzyme
and substrate atom fixed. Settings:GB8,200 cycles,drms 0.001, two N/OG1–Ow distance
bounds 2.6–3.2 A (rk 2=rk 3=10). This prepares coordinates; it does not validate
classical TI1 energies. Sander's GB belly implementation requires the moving
atoms first: retain a temporary water-first topology, per-atom parameter check,
explicit atom permutation and inverse coordinate mapping. Nonwater positions
must remain within 1e-6 A. Attempt 03 retains the GB atom-order rejection; attempt 04
retains a ParmEd AmberParm copy/LJ-table error before engine launch. Converting
to a parameterized generic Structure before reordering corrected the copy issue;
the mapping/fixed-coordinate regression passed with 94% pre-minimizer coverage.

DFTB3/MM follows in three stages:100 cycles with nonwater heavy atoms restrained
at 5 kcal/mol/A^2;300 cycles allowing nearby residues to move with outer heavy
atoms restrained at 5;300 cycles with bridge restraints released. Bridge force
constants 10,5,0; gradient targets 0.5,0.05,0.05. Nearby means whole residues with
a heavy atom within 5 A of Thr 233/ligand/water, plus beta residue 2 (29 residues
here). This allows adjustment, not deliberate pocket widening. All positional
restraints use the original prepared reference. QM mask
`:233,329,330|(:234&!@C,O)`, charge 0,109 atoms,364 electrons including link H;
sole boundaryCA234–C234. Numerical SCC controls:dftb_maxiter 200,itrmax 2000,
scfconv 1e-8; electronic temperature unchanged. Explicit inaccurate-force SCC
warnings stop the stage automatically.

Water parameters remain WAT/O,H1,H2,OW/HW,q=-0.834/+0.417/+0.417,
oxygen Rmin/2=1.7683 A,epsilon 0.152 kcal/mol. Attempt 05 is the corrected physical
run under `output_nyl12_local_repair_20261021/pa 6_l 4_alpha_beta/`.
