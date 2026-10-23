# Truncated native-state TI1 pilot — 2026-10-23

The user authorized replacing the 103-atom full-PA6-L4 pilot with a local
substrate model and increasing the CPU allocation on violet. The old worker
received a targeted SIGTERM; its launcher records exit status 1 and an end
time. This is a user-requested supersession, not a chemical acceptance failure.
Inputs, logs, checkpoint and scratch were not manually deleted.

## Model and atom mapping

The new model contains Thr NH3+, Thr-NME, the covalent tetrahedral oxyanion,
and separate explicit H2O. The substrate fragment is
`OG1-C*(O-)(CH2CH3)-NH-CH2CH3`: formula C10H25N3O4, 42 atoms including 25
explicit hydrogens, net charge 0, singlet. PA6 bonds C70-C78 and C75-C79 are
cut. Retain C54/O49/N59/C62/C70/C67/C75 and their attached original H atoms;
add one H at C70 and one at C75 along the removed bond directions at 1.09 A.
The nearest substrate methylenes, C62 and C67, are retained. No ACE is added.

All 40 retained atom coordinates match the source seed exactly. C54 remains
R and Thr CA/CB remain S/R. Initial graph, stereochemistry and short-contact
checks pass. The water comes from the earlier physical TIP3P water-only
preconditioning, including both water hydrogens. No manual water displacement
or failed DFTB endpoint is used. This local parameter model does not change
the intended full PA6 L4 substrate for subsequent enzyme calculations.

## Calculation and acceptance

Active tmux: `nyl12_truncated_ti1_16core_20261023`.
Directory: `output_nyl12_truncated_reference_20261023/truncated_pilot_01/`.
Gaussian runs a geometry optimization only: wB97XD/def2TZVP, Tight optimization
and SCF, UltraFine integration, XQC, MaxCycles=150 and MaxStep=5. It uses 16
shared-memory processors and 32 GB RAM, pinned to CPUs 0-15, which represent
16 distinct physical cores out of 24. Gaussian confirms 16 processors and
677 basis functions, versus 1,693 in the 103-atom pilot. Thread count alone
does not establish scaling or a completion time.

The same conditional pilot policy applies: 12 heavy atoms frozen, 30 atoms
mobile, plus one AF3-seed chi1 constraint. Thr OG1, PA6 C54/O49/N59, water O,
and all H atoms are mobile; reactive bonds are unrestrained. The AF3-seed
torsions remain explicitly distinct from the n12c parameter-reference values.
No extra dispersion correction or solvent is added.

This is a state-retention pilot, not an accepted parameter reference. After
a converged candidate passes the chemistry and water checks, separately
obtain its exact Hessian and gradient and assess the constrained curvature.
No Hessian job is automatically chained. The original R cap-sensitivity and
charge/parameter-transfer gates remain unresolved. Effects of the remote PA6
chain and charged termini are outside this cluster; full-system transfer
validation remains necessary. No MD or new AF3 calculation has started.

Eight focused chemistry and corruption tests passed in tmux: 97% overall
coverage, 94% for preparation and 95% for the assessor. Ruff lint and format
checks passed. The unchanged repository suite last passed 447 tests with 13
optional fixture skips. No new production dependency or module was added.
Exact inputs, source snapshots, hashes, process logs and separate scientific
acceptance records are retained. The older 103-atom record below is history.

# Native TI1 seed and frozen-scaffold pilot — 2026-10-23

The user approved the four inspected capped reference geometries and authorized
continuation. This does not reverse the failed R cap-sensitivity gate. Map the
existing enzyme-bound seed's attack stereochemistry before choosing a reference.

Build the already specified 103-atom native graph: N-terminal Thr NH3+, full
zwitterionic PA6 L4, covalent tetrahedral oxyanion, NME cap, separate H2O; net
charge 0, singlet. All 60 hydrogens are explicit. No ACE cap and no native phi.

Use attempt05's unshifted AF3-derived preparation for every retained nonwater
atom. Use its physically minimized water-only `original_order.rst7` for all
three water atoms. This earlier, authorized Amber/TIP3P stage moved only water;
it is not the rejected manual oxygen displacement. Its 200-cycle endpoint did
not meet the minimizer convergence criterion and is only a starting seed.
Never use the subsequent failed DFTB3/MM endpoint. Preserve source atom maps,
full hashes and the distinction between generated and retained cap hydrogens.

NME N/H/CH3/one methyl H come from the following beta residue N/H/CA/HA.
Replace its two cut CA–C and CA–CB bonds with methyl H along their original
directions at 1.09 A. This is declared boundary construction, not optimization.

The first native Gaussian calculation is a **frozen-scaffold state-retention
pilot**, not a charge-fit or force-field reference. Use wB97XD/def2TZVP,
UltraFine, Tight SCF with XQC, Tight redundant-coordinate optimization,
MaxStep=5, at most150 cycles, 8 CPUs/32GB; no additional dispersion or solvent.
Freeze every heavy atom except Thr OG1, PA6 C54/O49/N59 and water O. All H are
mobile. Retain this seed's chi1 through one explicit dihedral freeze. Psi is
fixed by the scaffold. These are explicitly the AF3 seed geometry, not a
replacement for the n12c.cif psi/chi1 reference used for the capped series.
No reactive bond or proton-transfer distance is restrained.

The purpose is to learn whether the native proton inventory and tetrahedral
center survive local QM relaxation without allowing charged PA6 ends to fold
onto the reaction center after removing the protein. Frozen coordinates do
not restore the missing protein electrostatic environment; failure may require
an explicitly defined embedding/expanded cluster. Success alone cannot license
RESP fitting, force-field transfer or protein refinement. This pilot has no
frequency/Hessian claim, and the 37-atom reference gates remain in force.

Predeclared post-run checks are normal termination, optimization convergence,
charge/multiplicity, graph/stereo, frozen Cartesian displacement <=0.0001 A,
chi1 error <=0.1 degree, no O49–H distance below1.20 A, and the existing four
reactive bond windows in docs/ti1_parameter_protocol_20261022.md. Report relay
water geometry separately. Retain failures without changing thresholds.

Every preparation/check/calculation uses scripts and tmux, immutable input/code
snapshots and hashes, explicit environment and logs, with process status
separate from scientific acceptance. Active journal:
output_agent_actions/20261021.jsonl. No AF3 fork change is part of this step.

## User proposal retained for later; no MD started

After TI1 force-field validation, explicit-solvent MD could supply water-bridge
snapshots for QM/MM instead of requiring AF3 to order water. Evaluate both
donor/acceptor distances and directional H bonds, bridge occupancy and residence
behavior, and several representative snapshots rather than one favorable frame.
Water identity may exchange. TI1 chemistry must survive sampling; apparent bridge
occupancy under restraints is conditional and is not unbiased occupancy. Keep
TIP3P MM atom naming/LJ when selecting the water into the QM region and record
the exact trajectory frame/atom mapping. The user explicitly said not to run
this yet (2026-10-23). No MD or new AF3 prediction has been started.

## Preparation results and active calculation

The source-mapped seed is R at attacked C54; Thr CA/CB remain S/R.
The 103-atom graph contains60 explicit H and netcharge0. Existing graph,
stereo and short-contact screens pass. Water O moved1.408215A in the earlier
physical TIP3P stage; nonwater displacement is0.0A. Source seed psi115.362574
and chi1-38.788954 degrees differ from the n12c-derived reference values;
this pilot does not silently replace that parameterization reference.

Gaussian installation/freeze preflight completed: fixed O/H displacement0,
mobile H displacement0.074405A. Six synthetic corruption tests passed; the
post-run assessor has95% line coverage, while overall helper-bundle unit
coverage is32%. Other preparation scripts were exercised on the real inputs.
No new production module or dependency was added.

Active tmux: `nyl12_native_ti1_pilot_20261023`. Calculation directory:
`output_nyl12_native_reference_20261023/native_frozen_scaffold_pilot_01/`.
Read optimization.log and job/status.json; scientific_acceptance.json is
written separately after normal engine exit. Gaussian has started, not
converged. It uses65 mobile atoms and38 frozen atoms, plus one chi1 constraint.
No native ESP/RESP, new force-field terms, protein relaxation or MD yet.
