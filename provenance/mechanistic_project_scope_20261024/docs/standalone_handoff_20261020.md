## Proposed separate WT PA6 mechanistic project — 2026-10-24

User supports considering a full chemical-cycle DFTB QM/MM free-energy study
as a standalone paper. Focus PA6, optionally poly(MAXA), without expanding PA66.
Use hybrid DFT only for selected PES optimizations/checks, not free-energy
sampling. QM/MM validation and convergence are substantive research tasks.
Prior antechamber/AM1-BCC WT noncovalent PA6 L4 MD (3 × 100 ns) already exists;
corrected records show PA66 mainly loses local geometry/register, without
necessarily dissociating. Pin valid topologies/trajectories before reuse;
do not repeat simulations or change charges without an explicit rationale.
NylB is a comparator, not an authorized pivot; its published 2014 reactive
QM/MM study is DOI 10.1021/jz500323y. Scope and validation gates, including the
distinction between a reactive QM core and state-specific MM parameters, are
in docs/nyl12_mechanistic_project_scope_20261024.md. No new simulation started.
A full-cycle mechanism would support a paper; it is not required for every
exploratory sequence. A PMF maximum is not automatically a PES saddle point.
Updated scope and references are included in the existing laptop pull package.

## Literature/manuscript checkpoint — 2026-10-24

User requested preserving relayrefs and recent2025/2026 Baker/othergroup
lessons for Nyl12 redesign. R23--R33 and new BibTeX records now retained in
manuscript/. Introduction/Methods/sourceledger updated; proposals remain
labelled. New docs/recent_enzyme_design_lessons_20261024.md compares final
2025serine/2026metal, Oct14proteases, Oct15RFD3, nativefold/PET alternatives.
Primary-access/preprint limits explicit. Pseudo-ring excluded only in tested
fixedmodel, not a universal/transientrelay. No new QM/design/MD campaign.
Recommended boundedchemistry models and multistatescreening with fulloligomer
conformers/registers/termini; trueTS1 validation remains required for saddle/
barrier claims, not a prerequisite to every exploratory structural screen.
No new dependencies/software adopted. Laptop export includes updated drafts,
references, machine-readable BibTeX and both research documents.

## Mechanistic interpretation extended — 2026-10-24

User reconsiders frontwater and requests literature beyondNylC. Found positive
primaryproteasome QM/MM example Wei2013 DOI10.1021/jp405337v and competing
water-free protontransfer model in penicillinacylase Chilov2007. NylC2023
FEBS DOI10.1111/febs.16755 proposes Asp-Asp-Thr network for maturehydrolysis;
its precursorwater/autocleavage figures must not be conflated with catalysis.
No definitiveNyl12 mechanism assigned. See newmechanism evidence doc.
Oppositehemisphere diagnostic also0/225000 full andheavy-only. No newseed
orGaussian. Assess currentAF3-derivedconformation andmapliteraturecatalytic
residues before forcingwater orfittingassumedNH3+/oxyanion state.

## Resumed feasibility search completed — 2026-10-24

User resumed frontwater placement and parameterization/testing if feasible.
No acceptable front oxygen in fixed complete protein/substrate pocket:
0/225000 samples at2.6--3.2A N/OG1 contacts; heavy-only0/225000 as well;
expanded3.3A diagnostic0/302760. Best fullmargin−0.684A (LIGH23), heavy-only
−0.217A (substrateC70 plusAsp272 contacts). No new Gaussian or parameterfit.
RotatingNH3 or waterH cannot rescue an oxygen-heavyatom clash. This is a
finite fixed-conformer search, not an impossibility proof for the enzyme.
See docs/front_water_feasibility_20261024.md and
output_nyl12_front_water_20261024/. Full-pocket viewer uses two rejected
oxygen probes, not seeds. Next resolve protein/substrate geometry or choose
another validated conformation without bypassing TI1-FF validation requirement.
Eight SMDacceptance corruption tests nowpass, including lostbridge and reverse
protontransfer, with93.94% assessor linecoverage. Source hashes, journal and
all failed attempts retained. Earlier PAUSED/ACTIVE notes below are history.

## PAUSED at user request — evening 2026-10-23

No Gaussian process or tmux calculation from this effort remains running.
Other three existing tmux sessions are unrelated and untouched. Stop work.
Latest user XYZ is preserved exactly and checked; it clears the substrate but
clashes with Asp topology272 backbone in the original full pocket. Resolve
full-pocket geometry before treating that water placement as enzyme-compatible.
No seed/reference/FF has been accepted. No automatic restart or continuation.
SMD remains a diagnostic proposal with one user-rejected-placement pilot stopped;
that interrupted run does not establish whether continuum stabilizes TI1.
See detailed notes below and retained output roots. Active journal remains
output_agent_actions/20261021.jsonl. Unrelated untracked cif/ is untouched.

## User-supplied manual water seed checked — 2026-10-23

Saved exact supplied XYZ: output_nyl12_user_water_seed_20261023/native_seed_jmp.xyz.
Script/check_job/seed_checks.json preserve measurements and input hashes.
63atom inventory, graph and S/R/R stereo pass. WaterO-N2.67154A,
waterO-OG12.90269A; waterH-OG1 angle166.53deg, best NH-water117.73deg.
Nearest substrateheavy3.629A. However FULL original AF3-derived alpha/beta
pocket has severe clashes: ASP topology272 CA1.44314A, HA0.81751A from
waterO. This is model-context evidence, not experimental steric assignment.
Do not treat the isolated-cluster placement as enzyme-compatible. Native
pocket geometry/possible rearrangement must be resolved separately. User
coordinates remain unchanged; no new Gaussian calculation was launched.

## STOP requested: SMD seed water rejected by user — 2026-10-23

User inspected seed and requires water IN FRONT of amino terminus, away
from the substrate-side placement shown. Identity-verified Gaussian worker
4055951 sent SIGTERM; interruption record in SMD pilot. No replacement run.
Check job/status.json for terminal record. Preserve all files. This rejects
the water geometry, not SMD as a method. Earlier ACTIVE note is superseded.

The deterministic search was insufficient: only water moved, so original
NH3 proton directions biased donor selection toward H3; geometry screening
used covalent-radius contacts in truncated63atom cluster, not full-pocket
van der Waals clearance. Distance/angle satisfaction alone is not acceptance.
Next place/check against FULL substrate/protein context, define exposed face
in an atom-based frame, and allow NH3 H reorientation with water placement.
Do not restart Gaussian with another unreviewed side/clearance assumption.
User's five-minute departure constraint remains relevant; no extended work
or autonomous continuation requested. SMD/RESP protocol distinction preserved.

## ACTIVE SMD/water-bridge diagnostic — 2026-10-23

User requested a decision and tmux launch before leaving in5min.
New root: output_nyl12_smd_water_bridge_20261023/.
ACTIVE tmux: nyl12_bridge63_smd_20261023; pilot svp_pilot_01/.
63atoms, neutral singlet, wB97XD/def2SVP SCRF=(SMD,Solvent=Water),
16cores32GB,20cycles maximum, same frozen21heavy plus chi1.
No reactive-bond/proton/water restraints during optimization; no Hessian/RESP.

Prior water was not a validated bridge. New deterministic search moves ONLY
water O/H/H, retains every nonwater coordinate, and screens graph/clashes.
1800 orientations evaluated; selected N-Owater=OG1-Owater=2.70A,
N-H3-Owater=151.38deg, Owater-H-OG1=180deg, OH0.96A/HOH104.5deg.
Initial63atom geometry passes. Final acceptance adds both water heavy-atom
contacts2.5--3.3A and best NH-water / waterH-OG1 angles>=140deg, alongside
existing NH3, oxyanion, graph/stereo and convergence checks. A retained
state would be a constrained aqueous-continuum cluster diagnostic, not a
native enzyme minimum or accepted FF reference. Water continuum is not a
measured enzyme dielectric. Do not autoextend or fit parameters.

Ruff passes and real seed/source geometry checks pass. New bridge-gate unit
coverage has not been measured under the user's5min deadline. Earlier12tests
apply to the preceding gas-phase pipeline, not these added water criteria.
Source snapshots, input hashes and action journal retained. Job manifest pins
launched assessor; separate assessment_source_hash_at_checkpoint.json discloses
preformat digest in earlier staging source_manifest. No method silently changed.

Literature: Meadows etal2026 JCTC22,4012--4025 DOI10.1021/acs.jctc.6c00162
explicitly uses implicit water to stabilize glycine zwitterion for GAFF/RESP
parameter variants: https://pmc.ncbi.nlm.nih.gov/articles/PMC13130867/ .
Enzyme SMD dielectric4 precedent: https://doi.org/10.3389/fchem.2018.00249
(singlepoint enzyme energies; not proof SMD Hessians transfer to Amber).
Recommendation: solvated geometry is defensible; ESP charge protocol and
bonded fitting remain separate decisions. Do not adopt fully solvated ESP
charges or raw continuum Hessian as automatic ff19SB/GAFF2 parameters.
No replacement RESP or force-field fit launched. Diffuse/larger-basis and
actual protein-environment checks remain pending if TI1 survives.

## Pilot ended: TI1 state loss; trajectory available — 2026-10-23

The63atom donor/water pilot ended17:10:33EDT at20cycle limit, exit1.
No convergence and no TI1 acceptance: OG1-C54=1.90553A, N-H3=1.40361A,
OG1-H3=1.03419A, consistent with reverse proton transfer. WaterO-N=3.30394A
and waterO-OG1=3.53239A; an acceptable pseudo-ring is not established.
No continuation is running. Earlier ACTIVE entries below are historical.

Exported21complete Input orientation tables,63atoms and37explicitH each,
with both water H. Workstation snapshot:
output_nyl12_oxyanion_cluster_20261023/trajectory_snapshot_01/.
The usual rsync pull now includes
qm/ti1_oxyanion_cluster/optimization_progress/view_optimization.pml,
optimization_trajectory.pdb and .xyz, CSV metrics and acceptance report.
Run @view_optimization.pml then mplay from that directory. This is a
snapshot, not a live feed; fixed PDB bonds do not prove connectivity survived.
Source log snapshot/hash, exporter and tmux records preserved. Ruff and
PDB frame/atom/H roundtrip checks pass. No new QM or force-field run launched.

## ACTIVE: oxyanion-donor/water cluster pilot — 2026-10-23

User approved N-methylacetamide for Phe113 backbone NH, acetamide for Asn186
side-chain amide, and relay water, then selected omegaB97X-D/def2-SVP.
New cluster: output_nyl12_oxyanion_cluster_20261023/seed_01/; 63 atoms,
C15H37N5O6, charge 0, singlet, 37 explicit H. Components: 39-atom TI1,
3-atom water, 12-atom NMA, 9-atom acetamide. Original 42 seed coordinates
unchanged; protein frame verified without alignment; R at C54 retained.
Phe H--O49=1.7133 A; Asn HD22--O49=2.1559 A. Seed geometry screen passes.

ACTIVE tmux: nyl12_oxyanion63_svp_20261023, started 16:40:57 EDT.
Job: output_nyl12_oxyanion_cluster_20261023/svp_pilot_01/job/.
Gaussian input/policy/log in parent svp_pilot_01/. Uses 16 physical cores,
32 GB, ReadAllGIC, MaxStep=2, at most 20 optimization cycles; no Hessian,
extra dispersion or solvent. Short state-retention pilot, not a force-field
reference. All 37 H and 5 reactive/water heavy atoms move. Freeze original
12 scaffold heavy atoms, 9 added donor heavy atoms at protein positions,
and source AF3 chi1. No reaction-bond/proton restraint. Donor heavy freezes
represent held protein geometry, not a freely optimized isolated cluster.
Remaining protein electrostatics are absent; enzyme TI1 is not established.

Twelve focused tests pass in scripted tmux; Ruff passes. Source line coverage:
assessor 93.44%, donor helper 93.02%, cluster builder 94.52%, pilot preparer
91.94%. New gate checks all three ammonium N--H (0.85--1.25 A) and their
OG1--H distances (>=1.20 A), in addition to graph/stereo/oxyanion checks.
Final assessment runs even if Gaussian exits unsuccessfully. A cycle-limit
exit is not convergence; retention over a short pilot is not a proven minimum.
Inspect actual job status and log before continuing. No RESP, FF acceptance,
protein refinement, MD, AF3 or redesign campaign was started. Earlier stopped
water-free39 data remain preserved and rejected as a TI1 reference.

Seed construction/preparation/test scripts, source snapshots, hashes and all
failed preparation attempts remain in the new output root. Routine fragment
scripting delegated to Sol High and reviewed. Active journal unchanged:
output_agent_actions/20261021.jsonl. Unrelated cif/ and other tmux untouched.
The laptop package includes qm/ti1_oxyanion_cluster/view_seed.pml and seed,
plus pilot/input.gjf and pilot/pilot_policy.json. Next: assess state retention,
not automatic extension, charge fitting or force-field transfer.

## STOPPED: water-free TI1 reverses to reactant-like state — 2026-10-23

User requested stop after observing reverse proton transfer. The verified
worker3978424 receivedSIGTERM; jobended16:20:27EDT/20:20:27UTC,exit1.
All files preserved. No Gaussian job from this effort is now running and no
replacement model/method was launched. Earlier active-job notes are history.

Stopped19-frame analysis confirms ThrN–H3 lengthens1.010->1.488A while
OG1–H3 shortens1.993->0.916A; OG1–C54=2.139A, C54–O49=1.220A lastframe.
This is not a TI1 fitting reference. See
output_nyl12_waterfree_reference_20261023/state_loss_review_01/report.json.

Source-protein candidate oxyanion donors: Phe113 backboneNH (H–O49=1.713A)
and Asn186 ND2-HD22 (2.156A). Proposed next model restores relaywater plus
these actual donor fragments, retains shortsubstrate, and defines protein
boundary/embedding treatment. Exact caps still need mapping; no new run.
See docs/ti1_state_loss_and_qm_options_20261023.md for the measured evidence,
proposed~63atomlocalcluster, proton-transfer screening improvement and QM-cost
options. Basis/functional pivot is not yet adopted. KeepexistingTZVPdata;
considerlower-costbasis with samefunctional and targetedlarger/diffusebasis
checks. No RESP/forcefield/MD accepted. Livejournal unchanged.

## Trajectory inspection warning — 2026-10-23

A 14-frame snapshot is exported at
`output_nyl12_waterfree_reference_20261023/trajectory_snapshot_01/`.
At its last frame OG1-C54=2.02549 A and C54-O49=1.21644 A; these exceed
our TI1 acceptance windows, and frames13/14 fail the graph geometry screen.
This is an in-progress snapshot, not a converged chemical-state conclusion.
Do not fit parameters from this candidate without resolving state retention.
Input graph connectivity in the PDB is for display only. The movie, XYZ,
per-frame metrics and plain-command PyMOL viewer are in the laptop package.
The user requested proportionate work and fewer tokens/credits for simple
tasks: use existing tools first and avoid expanding small tasks unnecessarily.

## Water-free 39-atom candidate active — 2026-10-23

The user authorized removal of water for the local parameter reference.
The former 42-atom job FAILED at 13:49 EDT: five energy/gradient points,
23 min 52 s, l103 “Error imposing constraints”, exit 1. Preserve its entire
truncated_pilot_01 directory and failure_assessment_20261023.json. There is
no accepted geometry; this is a numerical constraint failure, not a conclusion
about TI1 stability. The old 103-atom job remains stopped by user request.

ACTIVE tmux: `nyl12_waterfree39_gic_16core_20261023`.
Directory: `output_nyl12_waterfree_reference_20261023/waterfree_pilot_01/`.
39 atoms, C10H23N3O3, 23 explicit H, net charge 0, singlet. Thr NH3+, O49-
and C54 R retained. Remove exactly water O/H1/H2 from original 42-atom seed;
all remaining coordinates unchanged. Do not restart from the failed endpoint.

Same method/resources: wB97XD/def2TZVP, 16 physical cores, 32 GB. Numerical
remedy: explicit Cartesian ReadAllGIC plus frozen chi1, MaxStep=2. Same 12
frozen heavy atoms, 27 mobile atoms, no reactive-bond restraint. AF3-seed
scaffold/torsions remain distinct from n12c reference values. GIC toy test
preserved frozen XYZ and torsion; eight focused tests pass (87% focused-run
coverage, preparer 94%, assessor 93%). Gaussian accepted the full GIC input.

Next inspect actual convergence and chemistry; a separate exact Hessian plus
gradient/constrained curvature follows only after geometry passes. No Hessian,
RESP fit, accepted force field, protein refinement, MD or AF3 job was started.
Water belongs to later enzyme/QM/MM validation with standard TIP3P in MM.
See docs/native_ti1_reference_20261023.md. Action journal remains
output_agent_actions/20261021.jsonl; unrelated cif/ preserved. Earlier active
42/103-atom notes below are historical and superseded.

## Truncated 42-atom pilot replaces 103-atom pilot — 2026-10-23

The user explicitly authorized truncation and more cores on violet. The old
worker PID 3922102 was identity-checked and sent SIGTERM; its job ended with
exit status 1. Preserve the old native_frozen_scaffold_pilot_01 directory,
including user_requested_interruption.json. It did not complete; this is
supersession rather than chemical rejection.

ACTIVE tmux: `nyl12_truncated_ti1_16core_20261023`.
Directory: `output_nyl12_truncated_reference_20261023/truncated_pilot_01/`.
The model has 42 atoms, formula C10H25N3O4, 25 explicit H, net charge 0 and
singlet multiplicity. It retains Thr NH3+, the oxyanion and explicit water.
C54 is R; Thr CA/CB are S/R. Retain two substrate carbons on each side; cut
C70-C78 and C75-C79 and add H at 1.09 A along the old bond directions.
All 40 retained atom coordinates are unchanged.

Gaussian uses 16 physical cores (CPU affinity 0-15), 32 GB RAM and 677 basis
functions. Method: wB97XD/def2TZVP optimization only, no CalcFC/Freq, extra
dispersion or solvent. Same AF3-seed chi1 and frozen-scaffold pilot policy:
12 frozen heavy atoms, 30 mobile atoms, reactive bonds free. This does not
replace the n12c torsion reference. Eight focused tests pass: 97% overall
coverage, 94% preparer and 95% assessor. See docs/native_ti1_reference_20261023.md.

Next assess actual convergence, chemistry and water geometry, then consider
a separate exact Hessian, gradient and constrained-curvature assessment.
No native RESP, force-field transfer or protein relaxation is accepted.
No MD: the user's proposal remains explicitly deferred. Active journal:
output_agent_actions/20261021.jsonl. Unrelated cif/ is preserved. The earlier
103-atom running-job note below is historical and superseded.

## Native TI1 pilot started — 2026-10-23

User inspected capped geometries and authorized the next step. Native seed
mapping establishes **R** at C54 (Thr CA S, CB R); do not select S to evade
R cap sensitivity. A source-mapped103atom/60explicitH charge0 native seed
passes initial graph/stereo/contact screening. Water comes from the earlier
physical water-only TIP3P restart, not the manual displacement or failed DFTB
endpoint. See docs/native_ti1_reference_20261023.md for exact construction.

Active tmux: `nyl12_native_ti1_pilot_20261023`; directory
`output_nyl12_native_reference_20261023/native_frozen_scaffold_pilot_01/`.
Gaussian wB97XD/def2TZVP,8CPU32GB;65mobile atoms,38frozen heavy atoms plus
seed chi1 constraint. This is a frozen-scaffold state-retention pilot using
explicitly labeled AF3-seed torsions, not a replacement for n12c-derived
parameter-reference torsions. No reactive-bond restraint, solvent or extra
dispersion. Engine started; no convergence or scientific acceptance yet.
Six checker corruption tests passed;95% assessor coverage. Tiny Gaussian
freeze test confirmed fixed coordinates unchanged. Preserve all source hashes.

The user suggested eventual explicit-solvent MD to find relay-water snapshots
for QM/MM, and explicitly said **do not do this yet**. Proposal documented;
no MD or new AF3 run. Parameter-transfer/protein-relaxation gate stays closed.
Next inspect pilot process and scientific checks; do not auto-start RESP,
production MD, or repeat the old DFTB repair. Active journal remains
output_agent_actions/20261021.jsonl. Unrelated untracked cif/ preserved.

# Standalone pipeline handoff — 2026-10-20

## Final follow-up checkpoint — 2026-10-23

All four capped-surrogate HF/MK ESP datasets completed and passed separate
scientific acceptance (8606/8843/8712/8897 grid points). No RESP fit or native
transfer occurred. The signed raw-Hessian R/−60 OG1–C* projection is
−5.88965426 kcal/mol/A² and was rejected as an MM parameter; it does not
overturn the separate constrained-stability result. See the final section
of docs/ti1_qm_results_20261023.md. No further calculations were launched
after the user requested a quicker wrap-up. The next work is larger/native
reference preparation and validated parameter transfer, not protein repair.

## Completed Gaussian references and next gate — 2026-10-23

User resumed and requested results, documentation and next steps. All four
37atom capped TI1 optimizations, separate frequencies and conditional-curvature
assessments finished successfully. Each has105 positive raw frequencies and
102 positive constrained internal directions; no proton transfer or graph/
stereochemistry failure was detected. These are constrained surrogate minima,
not native enzyme TI1 or force-field validation. See
`docs/ti1_qm_results_20261023.md` for all four bond/frequency values and scope.

R cap sensitivity FAILS: OG1–C* shifts0.05327A and C*–O− shifts0.03166A,
above the predeclared0.03A threshold. S length sensitivity passes. Do not select
the passing configuration or average references to hide failure. Short capped
Thr N–H/oxyanion contacts persist. Larger/native reference geometry and native
charge/transfer validation remain necessary before protein TI1 refinement.

Next work: independent frozen-output audit/final inspection exports, diagnostic
raw-Hessian bond-stiffness comparison and separate capped-surrogate HF/MK ESP.
No native charge transfer, force-field fit or protein relaxation is certified.
New work is under `output_nyl12_parameter_validation_20261023/`; fresh surrogate
ESP uses separate `surrogate_esp_01` calculation directories when launched.
Read their process/acceptance records to distinguish preparation from execution.
The final standalone ESP checker passed27 tests/100%line coverage; full regression
passed447 tests with13 optional workstation-fixture skips, plus Ruff/format.
Final source-pinned verification: `output_nyl12_parameter_validation_20261023/full_verification_02/`.

The prior AF3-water discussion is now consolidated in
`docs/af3_water_placement_decision_20261023.md`: both minimal AF3 water tests
failed geometry acceptance; physical preconditioning removed crowding but the
subsequent DFTB3/MM diagnostic failed TI1 preservation. Better AF3 restraints
are useful parallel work but not the present prerequisite. Leave the user's
AF3 fork untouched; no new AF3 run or unvalidated protein refinement is started.

The active append-only journal remains `output_agent_actions/20261021.jsonl`.
Retain permission-review timeouts/bwrap failures and all failed attempts.

## Ranking and LigandMPNN symmetry clarification — 2026-10-22

New decision record: `docs/design_ranking_and_iteration_20261022.md`. The user
reports single-mutant confidence/strain degradation relative to WT, and expects
coordinated LigandMPNN redesign to improve fold/context compatibility. Preserve
the observation without treating confidence as affinity or all WT-relative
deformation as bad. Proposed ranking: chemistry gates, independent ES geometry
robustness and multistate/access compatibility, then local confidence/physical
quality/strain and noncovalent interface_delta as support. Sequence convergence
is diagnostic. No ranker, matched no-MSA arm or design campaign is implemented.
At the user's wrap-up clarification, primary geometry reproducibility refers to
unguided noncovalent ES. Covalent TI1 geometry and TS saddle validation are
separate gates. Interface_delta is optional pending demonstrated incremental
value/accuracy; it may be omitted. Leave existing QM tmux jobs undisturbed.

Official LigandMPNN README/run.py confirm automatic homooligomer mode ties all
parsed chains. Mature Nyl12 needs explicit separate alpha and beta ties with
equal weights and homo_oligomer0, as already used in the recent October16
protected preparation. Sequence ties do not select a geometric point group;
D2-like scaffold evidence and experimental assembly verification remain distinct.
Methods and the laptop inspection export are updated at this checkpoint;
scientific job state must still be read from current status/acceptance records.

## Redesign strategy/accessibility checkpoint — 2026-10-22

The user supplied the more recent `/data/jpc/enzyme-substrate-modeling/` source.
October16 focused redesign uses explicit separate alpha/beta ties, equal
weights and homo_oligomer0. Use its protected_inputs, not superseded root/pilot
snapshots; these were preparations, not an executed design campaign. Exact
source/provenance and caveats are recorded in the new strategy document.

The user's latest requirements and careful options comparison are recorded in
`docs/redesign_backbone_symmetry_and_access_20261022.md`. Recommended first
campaign: focused LigandMPNN on WT backbones/WT mature lengths, then independent
AF3 and relaxation. RFD3 is optional; if needed, initially remodel only selected
short pocket/entrance loops at exact WT lengths. Exact loop ranges are not yet
selected. Preserve catalytic/maturation/core/interface context. Indels remain
disabled until a per-design correspondence map and all postprocessing tests
cover them; WT-specific exporters currently hard-code232/96 and Amber indices.

The enzyme is (alpha-beta)4, not eight identical chains. Tie alpha copies and
beta copies independently. User has working LigandMPNN symmetry examples and
will point to them when needed. Do not assume C4 from tetramer stoichiometry:
a source-pinned tmux audit of the retained scaffold is D2-like (three near180deg
rotations, orthogonal axes, CA RMSD0.051–0.116A). Verify experimental assembly
operators before production RFD3. Raw audit is
output_nyl12_parameter_validation_20261022/repeat_symmetry_audit_01/.

Access screening must include open-compatible apo/ES and catalytic states,
actual substrate/product sterics, matched WT controls and state-aware sampling.
Sequence symmetry is distinct from identical dynamics/occupancy. AF3 confidence
and minimization do not establish opening kinetics; no opening in a short run
is unresolved. LigandMPNN multiple-PDB mode is batching, not joint ensemble
design. RFD3 exposure hints are static and do not guarantee channel access.
No RFD3/LigandMPNN campaign or new AF3 run has launched here. User AF3 fork is
untouched. Independent standard Amber typing/coverage preparation can proceed
without fitting charges or relaxing a protein with an unvalidated force field.

## Amber workflow clarification and parallel checks — 2026-10-22 evening

The user asked whether this uses standard Amber parameterization with
antechamber/Gaussian and ff19SB plus GAFF/GAFF2. Source audit confirms the
legacy scaffold uses original GAFF, not GAFF2, and provisional AM1-BCC or
Gasteiger charges. The new coupled TI1 Gaussian/RESP fit has not run.
Do not describe the existing topology as Gaussian-RESP parameterized.

New candidate plan: retain ff19SB for unmodified protein, use GAFF2 for
ordinary substrate chemistry, and fit separate HF/6-31G(d) MK ESP/RESP charges
with AmberTools utilities and explicit coupled-state/cap constraints. Use
antechamber for typing/file preparation and parmchk2 for parameter inventory;
validate their assignments. First test the standard candidate; refit only
missing or demonstrably inadequate reactive terms. Modified-Seminario is an
optional targeted tool, not an antechamber replacement or an automatic refit
of all bonds/angles. No legacy topology was migrated and no new FF was fitted.
Full detail: docs/ti1_parameter_protocol_20261022.md, final section.

While Hessians run, installed MCPB angle routines passed 96 manufactured
pair-spring/compliance checks plus three shared-angle tests, covering units,
AMBER harmonic convention, scaling, endpoint reversal and rotation. This is
API/convention verification, not recovery of an arbitrary physical molecular
Hessian or validation of TI1 parameters. Exact script, process log and report:
output_nyl12_parameter_validation_20261022/mcpb_angle_conventions_01/.
Compact versioned evidence: provenance/ti1_angle_conventions_20261022/.

The prospective independent validation plan is now
`docs/ti1_heldout_validation_plan_20261022.md`: 22 initial bond/angle distortions
per accepted reference, fixed baseline torsions, exact-coordinate QM/full-MM
comparisons and explicit training/validation separation. No scan inputs or
Gaussian scan jobs have been generated. Torsions and native-state transfer
still need separate coverage. Existing scientific gates and live tmux jobs
remain unchanged; user AF3 work remains untouched.

## Resumed checkpoint — 2026-10-22, 20:36 EDT

The user explicitly resumed; the previous 21:00 UTC stopping point is historical.
Two optimizations have now converged normally and passed the geometry gate:
`attack_S_cap_phi_minus60` and `attack_R_cap_phi_minus60`. Their fresh frequency
calculations are active in the existing tmux controllers. R/-120 and S/-120
remain unconverged at the time of inspection. No completed TI1 Hessian or
conditional-stability assessment has yet been accepted. Do not fit yet.

An independent frozen-log S/-60 audit reproduced its accepted geometry in
`output_nyl12_parameter_validation_20261022/completed_reference_audit_02/`.
Its `result/` contains assessment.json and final_structure.xyz/.sdf with all
37 atoms and 20 explicit H. Reactive bonds: OG1-C* 1.508577 A, C*-O- 1.302946 A,
C*-N 1.472865 A, C*-C 1.547196 A. No proton transfer detected. The capped Thr
backbone N-H contacts the oxyanion at 1.581748 A; retain this possible cap bias
for comparison against the other cap orientation. This is a constrained
bonded surrogate, not the native charged TI1; relay water is absent by design.

The installed pymsmt MCPB bond extractor passed eight manufactured quadratic
checks for bohr/Hartree units, one-based indexing, rotation, block averaging,
frequency-scale squaring and the AMBER K*(r-r0)^2 factor of one-half. Actual
TI1 fitting and modified-angle validation remain pending. Do not substitute
an arbitrary installed Seminario wrapper. Attempt01 emitted a complex-cast
reporting warning; attempt02 explicitly recorded real/imaginary components,
all with zero imaginary residuals. Both are retained. Independent geometry
audit01 failed before execution due to missing PYTHONPATH; fresh audit02 used
explicit project src/root paths and passed without changing chemistry gates.

Compact outputs, exact helpers and both failures are versioned under
`provenance/ti1_completed_reference_20261022/`. Production code/dependencies
are unchanged; the previous full-suite 433-test result still applies to that
source revision. New checks ran as scripted tmux jobs with source hashes.
Active action journal remains `output_agent_actions/20261021.jsonl`.

Next: finish the four geometry/Hessian/conditional-curvature gates, compare cap
hydrogen bonding and sensitivity, then fit and independently validate bonded
terms. Native 103-atom geometry and separate HF ESP/RESP/held-out validation
follow. AF3-model MM relaxation remains gated on the full TI1 force field,
then DFTB3/MM optimization. The user's AF3 fork remains untouched.

## Parallel preparation checkpoint — 2026-10-22, 20:36 UTC

The user visually inspected the four interim optimization trajectories and
judged them good. This is visual approval only: convergence, chemistry,
cap-artifact and native-transfer gates remain unchanged. Work resumed for
approximately one hour from 20:06 UTC; arrange a stopping point by 21:00 UTC.
Existing authorized Gaussian jobs continue in tmux during the user's absence.

At 20:35 UTC all four optimizers remain active and unconverged, around steps
11/11/11/12 in R/-120, R/-60, S/-120, S/-60 order. No completed TI1 Hessian,
ESP fit, accepted force field or new protein refinement is available.

Completed while waiting:

- Implemented conditional stability from the actual Gaussian gradient/Hessian,
  including the nonlinear constraint-force (Lagrangian) correction. Three
  frozen torsions plus six removed rigid motions leave 102 internal directions.
  Raw energy-Hessian projection or raw frequencies alone are insufficient.
  Both finite-difference steps must pass the prospective numerical criteria in
  docs/ti1_parameter_protocol_20261022.md. Missing gradients fail explicitly.
- Four tested analysis controllers now wait in detached tmux sessions
  `nyl12_ti1_curvature_<case>_20261022`. Their retained commands/statuses live
  under `<case>/conditional_curvature_controller_01/job/`, beside the existing
  frequency controllers. After successful frequency-controller completion,
  the pinned analyzer independently rechecks source/input identity, process
  acceptance, geometry and torsions, then writes a fresh
  `<case>/conditional_curvature_01/assessment.json`. Failed prerequisites,
  negative/unresolved curvature or nonstationarity cannot certify parameters.
  The source snapshot is `pa6_local_bonded/conditional_curvature_source/`.
  Waiters use a 48-hour deadline and 60-second polling, one BLAS thread, and
  launch no new Gaussian jobs. Original pinned frequency sources are unchanged.
- Prepared R/S native graph crosswalks (14 heavy-atom correspondences) and
  explicit proposed RESP constraints for the native 103-atom charge-zero
  Thr-NH3+/PA6-L4/oxyanion/NME/water graph. No native 3D geometry or ESP exists
  yet. No surrogate charge transfer or arbitrary residue integer charges.
- Independently sourced ff19SB/TIP3P with tleap and verified all 11 proposed
  fixed atoms against clean exported NTHR/NME/TP3 units. Peptide C/O are the
  NTHR values +0.6163/-0.5722. Native free-atom charge target is -0.0441 e.
  This is a source-consistent proposal, not an executed/accepted charge fit.
  Whole-library CHCL3BOX parser warnings remain in the first retained log;
  clean target-unit extraction/verification provides the source acceptance.

Verification: **433 full-suite tests passed**, Ruff and changed-file formatting
passed. New/modified module coverage is 84–96% (combined 89%). Separate waiter
safety tests passed 14/14, 94% coverage, without invoking engines. All checks
ran as scripted tmux jobs. Compact evidence, exact helpers and source hashes:
`provenance/ti1_parallel_preflight_20261022/`; larger originals are under
`output_nyl12_parameter_validation_20261022/`. The active append-only journal
remains `output_agent_actions/20261021.jsonl`. No dependency was added.

Next: inspect optimizer/frequency/conditional-curvature acceptance, then assess
cap H bonding and cap sensitivity before extracting bonded terms. Validate
bonded candidates with independent geometries/scans; build accepted native
charge-fit conformers and execute separate HF/6-31G(d) ESP/RESP/held-out checks.
Only after the full TI1 force-field gate passes may AF3-model classical repair
resume, followed by DFTB3/MM. The AF3 fork remains the user's independent work.

## Interim TI1 validation — 2026-10-22, 19:45 UTC

User requested validation and inspectable Gaussian outputs. All four optimizers
remain active around steps 6–7, with force criteria unmet; no frequency or final
parameter acceptance is available. Frozen live-log prefixes were audited in
tmux under `output_nyl12_parameter_validation_20261022/progress_audit_01/`.
Its `result/` contains 7/7/7/8-frame XYZ trajectories, latest XYZ/SDF, input
and atom maps, all-frame geometry checks, CSV metrics and PyMOL viewers.
Every frame has 37 atoms and 20 explicit H. Latest frames pass the declared
geometry/stereo/torsion/proton screens, but remain unconverged.

Both S cases show capped Thr backbone N–H/oxyanion contacts (1.77–1.93 A),
so cap H bonding needs explicit review before accepting any fitted reference.
The earlier S/-120 frame 4 had C*–O 1.225 A and failed the predeclared bond
screen; latest C*–O is 1.291 A. Retain all trials and unchanged final gates.
Do not interpret the latest-screen pass as convergence, a free minimum,
validated parameters or proof of negligible cap artifacts. Inspect S/-60
against R/-60 first. Relay water is intentionally absent from these small
bonded surrogates, not missing from the enzyme refinement model.

Actual PyMOL checks verified every exported state, 20 visible H, 36 intended
display bonds, coordinate agreement below 2.4e-7 A, and repeatable viewer reloads.
Checks: 371 full-suite tests passed; 34 focused export/converter tests passed,
97% packager coverage. Raw check outputs are `progress_viewer_check_03/` and
`inspection_full_verification_01/` under the same parameter-validation root.
The laptop package adds `qm/ti1_bonded_reference_progress/`; the existing
zero-argument pull script is unchanged. Original package revision 07 is retained
when the refreshed bundle is promoted. The next real calculation-stage gate
remains completed-optimization acceptance, then the existing automatic Hessian
continuation. Fitting/native-state transfer validation cannot use these interim
snapshots as final references.

## Active TI1 parameter checkpoint — 2026-10-22, afternoon

The user explicitly resumed after the earlier immediate stop. Four actual
Gaussian constrained optimizations are now running in detached tmux sessions
on violet, each with 8 cores and 16 GB. Inspect their retained status/logs;
this statement describes launch, not successful completion or fitted parameters.

The user selected `cif/n12c.cif` as the geometry reference. Author chain A,
Thr233 gives psi 134.110276428 degrees and chi1 -58.964417878 degrees;
native phi is undefined at this N-terminal Thr. The prescribed small bonded
surrogate is ACE-THR(-O-C*(O−)(CH3)(NHCH3))-NME, 37 atoms, charge -1,
singlet. Compare explicit ACE-cap phi choices -120/-60 degrees and both
attacked-carbon R/S configurations. Only phi, psi and chi1 are frozen;
reactive bonds and hydrogens are free. This capped anion does **not** supply
transferable charges for the native NH3+/oxyanion TI1.

All four deterministic starting geometries passed the predeclared geometry
screen. Prepared inputs, exact source/CIF snapshots, hashes and per-case jobs:
`output_nyl12_ti1_parameters_20261022/pa6_local_bonded/`.
Cases are `attack_R_cap_phi_minus120`, `attack_R_cap_phi_minus60`,
`attack_S_cap_phi_minus120`, `attack_S_cap_phi_minus60`; each actual run is
`optimization_01/`, containing `input.gjf`, `optimization.log`, retained
launcher and `job/status.json`. Tmux session names are
`nyl12_ti1_<case>_20261022`. Use the existing passwordless `sg g16 -c ...`
wrapper, explicit Gaussian profile and descriptive scratch directories.
No sudo or permission changes were made.

Optimization method is wB97XD/def2TZVP with UltraFine grid, Tight SCF/Opt,
CalcFC and XQC, no extra dispersion. Do not launch frequencies/ESP until
normal completion **and** TI1 protonation/connectivity/stereochemistry pass.
Four source-pinned frequency continuations also run in tmux as
`nyl12_ti1_freq_<case>_20261022`; their status/logs are under
`<case>/frequency_controller_01/job/`. They wait for the optimizer's final
process record, verify unchanged input hashes and apply the tested scientific
gate. Only an accepted reference proceeds to fresh `<case>/frequency_01/`,
with a separate Hessian checkpoint, formchk output and geometry cross-check.
Rejected references stop with retained acceptance records. There is no
automatic ESP fit, parameter certification or protein refinement.
Frequency and HF/6-31G(d) MK inputs use separate checkpoint names. Raw
frequencies of this constrained reference do not certify a free minimum;
conditional curvature, bonded fitting, held-out validation and native-state
RESP/transfer checks remain outstanding. The protein-refinement gate is closed.
See `docs/ti1_parameter_protocol_20261022.md` for the prospective acceptance
criteria and scope.

Latest full verification passed **368 tests**, Ruff and 87% combined coverage;
each new scientific module/entry point measured 82–99%. The separate
continuation safety tests passed **9 tests**, Ruff and 93% coverage without
invoking engines. Retained outputs are `full_verification_01/` and
`controller_verification_01/` under `output_nyl12_parameter_validation_20261022/`.
A real neutral-reference checkpoint conversion succeeded
in tmux under `output_nyl12_parameter_validation_20261022/real_fchk_probe_01/`;
the first parser check found a fixed-width character-array bug and was retained
as a failure. The corrected reader handles those records, then correctly
rejects this final-HF checkpoint because it has no Cartesian Hessian. No real TI1
Hessian has yet been accepted. Active append-only journal remains
`output_agent_actions/20261021.jsonl`; all original artifacts are preserved.

## Resumed 2026-10-22: viewer correction and TI1 validation

User is independently improving the AF3 contact-restraint fork and will report
when ready. Leave that fork untouched; no new AF3 tests until their update.
Water hydrogens exist in all 100 saved frames. Viewer scripts hid them, and a
semicolon in a comment caused the reported PyMOL SyntaxError. Corrected viewers
were tested twice in actual PyMOL 3.1.0: 100 trajectory states/one final state,
O/H1/H2 inventory, both O–H bonds, two visible water H, and no parser errors or
state duplication. Coordinate files are unchanged. Raw before/after reports and
PyMOL console logs are under `output_nyl12_viewer_fix_20261022/`.
Next scientific work remains coupled TI1 force-field validation. Source audit
confirms guessed junction terms, unfitted reactive charges and arbitrary O49
charge balancing; neutral free-NH2 Thr-NME QM artifacts are reusable provenance
but do not validate the covalent NH3+/oxyanion state.

The user clarified substrate termination on 2026-10-22: retain a real charged
PA6 terminus to support L2 release; a capped XA/MA fragment may be appropriate
for local copolyamide chemistry, recalling the Dishner JACS model. This concerns
substrate caps, separately from the peptide cap on Thr. See
`docs/substrate_termini_decision_20261022.md`. Preserve the existing full PA6 L4
reference and its terminal chemistry; source inspection confirms NH3+/COO−
termini. Do not silently cap it. The user subsequently supplied the original
capped copolyamide SMILES, preserved verbatim in the decision note. Scripted
RDKit/tmux audit confirms neutral acetyl-X-A|M-A-N-methylamide, C29H47N5O5;
slash cut C17/O36/N18 in zero-based original SMILES atom order. Original input,
script, report and job records are in
`output_nyl12_parameter_validation_20261022/substrate_termini_01/`.

Latest user clarification: PA6 L4 and L2 were both detected by MS, supporting
the hypothesis that soluble L4 can be a substrate yielding two L2 molecules.
Keep this first PA6 reference, but detection alone does not establish turnover.
Both charged termini represent a small oligomer; a long-polymer end model
would instead retain one real N- or C-terminus with a neutral distal boundary.
That is a separate future comparison. Only certain soluble low-molecular-weight
oligomers are detectable in the user's MS activity assay; do not infer the full
cleavage mechanism or total polymer conversion from that readout alone.
No coupled TI1 reference, QM calculation or AF3 refinement has been generated
in this audit. Next define its stabilization before generating QM fitting data.

## New scientific prerequisite — TI1 force-field validation

The user now requires validated TI1 parameters before classical AF3-model
relaxation, followed by DFTB3/MM. Read
`docs/ti1_forcefield_validation_gate_20261021.md` before any further refinement.
The current reactive topology is provisional and does not satisfy this gate.
The user suspects the poor initial geometry led the QM/MM structure toward a
reactant-like state; this is plausible but is not an accepted ES assignment.
The original failed run has been reproduced in a separate trajectory-only
replay:100 full-system frames, endpoint coordinate difference0.0 A.
No new model refinement was run after this prerequisite was established.

## Latest evening checkpoint — 2026-10-21

**Terminal update before disconnect:** attempt05 stopped after its first 100-cycle
QM/MM stage at 22:37 Eastern. No calculation from this work remains active.
Water contacts pass (maximum shortfall 0.029899 A), but OG1–C54 stretched to
2.195122 A, the directional water bridge failed, and RMS 0.63057 exceeded 0.5.
Later stages were withheld. This is a retained failed TI1 repair.
Diagnostic canonical-chain PDB:
`output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation/01_hydrogens_water/optimized.pdb`.
Next review state-preservation during early water relaxation and the DFTB3
Hamiltonian/starting geometry before defining another experiment. Do not rerun
unchanged or relax acceptance. The former running-state notes below are historical.

Read `docs/local_repair_checkpoint_20261021.md` first. The user authorized work
for about45 minutes after21:56 Eastern, with wrap-up around22:41. Active tmux
job: `nyl12_local_preconditioned_20261021`, attempt05 under
`output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/`.
Original AF3 oxygen is preserved during hydrogen construction. User-approved
TIP3P water-only pre-minimization completed with all nonwater atoms fixed exactly;
water O moved1.408215 A by minimization. It exhausted200 cycles above its gradient
target and is labeled coordinate preconditioning. DFTB3/MM subsequently started
with converged SCC evaluations. Inspect actual job/stage status on resume;
no accepted local minimum is yet claimed. Earlier attempts and failures retained.

## Active work — 2026-10-21, after evening resume

User authorized AF3 rough-water placement -> explicit hydrogens -> restrained
local DFTB3/MM repair. Begin with the two-restraint alpha/beta PA6 L4 model.
Both AF3 predictions completed; both independent water screens failed.
Attempt_01 has a 1.673 A water–C62 clash despite satisfying N/OG1 bounds;
attempt_02 with 28 additional exclusions breaks the N bridge (4.847 A).
See `docs/relay_water_repair_20261021.md` for the newly authorized separate
repair stage and water naming/Lennard–Jones requirements. Preparation is in
progress; no local minimization has been launched at this writing.
The earlier 18:00 pause was followed by explicit user authorization to resume.
Active journal remains `output_agent_actions/20261021.jsonl`.


## Latest checkpoint — 2026-10-21, before 18:00 Eastern

- User authorized work until 18:00 Eastern and a minimal AF3 water test that may
  run after disconnect: one alpha A, beta E, covalent PA6 L4 TI1 substrate I,
  and explicit HOH ligand J. Full-sequence Thr233 is E1.
- Launched `nyl12_af3_water_pa6_l4_20261021` in tmux on GPU0. Artifacts:
  `output_nyl12_af3_water_test_20261021/pa6_l4_alpha_beta/attempt_01/`.
  Check `job/status.json`, `job/job.log`, `feature_acceptance.json`, and
  `prediction/` on resume; do not infer success from session absence.
- AF3 inference also completed and wrote one predicted structure; the output
  water/covalent geometry has NOT yet been independently audited.
- Actual feature preflight PASSED: one water O and both contact restraints
  retained. Explicit HOH ligand has a non-polymer chain type. Restraints connect
  E1 N and OG1 to J1 O, range 2.6–3.2 Å, water endpoint mobile. One seed
  20261021, one sample, three recycles; guidance 0.5/interval1/20steps.
  Existing MSAs/templates/CCD reused. Source/input snapshots and hashes retained.
  Heavy-atom feature success does NOT establish NH3 protonation, water H-bond
  directions, TI1 geometry, or catalytic accuracy. No QM/MM refinement launched.
- Scripted local relay repair remains experimental and has NOT run on real
  models. Finish search-path tests/review before launch. Geometry/helper and
  CLI/export tests exist; 249 project tests passed at this checkpoint before
  final newly added tests. Full new-code coverage is not yet established.
- Canonical inspection labels are `pa6_l4` and `maxa_maxa`; original identifiers
  remain provenance aliases. Viewer source fixes `resi 113+186` and adds
  `set valence, 0`. Existing protein coordinates remain unrefined/clashy until
  a separately accepted replacement exists. Do not describe renaming as repair.
- Active journal: `output_agent_actions/20261021.jsonl`. Repeated automatic
  approval-review timeouts delayed tools; no safety rejection was inferred.
  See live journal for final export/commit/push outcomes.

## Resume checkpoint — 2026-10-21

This section supersedes historical job-state statements below.

- Both overnight neutral free-NH2 L-Thr-NME Gaussian calculations finished
  normally and passed independent rechecking: 21 atoms, CA S/CB R, 57 positive
  frequencies, minimum frequencies 45.9612 and 42.2567 cm^-1. Opt/Freq used
  wB97XD/def2TZVP; separate HF/6-31G(d) ESP grids contain 5636 and 5767 points.
  They are QM prerequisites, **not fitted force-field parameters**.
- The read-only paired audit and optimized SDF/XYZ exports are under
  `output_nyl12_resume_20261021/qm_conformer_audit_v2/`. Conformer 02 is
  4.6017 kcal/mol higher in electronic DFT energy; same-order heavy-atom RMSD
  is 1.3993 A. No thermal populations are inferred. The first audit attempt is
  retained: it failed on repeated optimization completion during Opt/Freq;
  parser correction and regression tests preceded the separate v2 attempt.
- `scripts/gaussian_log_to_xyz.py` converts logs without third-party Python
  packages. Both current logs were converted in tmux under
  `output_nyl12_resume_20261021/qm_trajectories/`: 64 and 24 Input-orientation
  snapshots, respectively, each with 21 atoms. Repeated frequency/ESP frames
  remain; these are not MD trajectories or one-frame-per-iteration promises.
- The finite-grid relay diagnostic ran in tmux under
  `output_nyl12_resume_20261021/relay_feasibility/`. Neither original model has
  an eligible oxygen site on the 5832-point grid (best maximum shortfall
  0.2510 A PA6 / 0.2679 A poly; gate 0.10 A). Optimistically excluding only
  Asn186 side-chain atoms and Thr929 N hydrogens exposes 9 / 5 eligible sites,
  respectively. No atoms were actually moved. This supports trying the
  documented rotamer/NH3 orientation repair; it does not prove a valid water
  orientation, a clash-free rotamer or an energy minimum.
- RESP still requires an explicit shared multi-conformer fit, cap constraints
  and transfer validation. Atom names repeat across THR/NME: use residue plus
  name or mapped indices. The current nonreacting beta-N charge -0.6254 is a
  charge-balance placeholder, not fitted evidence. Do not change only that
  charge or transfer neutral-template charges onto TI1 NH3+.
- User requested descriptive cache directories, a small laptop rsync package,
  a Gaussian trajectory converter, complete action logging, and commit/push.
  See `docs/laptop_inspection_20261021.md` and
  `docs/action_logging_and_rollback.md`. Logs and original scientific artifacts
  are preserved. Cache identity hashes remain metadata, not normal directory
  names. The inspection package distinguishes optimized small-molecule QM
  models from the still-unrefined protein TI1 starting models.

Next scientific work: explicit Asn186 rotamer/NH3 orientation sampling with
full-environment clash checks, followed by unchanged TI1 preflight before any
Sander launch; independently finalize RESP cap/transfer constraints. No repaired
TI1 minimum, TS1, ES model, reaction path or fitted FF is claimed by this resume.

Verification at this checkpoint: **212 tests pass**; Ruff lint, changed-file
formatting, shell syntax and diff checks pass. New measured modules have
80–97% coverage (89% combined); see
`provenance/resume_20261021_coverage.json`. Both cache directories were migrated
with all file and manifest digests unchanged; exact old/new paths are in
`provenance/cache_directory_migration_20261021.json`. Final code and documentation
are to be committed/pushed to the existing `codex/initial-import` branch at the
user's explicit request; the live action journal records the actual outcome.

## Evening stopping point — 2026-10-20, America/New_York

User requested a wrap-up for the night and authorized commit/push. Gaussian
jobs from the next section remain running independently in tmux on violet;
the laptop can disconnect. No new scientific jobs were launched in this
documentation/TS-capability assessment. Inspect current job records on resume.

- Started `manuscript/INTRODUCTION_DRAFT.md`, `METHODS_DRAFT.md`, `REFERENCES.md`
  and `SOURCE_LEDGER.md`. Existing Nyl12 manuscript sources were reviewed and
  left untouched. Completed preparation and proposed methods are separated;
  no TS, fitted parameter set or improved enzyme is claimed. Source hashes
  and revisions are in `provenance/manuscript_sources_20261020.json`.
- Read `docs/ts_search_and_symmetry_decisions_20261020.md` before TS work.
  Installed Sander supports native DL-FIND via `imin=1,ntmin=5`, with P-RFO,
  two-point numerical Hessians, Bofill updates and Hessian export. Source and
  library linkage were inspected; no runtime saddle-search test has been run.
  The interface follows the lowest mode and disables microiterations. Require
  a fresh final Hessian and endpoint tests in the declared mobile subspace.
- DL-FIND is freely available under LGPL; the user accepts it on that basis
  and requested newer alternatives. Sella is the preferred newer benchmark
  candidate (iterative Hessian diagonalization/RS-P-RFO); pysisyphus and
  geomeTRIC were also assessed. No new dependency was added or optimizer
  silently substituted. A benchmark proposal is documented, not executed.
- AF3 contact guidance already exists in
  `/data/jpc/code/alphafold3-contact-restraints` at
  `ef5dc91d65b4845f63e95ea03316f1613c6dfb5e`. Tests/source and actual historical
  production manifests/geometry summaries were located. Reuse this capability
  if needed; guided and independent prediction arms must remain distinguishable.
- Unrestrained ES ligands may occupy any equivalent site AE/BF/CG/DH.
  Score all four using a fixed register, permit symmetry-equivalent matches,
  and treat sites as nested within each model. The historical four-ligand
  analyzer's fixed I/J/K/L assignments and old atom names must not be applied
  blindly to single-ligand ES or current TI1 models. Symmetry-aware ES scoring
  remains to be implemented.

Immediate restart: check both running Gaussian jobs and acceptance records,
then continue parameter validation and TI1 geometry repair. Do not restart
duplicate calculations or claim the old water clashes have been resolved.

## Active Gaussian jobs and design clarification — 2026-10-20, evening

This section supersedes the job-state statements in earlier checkpoints.
See [method, substrate audit and multistate design plan](qm_method_and_substrate_checkpoint_20261020.md).

- User approved **wB97XD/def2TZVP** geometry/frequency calculations, with no
  additional dispersion correction, and separate HF/6-31G(d) ESP for RESP.
  Two neutral free-NH2 L-Thr-NME conformers (CA S, CB R, charge 0/singlet)
  are running on violet in tmux: `nyl12_thr_wb97xd_01_20261020` and
  `nyl12_thr_wb97xd_02_20261020`, each four cores/4 GB. Gaussian has produced
  SCF output; optimization/frequency/ESP acceptance is not yet established.
  Output root: `output_nyl12_resume_20261020/neutral_thr_nme_wb97xd/`.
  Check `jobs/conformer_*/status.json`, `calculations/conformer_*/gaussian.log`
  and eventual `process.json`/`acceptance.json`. Do not rerun blindly.
- Preparation, source snapshots, input/tool hashes, software versions, logs,
  separate scratch and content-addressed cache are implemented. Accepted QM
  data are not fitted FF parameters. RESP/cap constraints, bonded/torsional
  validation and production use remain unfinished. Existing unrelated tmux
  sessions were left untouched; no other workstation was needed.
- Completed ligand audit in tmux:
  `output_nyl12_resume_20261020/confidence_audit_v2/report.json`.
  Literal sequence is **MAXA|MAXA**, producing two MAXA fragments. Stable target
  ID remains `poly_xama_maxa`. Ligand mean pLDDT is 29.64 (59/69 atoms <50),
  terminal monads mean 13.36; reactive C/O/N 83.84/54.43/60.16. PA6 ligand mean
  60.84; reactive C/O/N 92.53/85.86/87.64. Retain full MAXA|MAXA as the
  product-consistent reference while sampling uncertain tails; no truncation
  or coordinate change occurred. Initial audit failed on rounded JSON versus
  full-precision CSV ranking scores and is preserved under `confidence_audit/`.
- Completed both descriptive metric baselines in tmux under
  `output_nyl12_resume_20261020/active_site_baseline/`:
  PA6 (d1,d2,d3,a_BD,phi) = (1.3320,2.7468,2.9584,101.9127,-99.1807);
  poly = (1.2864,2.9455,2.9865,103.3262,87.9888), distances A/angles degrees.
  Signed projected Herschlag convention matches the original reference math.
  These are measurements of original unrefined starts, not acceptance results.
- Recommended downstream plan: validated TI1 -> TS1 -> ES; focused LigandMPNN
  initially, optional RFD3 local backbone diversification; independent covalent
  AF3 and unforced noncovalent ES predictions; cheap multistate checks including
  acyl enzyme/water and TI2 geometry; reaction barriers for a shortlist.
  Map the reference cycle coarsely so TS1 is not assumed rate-limiting. This is
  a documented plan, not completed redesign or a reaction profile.
- All **149 tests pass**, Ruff passes. New-code coverage: geometry 97%, metric
  CLI 95%, template preparation 94%, Gaussian acceptance 100%, Gaussian runner
  98%, confidence audit 82%. Evidence:
  `provenance/qm_preparation_20261020_coverage.json` and
  `provenance/qm_design_checkpoint_20261020.json`. Routine scripting delegated
  to GPT-6 Luna and reviewed. No dependencies added; changes remain uncommitted.

TI1 water-contact failures remain unresolved. No new Sander optimization,
accepted TI1 minimum, transition state, ES model or validated parameter set
has been produced. Next: inspect the running Gaussian jobs without duplicating
them; fit/validate parameters after scientific acceptance; implement local
TI1 geometry repair before starting reaction paths or redesign.

## Resumed checkpoint — 2026-10-20, 17:50 America/New_York

This section supersedes the earlier pause/runner status below. The previous
bwrap-blocked resume attempt made no changes. This session also encountered
`bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted`; approved
escalated execution worked. No infrastructure settings were changed.

- Integrated the expanded audited mask into the TI1 input writer. Added
  `scripts/run_nyl12_ti1_qmmm.py`, `ti1_attempt.py`, and `ti1_checks.py` for
  explicit fresh attempts and terminal-checkpoint continuation. Legacy implicit
  build-and-run is rejected before writing. No dependencies were added.
- New attempts preserve actual inputs, code snapshots, environment/version
  metadata and hashes, original restraint references, and separate process and
  scientific acceptance records. Continuation rejects mismatched topology,
  reference, target/method, executable, DFTB data, or checkpoint digest.
- Both real-model preflights ran in tmux and have terminal wrapper exit 1:
  `/data/jpc/covalent-enzyme-design/output_nyl12_resume_20261020/jobs/preflight_pa6_l2/`
  and `jobs/preflight_poly_xama_maxa/`. Attempts are under the same output root
  at `{pa6_l2,poly_xama_maxa}/starting_model_preflight/`.
  Both passed expanded-QM and proton/connectivity inventories, then failed the
  original water-placement contact gate: maximum shortfalls 0.686619 and
  0.661263 A, respectively, against 0.10 A. These are the original LEaP seeds,
  not the improved but still failed rigid-search candidates from earlier jobs.
  Sander did not start; there is no new optimized geometry or restart.
- Routine CLI/test scripting was delegated to GPT-6 Luna and reviewed. All
  89 tests pass; Ruff lint passes. New-module line coverage: attempt 84%,
  chemistry/convergence checks 92%, CLI 96%. Full coverage evidence is
  `provenance/ti1_resume_20261020_coverage.json`. The new files are uncommitted.
- The concrete next work is local conformer/geometry repair as described in
  `docs/ti1_attempt_protocol_20261020.md`. No repair success is claimed.
  Gaussian/RESP fitting, separate neutral-NH2 ES, and parameter validation
  remain unfinished. Do not lower the water gate or repeat an unchanged job.
- A completed tmux export produced both full hydrogenated, canonical-chain
  starting models for PyMOL:
  `/data/jpc/covalent-enzyme-design/output_nyl12_resume_20261020/nyl12_pymol_starting_models.tar.gz`.
  It includes the viewer, residue maps, original audits, provenance, and hashes.
  The viewer now ends with the user-requested commands
  `show sticks, resi 113 +186 and chain A` and `util.cnc`.
  The original archive/viewer/manifest remain in `inspection_revisions/v1/`
  under that output root; the revision script and source export snapshot are
  in `inspection_source/`. No enzyme coordinates changed during this revision.

No jobs from this resumed checkpoint remain active. The two unrelated tmux
sessions listed in the historical pause note remain untouched. Code tests and
passing inventory checks do not certify physical models or MD readiness.


## Paused at the user's request — 2026-10-20

Ready to resume; no calculations from this work remain active. All tracked
job attempts have terminal exit records. Failed water searches and superseded
diagnostic attempts remain preserved alongside successful later diagnostics.
All delegated agents finished. Existing unrelated tmux sessions
`nyl12_effective_strain_20261015` and
`nyl12_focused_redesign_preflight_20261016` were left untouched.

Implementation checkpoint: `ab85811962f40f8988f9113c3c9e7b9edcc0e13c` on
`codex/initial-import`, pushed to the private remote. This pause update changes
documentation only; the last checks remain 51 passing tests and passing Ruff.
No new jobs were launched for the pause.

To resume, open `/data/jpc/covalent-enzyme-design` and read this note first.
Start with the **Next scientific work** section below. Do not rerun the old
TI1 runner unchanged: it still uses the earlier QM boundary and lacks safe
checkpoint continuation. Geometry repair and charge fitting remain unfinished.

Suggested resume request:

> Read `/data/jpc/covalent-enzyme-design/docs/standalone_handoff_20261020.md`
> and resume Nyl12 TI1 work. Delegate routine scripting to less expensive
> models, run scientific jobs in tmux, and preserve manuscript provenance.


## Authoritative locations

- New development checkout: `/data/jpc/covalent-enzyme-design`.
- Private remote: https://github.com/parksjm/covalent-enzyme-design .
- Original source and scientific artifacts: `/data/jpc/enzyme-substrate-modeling`.
- Current calculation artifacts: `/data/jpc/enzyme-substrate-modeling/output_nyl12_ti1_refinement_20261020`.
- Full scientific history: `docs/nyl12_covalent_resume_20261019.md` (historical relative output paths refer to the original workspace).

Continue code development in this standalone checkout. Use explicit absolute
input/output/tool paths while workstation-specific defaults are being removed.
Existing jobs and results were not relocated. The import namespace and CLI
names remain unchanged. The older workspace is retained for data and provenance.

## Extraction provenance

`provenance/source_inventory.json` and `provenance/import.json` record the
42-file source snapshot and archive digest. `provenance/import_snapshot.py`
verified the complete archive inventory and each file digest before writing.
The archive, original scripts, and tmux export/import records remain under
`output_nyl12_ti1_refinement_20261020/repository_export` and `jobs/repository_*`
in the original workspace. The imported inventory generator is a migration
utility with the original explicit allowlist, not a manifest of every future
file added to this repository. Git commits and the uv lock now version further
source changes.

Root README/AGENTS/pyproject/gitignore started from the recorded templates.
Post-import changes include the uv lock, pytest-cov as a development-only
dependency, and pytest's root import path for the scripts namespace. Production
dependencies remain NumPy 1.26.4, RDKit 2025.3.6, and ParmEd 4.3.0.

## Scientific status

Both TI1 topology inventories pass; neither geometry is an optimized minimum.
The rigid-pocket relay-water search has no accepted candidates. Contact
reports diagnose heavy-atom crowding as well as hydrogen contacts. Avoid
interpreting this finite rigid search as evidence against the reaction mechanism.

The expanded QM selection passed structural audits: 109/181 atoms, formal
charge 0, selected partial sum -0.0294, link-H electron counts 364/628, and a
single CA930–C930 cut. **The old QM/MM writer still needs this new selection
integrated.** This audit does not validate electrostatic embedding or DFTB3.

Gaussian 16 C.02 on violet passed a real water HF/STO-3G single-point preflight.
Use `/usr/local/g16`, `g16root=/usr/local`, and scratch
`/data/scratch/gaussian_scratch/jpc`. The current long-lived process needs
`sg g16` to acquire the already granted group; a fresh login should acquire it
normally. The preflight is not a charge-fitting result.

## Next scientific work

1. Integrate the audited QM selection and explicit new-attempt/checkpoint paths
   into the TI1 runner; require convergence and post-run chemistry checks.
2. Develop a documented local geometry-repair strategy for crowded relay-water
   sites, with subsequent contact and protonation checks. Retain failed seeds
   and distinguish restrained starting models from stationary points.
3. Prepare a documented Gaussian/RESP neutral beta-N template protocol and
   cache validated artifacts by chemical state, conformer, method, and versions.
4. Build/audit the requested neutral-NH2 ES state separately; TI1 has NH3+ at
   reacting Thr233 and an intact substrate C–N bond.
5. Complete full bonded/charge parameter validation before classical MD.
   Reaction paths and enzyme design follow reliable TI1 structures.

The measured-product gates remain PA6 L2 and poly(X,AMA) MAXA. New covalent AF3
models preserve custom atom names; no affinity calculation or naming increment
is part of this workflow.

## Skills

`scientific-tmux-jobs` is installed at
`/home/jpc/.codex-hosts/violet/skills/scientific-tmux-jobs` and versioned under
`tools/skills/`. Its installed pipeline adapter already points to this checkout.
The covalent-modeling skill is authorized but awaits a validated refinement
workflow. Consider parameter-fitting/cache, canonical inspection, and manuscript
provenance skills when repeated use justifies separate routing.

## Final import checks

51 tests pass in the standalone uv environment; Ruff passes. Overall source
line coverage is 38%, with older chemistry/topology/prediction code untested.
New workflow modules/wrappers each have at least 80% line coverage. This
coverage does not validate the physical models. The installed tmux skill now
points to this checkout. Full machine-readable coverage, source differences,
and external artifact hashes are retained under `provenance/`.

## Model-switch handoff — 2026-10-21, evening

User wants GPT-6.1 Sol High for continuation. The rough-water restrained local
minimization remains authorized; no repair minimization has launched. Preserve
all current uncommitted work. New draft modules are `local_repair_model.py`
(explicit octamer topology subset -> 330-residue alpha/beta model; AF3 heavy
coordinate mapping; retained TIP3P parameters), `local_water_seed.py` (bounded
1.5 A oxygen displacement and water H orientation), and `local_hydrogens.py`
(local bonded-frame transfer from original topology). All under
`src/nyl12_covalent/`. These are development code, not accepted models.

Luna hydrogen helper agent added `tests/test_local_hydrogens.py`: initial run
reported 93% helper coverage but two incorrect expected-coordinate assertions;
agent corrected them, then could not rerun because tool approval timed out.
Water-seed test agent could not execute any command; no tests were added.
Root has not yet reviewed or tested these modules. Source/import integration,
terminal OXT handling, exact atom mapping, covalent inventory, QM boundary and
water parameters must be checked before use. No new dependencies were added.

Next: finish and test preparation, implement a separate local DFTB3/MM repair
runner and geometry audit, then run in tmux. Proposed stages: initial H/water
relaxation with protein/substrate heavy atoms restrained; local pocket relaxation
with outer protein restraints; then a separately recorded bridge-restraint
release. Use only two broad N/OG1–water distance restraints. Record flexible
selection and constants explicitly. These stages are not implemented yet.
Dimer QM mask is `:233,329,330|(:234&!@C,O)`, charge0, expected109 atoms,
sole link CA234–C234; verify before launch. Existing runner/audits use octamer
indices and must not be applied unchanged to this dimer. Do not silently relax
legacy scientific acceptance to force a launch.

Sources: topology/coordinates in
`/data/jpc/enzyme-substrate-modeling/output_nyl12_covalent_qmmm_20261019/pa6_l2_ti1_af3_neutral_beta/ti1.prmtop`
and sibling `ti1.inpcrd`; AF3 source is
`output_nyl12_af3_water_test_20261021/pa6_l4_alpha_beta/attempt_01/prediction/nyl12_pa6_l4_ti1_alpha_beta_water_contact/nyl12_pa6_l4_ti1_alpha_beta_water_contact_model.cif`.
Extract original residues1–232,929–1024,1313–1314 without changing parameters.
Sander: `/data/jpc/code/ambertools26/venv/bin/sander`.

AF3 audits are complete under
`output_nyl12_af3_water_audit_20261021/pa6_l4_alpha_beta/attempt_{01,02}/result/`.
Both failed water acceptance. Audit/export and optional exclusion-input code
and tests remain uncommitted from prior work. Laptop export has not yet been
refreshed with these diagnostic CIFs. Commit/push remains authorized at a
verified stopping point; do not claim it was completed.

Frequent automatic approval-review timeouts precede process creation; they are
not safety rejections. Default sandbox still fails bwrap RTM_NEWADDR. Record
failures honestly in `output_agent_actions/20261021.jsonl`; repair logging before
mutations. Neither unrelated tmux session was changed. Current tools may work
on a short single retry; do not repeatedly escalate unchanged rejected actions.
