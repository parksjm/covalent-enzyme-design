# Methods — working draft

Status: 2026-10-25. Sections explicitly distinguish completed preparation from
planned methods. Unresolved numerical settings must be filled from accepted
run manifests, not inferred from a generic protocol. References are in
[REFERENCES.md](REFERENCES.md); local evidence is in [SOURCE_LEDGER.md](SOURCE_LEDGER.md).

## Molecular identities and assembly — implemented preparation

Nyl12 was represented as a mature assembly of four alpha chains (A–D) and
four beta chains (E–H), with catalytic sites AE, BF, CG and DH [S1,S4]. The
catalytic threonine is full-sequence residue 233, corresponding to residue 1
of a processed beta chain. The current AE topology places this residue at
one-based Amber residue 929. Atom and residue correspondence files accompany
canonicalized models; identifiers from different prediction pipelines are not
interchanged without checking the underlying molecular graph.

Substrates were specified by explicit molecular graphs and intended hydrolysis
products. The PA6 substrate is L4 (`pa6_l4`), cleaved in an L2-producing
register; `pa6_l2` is a historical input alias. The copolyamide substrate is
MAXA|MAXA (`maxa_maxa`), historically labeled `poly_xama_maxa`. Its literal
modeled sequence is MAXAMAXA. M denotes hexamethylenediamine, X denotes
p-xylylenediamine and A denotes adipic acid [R1]. Graph cleavage yields two
MAXA fragments, each C26H42N4O5. This soluble fragment represents a selected
poly(X,AMA) register; it is not evidence for the activity of Nyl12 on separately
prepared poly(MAXA). Substrate length, terminal groups and protonation are
recorded independently of the polymer family label.

The intended first tetrahedral intermediate (TI1) retains the substrate C–N
bond and contains a new Thr OG1–scissile-carbon bond, an oxyanion, and a
protonated catalytic beta-chain N terminus. Its precursor ES has a neutral
beta-chain NH2 group and Thr OG1H. A neutral relay water was included in the
earlier working mechanism. The current coordinate-driving pilot instead tests
an Asp272-mediated relay without a designated QM relay water; bulk MM solvent
is retained. These state definitions are audited hypotheses; predictor
coordinates alone do not establish proton transfers or intermediate stability.

## Structural preparation and confidence audit — executed

The existing AF3 TI1 predictions were canonicalized with retained atom maps
and used to prepare topology/coordinate pairs. Topology audits check covalent
connectivity, formal/partial charge inventories, explicit hydrogens and QM
boundary membership. The expanded proposed QM selection is
`:929,1313,1314|(:930&!@C,O)`, including the catalytic residue, ligand, relay
water and part of the following residue; its formal charge is zero. Source
topologies, force-field settings and the link boundary are retained in the
calculation records. These checks do not certify parameter quality.

Original AF3 atom-level pLDDT and directional ligand–protein PAE were analyzed
separately from whole-complex ranking scores. Polyamide monad assignments were
derived from the molecular graph. Full-precision CSV ranking scores and rounded
summary-JSON values were retained, with an explicit two-decimal consistency
check. Confidence measurements characterize prediction uncertainty, not binding
free energy. The audit and original failed precision-check attempt are preserved
[S6]. No substrate truncation was performed.

Descriptive geometry comprises Thr OG1–scissile C (d1), Phe113 N–substrate O
(d2), Asn186 ND2–substrate O (d3), angle OG1–C–O (a_BD), and the existing signed
projected Herschlag phi convention [S4,S7]. For phi, carbonyl C is projected
onto the plane through substrate O, amide N and adjacent C; the signed dihedral
is evaluated as OG1–O–projected-C–adjacent-C. Degenerate constructions fail
explicitly. ES and TI1 require different geometry interpretation; ES near-attack
thresholds are not automatically applied to the tetrahedral state.

Both original TI1 starts failed the relay-water contact gate after passing
connectivity and QM-inventory checks. Sander did not run in these preflight
attempts. No accepted TI1 minimum or TS1 is available at this checkpoint [S5].

## Force-field parameterization — reference preparation executed; fitting planned

### Force-field components and AmberTools preparation

The candidate classical model will retain ff19SB for unmodified protein
residues [R8] and use GAFF2 for ordinary substrate chemistry [R15,R18].
Antechamber and parmchk2 in AmberTools26 will assign and audit substrate atom
types and identify missing or estimated parameters before assembly in LEaP
[R16,R17]. The installed GAFF2 parameter file identifies version 2.2.30
(October 2025); its digest and the loaded leaprc/library files will accompany
each generated topology. The original GAFF paper [R15] documents the parent
force field; the specific GAFF2 release is identified by the Amber
reference documentation and installed parameter file [R18]. Standard protein
terms will be retained rather than assigning GAFF2 types to the entire enzyme.
The modified catalytic Thr and its substrate linkage will be represented by
explicit nonstandard residue definitions with audited atom correspondence.

Relay water will retain TIP3P atom names, charges and Lennard-Jones parameters
when represented classically [R19]. Its explicit hydrogens and correspondence
to the water included in the proposed QM subsystem will be preserved. A water
included in the QM region is treated by the QM/MM Hamiltonian rather than by
adding a duplicate MM electrostatic contribution. Water/ion choices and the
combined model will require validation; ff19SB/GAFF2 selection alone does not
establish accuracy for this covalent active site.

The existing development topology uses original GAFF and provisional
AM1-BCC or Gasteiger charges. It has not undergone the coupled-state Gaussian
RESP fit described below and is excluded from validated TI1 refinement.

### Chemical-state-specific electrostatic parameters

The native TI1 charge model will retain the catalytic Thr NH3+ group, the
substrate oxyanion and intact scissile C–N bond, full PA6 L4 with its NH3+/COO−
termini, the peptide-boundary cap and a separate neutral relay water. This
scripted cluster has 103 atoms, total formal charge zero and singlet
multiplicity. It represents the soluble L4 substrate hypothesis; it does not
represent an indefinitely long polymer or establish experimental L4 turnover.
The 37-atom, charge-minus-one ACE/Thr/NME bonded surrogate described below
has different terminal chemistry and will not supply native-state charges.

Gaussian 16 revision C.02 [R20] will provide HF/6-31G(d) Merz–Kollman ESP data
at the separately obtained DFT reference geometries [R21,R22]. Native ESP calculations
await accepted native conformers; completed neutral free-Thr ESP calculations
are separate development references. AmberTools antechamber/espgen/respgen/resp
utilities will be used for atom mapping, ESP preparation and RESP fitting
[R7,R16,R17]. Explicit control files will implement a shared multi-conformer
charge set, fixed boundary charges and declared equivalences; a default
single-ligand antechamber invocation is insufficient to impose this model.

The initial boundary proposal preserves the installed N-terminal Thr peptide
C/O charges and the NME cap charges, while fitting the reactive Thr ammonium/
side chain jointly with the substrate. TIP3P water charges are provisionally
fixed; constrained and held-out ESP errors will determine whether that
boundary treatment is acceptable. These source values were independently
checked against clean LEaP exports. Exact total charge will be enforced,
without arbitrary integer partial-charge constraints on the two covalently
joined residues. Equivalences, conformer identities, RESP controls, fitted
charges and fit/held-out residuals will be archived. ES and other intermediates
will require their own proton inventories and parameter eligibility checks.

### Bonded parameters and independent validation

Standard ff19SB/GAFF2 assignments will be evaluated first. Missing or
inadequate reactive-region bond, angle, torsion or improper terms will require
targeted QM-based refinement; parmchk2 estimates will not be accepted solely
because a complete topology can be written. If needed, modified-Seminario
projection [R14] will provide candidate bond/angle constants from an accepted
DFT Hessian. The implementation's Amber harmonic conventions,
E_bond = K_b(r-r0)^2 and E_angle = K_a(theta-theta0)^2 with angles in radians,
will be checked explicitly. Such projection does not by itself determine
reliable proper-torsion or improper parameters or validate transfer to the
native enzyme.

Reference acceptance will require normal completion, retained atom order and
stereochemistry, intended TI1 bonds/protons, agreement with frozen torsions,
and local stability under those constraints. Nonlinear constraint-force terms
will be included in the reduced Lagrangian-Hessian assessment; raw constrained
frequency signs alone do not establish a free minimum. Cap hydrogen bonding
and sensitivity to the two declared ACE orientations will be assessed before
using the reference for parameter estimation.

A prospective initial validation panel comprises 22 independent distortions
per accepted reference: four reactive bonds at +/-0.02 angstrom, the six
angles at the tetrahedral carbon at +/-3 degrees, and the CB–OG1–C* angle at
+/-3 degrees. These energy points are reserved for validation rather than
Hessian-based fitting. DFT and full surrogate MM single-point relative
energies will be evaluated at identical coordinates with the reference's
frozen torsions preserved. The surrogate will use its own independently
fitted charges; native charges will not be transplanted onto the capped
anion. Distortions changing connectivity, protonation or stereochemistry will
be reported as outside the declared TI1 basin rather than silently discarded.
Separate torsional scans and native-state transfer tests remain necessary.

Project acceptance criteria are bond errors <=0.03 angstrom, angle RMS error
<=3 degrees and maximum angle error <=5 degrees, and held-out local relative
energy RMSE <=1 kcal/mol. Cap comparisons require local bond targets within
0.03 angstrom and candidate force constants within 20%; excessive sensitivity
requires reconsidering the reference model. These are prospective project
criteria, not universal literature bounds. Numerical errors by coordinate
family, complete parameter coverage, topology identity and native-state charge
validation must pass before AF3-model classical relaxation and subsequent
DFTB3/MM optimization. At this checkpoint no new coupled TI1 RESP fit, GAFF2
topology or complete force field has been accepted.

## Parameter-model quantum calculations — executed; RESP pending

### Capped anionic TI1 bonded reference — staged QM preparation and validation

A local bonded-parameter surrogate was constructed as
ACE-THR(-O-C*(O−)(CH3)(NHCH3))-NME, retaining the tetrahedral carbon and intact
leaving-group C*–N bond. Its explicit-hydrogen graph contains 37 atoms
(C10H20N3O4−), charge -1 and singlet multiplicity; Thr CA/CB configurations
are S/R. Methyl groups truncate the reacting substrate amide. The surrogate
has neutral peptide caps and does not represent the native catalytic NH3+
group, full substrate termini or relay water. Its charges cannot be transferred
directly to the native TI1 model.

The user-selected protein reference `cif/n12c.cif` comes from a coauthored Nyl12
manuscript currently in revision (user communication, 24 October 2026;
formal citation pending publication). It was parsed with an
explicit author-chain-A, Thr233 selection. Measured psi and chi1 were
134.110276428 and -58.964417878 degrees. Because native phi is undefined at
this chain terminus, ACE-cap phi values -120 and -60 degrees were declared
model choices. Each was combined with an explicit R or S attacked-carbon
configuration, yielding four references without assigning a preferred attack
face. ETKDGv3 coordinates, seed 20261022, were adjusted to these three
dihedrals and passed a conservative starting connectivity/contact screen;
embedding was not an energy minimization.

Four gas-phase Gaussian 16 wB97XD/def2TZVP optimizations were launched with
UltraFine integration, Tight SCF/optimization, CalcFC, XQC, no symmetry and
no additional dispersion [R4,R5]. Each uses 8 cores and 16 GB. Only phi, psi
and chi1 are frozen; reaction-center bonds and hydrogens remain mobile.
Completion, electronic state, atom order, stereochemistry, reactive-bond
lengths, proton inventory and frozen torsions must pass before a separate
DFT frequency/Hessian calculation. Independent checkpoint names preserve
optimization, Hessian and HF/6-31G(d) MK ESP data separately.

Constrained references are not free minima. Gradients, constraints and raw
Hessians must support conditional-curvature assessment before accepting
bonded parameters. Existing AmberTools modified-Seminario routines are being
evaluated for candidate bond/angle terms [R14]; no fitted terms are accepted
at this checkpoint. Cap sensitivity, held-out QM/MM comparisons and a separate
native NH3+/oxyanion RESP/transfer validation remain required. Inputs, source
and CIF hashes, preparation screens and optimization process records are under
`output_nyl12_ti1_parameters_20261022/pa6_local_bonded/`; prospective acceptance
criteria are in `docs/ti1_parameter_protocol_20261022.md`.

By 2026-10-23 all four optimization/frequency references passed their
predeclared geometry and constrained-stability assessments. Each raw report
contained105 positive frequencies; the corrected constrained Lagrangian
Hessian was positive in all102 internal directions at both derivative steps.
However, the R cap pair changed OG1–C* and C*–O− lengths by0.05327 and0.03166
angstrom, exceeding the prospective0.03angstrom cap-sensitivity limit. The S
pair passed the bond-length comparison. Short backbone N–H/oxyanion contacts
persisted. Consequently, these results do not qualify parameters for native
transfer; larger/native references and independent validation remain required.
Separate capped-surrogate HF/MK ESP and raw-Hessian bond-stiffness diagnostics
are preparation for surrogate evaluation, not a native charge fit. Full
results and source paths are in `docs/ti1_qm_results_20261023.md`.

An interim audit froze complete prefixes of the four active optimization logs
at approximately 19:45 UTC on 2026-10-22. All complete Input-orientation tables
were retained, yielding 7/7/7/8 snapshots of 37 atoms with 20 explicit hydrogens.
These include trial geometries and do not establish convergence. Latest frames
passed the predeclared geometry screens; both S references displayed capped-Thr
backbone N–H/oxyanion contacts of approximately 1.77–1.93 angstrom. Such contacts
require cap-sensitivity review before fitting. A transient S/-120 frame had
C*–O 1.225 angstrom and failed its bond screen; the failed frame was retained.
PyMOL checks verified atom/hydrogen inventories and coordinate fidelity in all
exported states. These observations are interim model diagnostics, not accepted
parameters or mechanistic results.

### Exploratory local relay-water refinement — failed diagnostic retained

The minimal alpha/beta PA6 L4 covalent AF3 model used explicit HOH and two
Thr233 N/OG1–water-O bounds(2.6–3.2 A). Both distances were satisfied but an
independent audit found a1.673 A water–substrate-carbon contact. Adding28
exclusion restraints lost the N bridge. The active preparation retains the
original AF3 oxygen, adds hydrogens preserving TI1 protonation/connectivity,
and retains TIP3P WAT/OW/HW parameters. Direct DFTB3/MM minimization from this
seed encountered SCC nonconvergence and was stopped. The approved alternative
uses fixed-environment, water-only Amber/TIP3P pre-minimization, followed by
staged local DFTB3/MM refinement with water in the QM region. Classical TI1
energetics are not validated by this coordinate-preparation step. Inputs,
failed attempts and actual outcomes are retained under
`output_nyl12_local_repair_20261021/`; full specifications are in
`docs/relay_water_repair_20261021.md`. No accepted minimum is yet claimed.

Two neutral free-NH2 L-threonine methylamide conformers were prepared with
RDKit ETKDGv3 using seeds 20261020 and 20261021 and initial psi values of -40
and +140 degrees. The C5H12N2O2 model has CA S/CB R stereochemistry, charge
zero and singlet multiplicity. Input embedding involved no energy minimization.
Gaussian 16 optimizations and harmonic frequencies use omegaB97X-D/def2-TZVP,
UltraFine integration, tight SCF/optimization and no added dispersion term
[R4,R5]. The functional has peptide-parameterization precedent [R6], although
that study used a different basis and higher-level single-point energies.

A separate HF/6-31G(d) Merz–Kollman ESP step follows each optimization, with
`Pop=MK IOp(6/33=2,6/42=6)`. Subsequent RESP fitting is planned [R7]. Acceptance
requires completed optimization, 57 positive vibrational frequencies, normal
termination, correct atom inventory/order and stereochemistry, conservative
connectivity/contact checks, and successful ESP extraction. Each job used four
cores and 4 GB memory. Both jobs completed on 2026-10-20/21 (local time) with
Gaussian 16 revision ES64L-G16RevC.02 and passed independent rechecking. The
minimum frequencies were 45.9612 and 42.2567 cm^-1, respectively; extracted ESP
grids contained 5636 and 5767 points. Original logs, independent acceptance
records and optimized-coordinate exports are retained [S8; 2026-10-21 source
ledger update]. No fitted charges
or new force-field parameters have yet been accepted. Compatibility with the
intended ff19SB environment requires validation of charge constraints, caps,
bonded parameters and relevant conformational energies [R8].

## Reaction paths and saddle-point validation — proposed, not executed

DFTB3/3ob QM/MM will first be assessed for local refinement and exploratory
paths [R9]. After TI1 repair, reverse bond formation and proton transfer will
be used to obtain ES/path guesses. Driven scans or NEB will provide TS starting
geometries, followed by Hessian-based saddle refinement. The installed Sander
interface supports DL-FIND P-RFO (`imin=1, ntmin=5`, `optalg='PRFO'`), central-
difference Hessians and Bofill updates [R3,S9]. Newer Sella and pysisyphus
implementations are comparison candidates, not adopted dependencies [R12,R13].

Final validation will require convergence in all movable coordinates, a fresh
Hessian with exactly one significant negative eigenvalue in the declared
active space, the intended reaction mode, and downhill connectivity to ES and
TI1. Finite-difference displacement and SCF convergence sensitivity will be
recorded. Driving restraints will be removed for saddle refinement. The
movable region, fixed environment and any remaining positional restraints
will be reported explicitly. An active-space saddle is not claimed to be an
unconstrained full-enzyme saddle. The Hessian must cover mobile MM coordinates
as well as the QM region. Formal IRC calculations and displaced endpoint
minimizations will be labeled distinctly. Higher-level benchmarks will assess
DFTB reaction energetics; single-point corrections do not validate stationary
points at the higher level.

## Design and independent prediction — proposed

The initial campaign will use focused LigandMPNN sequence redesign on retained
WT backbones, followed by independent AF3 prediction and physical relaxation
[R2,R10]. WT mature alpha/beta lengths (232/96 residues per repeat) will be
retained. Fixed-backbone sequence design does not constrain the subsequent
prediction/relaxation backbone. Maturation architecture and catalytic residues
will be protected. RFD3 will be an optional later comparison for selected short
pocket/entrance loops, initially using exact WT loop lengths and retained core/
interface context [R11]. Exact loop ranges await accepted structural/access
analysis; whole-fold backbone redesign is not the initial plan.

The assembly will be treated as four alpha/beta repeats rather than eight
identical chains. LigandMPNN will tie equivalent alpha positions and equivalent
beta positions separately, using explicit tied groups and equal weights with
`homo_oligomer=0`. Automatic homooligomer mode ties all parsed chains and is
inappropriate for the mature eight-chain alpha/beta representation. A geometric
audit of the retained scaffold is
D2-like; the experimental assembly/operators will be checked before a production
RFD3 symmetry specification. Symmetric sequence design does not require
identical instantaneous loop conformations or ligand occupancy. The exact
RFD3 combination of alpha/beta repeats, retained context and catalytic atomic
motifs still requires an input preflight. Variable-length variants are deferred
until explicit indel-aware mapping can regenerate residue selections, symmetry
groups, maturation checks and Amber/QM indices; current WT exporters are not
certified for indels.

Candidates will be evaluated for compatibility with both catalytic and
open-compatible apo/ES conformations, including substrate/product entry and
exit. Cross-backbone comparison is a project-level screening step, not a claim
that LigandMPNN's multiple-input option implements joint multistate design.
Matched WT controls, geometry/steric checks and, after parameter validation,
independent sampling from distinct starting conformations will assess gate
behaviour. Neither confidence, minimization nor failure to observe opening in
a short trajectory establishes an opening barrier or rate. RFD3 exposure hints
can discourage inappropriate ligand burial but do not guarantee an access path
or flexible loop. Full design/access requirements and source identities are
recorded in docs/redesign_backbone_symmetry_and_access_20261022.md.

Independent noncovalent ES
predictions will omit enzyme–substrate bonds and imposed catalytic contacts.
Every ligand will be evaluated against all four equivalent active sites using
a fixed register and common criteria. All site outcomes and symmetry mappings
will be retained; model-level any-site recovery will count a model once.
For multiple ligands, assignments must avoid duplicate occupancy attribution.
Site observations within one model are nested, not independent replicates.

The existing contact-restraint AF3 fork is available for a separately labeled
guided arm. Historical production used fork commit
`ef5dc91d65b4845f63e95ea03316f1613c6dfb5e` and generated 90 models with four
ligands each [S4]. Guidance diagnostics do not establish independent binding
or catalytic success. This campaign uses older PA6 names C12/O4/N3/C11, which
require an explicit graph mapping to current C54/O49/N59/C62 selections.

Candidates will also be screened in apo, TI1, acyl-enzyme/water and
TI2-compatible states. Detailed acylation/deacylation paths are reserved for
a shortlist, after coarse reference-cycle mapping. Barrier comparisons require
balanced atom/proton inventories and corresponding reactant basins. Electronic
barriers, sampled free energies and experimental turnover will remain distinct
observables; no activity improvement is inferred from confidence alone.

Design ranking will use chemistry/maturation/assembly gates followed by a
diverse Pareto shortlist. Primary evidence will be independent ES geometry
recovery across seeds and compatibility with catalytic and accessible states.
Local AF3 confidence, pose variability, physical quality, effective strain and
same-state noncovalent Rosetta interface_delta will provide supporting evidence.
Interface_delta remains optional: its incremental value and accuracy will be
assessed against matched WT/mutant controls before inclusion. ES reproducibility
specifically means noncovalent near-attack recovery without imposed catalytic
contacts; retaining TI1's covalent bond is not equivalent evidence. TI1 requires
state-specific tetrahedral, oxyanion, protonation and relay-water checks, and
TS1 requires independent saddle-point and endpoint validation.
All attempted samples, failures and missing outputs will be recorded; equivalent
sites within one prediction will not count as independent replicates.
Effective strain measures structural deformation rather than mechanical
energy, and low WT-relative deformation will not be required for coordinated
redesigns. Conventional separated-interface scoring will not be applied to
covalent intermediates or interpreted as binding free energy.

Independent no-MSA prediction is a planned sequence–structure consistency arm,
with WT, single mutants and designs evaluated under matched settings. MSA- and
template-assisted predictions, if used, will be labeled separately. LigandMPNN
does not require an MSA for sequence generation [R10]; favorable model
confidence does not establish folding stability, substrate affinity or activity.
Sequence convergence will monitor diversity and iteration behavior rather
than define scientific fitness. A bounded pilot will establish the evaluation
budget and plateau rule, followed by an independent final assessment. Full
definitions and implementation limits are recorded in
`docs/design_ranking_and_iteration_20261022.md`.

## Reproducibility — implemented infrastructure

Scientific preparation and calculations run through retained scripts in tmux.
Records preserve commands, explicit environments, seeds, input/source/tool
SHA-256 hashes, source snapshots, software versions and logs. Process exit
status and scientific acceptance are separate. Content-addressed caching is
restricted to accepted artifacts with matching identities. New attempts use
new directories; rejected attempts and original coordinates remain available.
Repository dependency versions are managed with uv. Manuscript source hashes
and decision records are retained alongside the draft.


## Constrained TI1 stability and native transfer — implemented checks

For the capped, three-dihedral-constrained reference, local stability will be
assessed from the Gaussian Cartesian gradient and Hessian using the reduced
Lagrangian Hessian, including multiplier-weighted second derivatives of the
nonlinear dihedral constraints. Six rigid motions are excluded. Central
differences at two step sizes check constraint derivatives; unconstrained
frequency signs alone do not establish constrained stability. Numerical
thresholds and analytic tests are retained in the parameter protocol and source
snapshot. At this checkpoint no completed TI1 Hessian has been accepted.

The native-state charge model is separately mapped to full PA6 L4, reacting
Thr-NH3+, the oxyanion and explicit relay water. A proposed RESP boundary
protects installed NTHR peptide C/O, NME and TIP3P charges while jointly fitting
the reactive region. Clean tleap exports independently verified these source
values. This proposal has no accepted native geometry, ESP or fitted charges;
fit/held-out ESP errors and transfer validation remain required. User visual
approval of interim trajectories is recorded separately from computational
acceptance. Neither this preparation nor test-suite success certifies the TI1
force field.

## Relay-water feasibility and staged design scope — evidence and proposal, 2026-10-24

A deterministic front-side oxygen search was performed in the original full
alpha/beta protein/substrate context. The user-specified oxygen defined the
front hemisphere in a camera-independent frame perpendicular to the terminal
N-to-OG1 axis. At 0.025 A distance and 0.5-degree azimuth sampling, none of
225,000 positions at 2.6–3.2 A from both atoms passed the declared steric
screen; heavy-only and opposite-side diagnostics also failed. These finite
fixed-conformation screens concern spatial compatibility, not activation
energies or mechanistic exclusion. Thresholds, atom maps, source hashes,
sensitivity analyses and failed placements are retained in the source ledger.

The subsequent design strategy is proposed, not executed. We will prioritize
a small set of chemically justified local-state models and explicit alternative
protonation/relay hypotheses, using the native structure and mapped catalytic
residues as constraints. Exploratory motif screening may use an intermediate
or transition-state-like template labelled with its actual provenance. A claim
of TS1 optimization or a quantitative barrier comparison will still require
the saddle, reaction-mode and endpoint tests specified above.

Multistate structural evaluation will cover independent noncovalent ES, TI1,
acyl enzyme and deacylation-compatible states, with matched WT controls and
apo/product access checks [R27–R33]. An initial focused LigandMPNN arm will
preserve WT mature lengths and explicit alpha/beta symmetry ties; local RFD3
backbone sampling is an optional arm rather than a replacement of the native
fold. Repeated structural predictions will be interpreted as model ensembles,
not equilibrium populations. For large polyamides, local QM truncation will
be separated from full oligomer specificity testing: distal conformers,
cleavage registers, chain ends and open-pocket geometry will be varied.
Neither activated-ester activity in a published design nor a predicted
polymer-bound pose establishes nylon turnover, L2 selectivity or processivity.

### Recovery of the WT PA6 ES ensemble and installation checks

On 24 October 2026, a scripted audit recovered the existing WT PA6 L4
noncovalent simulations (Boltz model 41). All three production logs reached
50,000,000 steps at 2 fs and reported normal completion; the declared thermostat
seeds were distinct. Original LEaP input specifies ff19SB, GAFF (first version)
and TIP3P. Retained AM1BCC intermediates and exact ligand-library/frcmod hash
matches support the reported charge provenance; the exact historical command
and software version remain to be recovered. The original catalytic Thr has
OH/NH2 connectivity after manual terminal neutralization. Small original
residue/ligand charge residuals were recorded without correction.

One geometry-eligible ES seed per replica was selected in each predeclared
20–45, 45–70 and 70–100 ns interval. Selection minimizes time from the interval
midpoint among historical Tier B candidates, rather than optimizing geometry.
Full-solvent coordinates, hydrogens and boxes were retained. The legacy CSV
uses a zero-origin time estimate; physical production-relative saved-frame
times have a 0.05 ns offset, preserved explicitly in the new records. Original
trajectories, inputs, parameters and source code were fingerprinted with SHA256.
Native residue identities were checked directly against the supplied n12c CIF.
These selected seeds do not constitute an unbiased population estimate.

Separate installation checks reproduced the shipped Amber DFTB3 water
minimization and short periodic N-methylacetamide/water energy series at printed
precision. These are engine regression tests, not validation of nylon reaction
energetics. No new enzyme QM/MM pathway, TI1 force field or production
free-energy profile resulted from these recovery and installation checks.

For a complementary hydration diagnostic, CPPTRAJ watershell analyzed all 2,000
saved frames per WT replica with periodic imaging and explicit WAT oxygen masks.
First-shell counts at 3.5 Å around Thr N, OG1 and their union were combined by
inclusion–exclusion to identify simultaneous contacts. A separate 4 Å union
shell measured local water availability. This analysis found no qualifying
contacts in the 6,000 saved frames. Results are conditional on the existing
classical ES/protonation model and 50 ps output interval, and do not exclude
short water visits, slow buried hydration or hydration following chemical change.
Independent seed geometry validation retained full-precision NetCDF metrics;
historical phi was reproduced from its actual rounded PDB source, with explicit
coordinate-precision bounds and hashes, rather than changing an angle tolerance.

An exploratory ES relaxation retained the complete 103,251-atom periodic
enzyme–PA6 L4–solvent system and its inherited ff19SB/GAFF1/AM1-BCC/TIP3P
parameters. The DFTB3/3ob-3-1 region contained 98 original atoms and eight
hydrogen link atoms (charge −2, singlet), including the catalytic N-terminal
Thr233, adjacent Thr234, the Phe113 backbone donor, Asn186 side chain,
Asp272/Asp274 carboxylates and the central substrate amide. No water was
forced into the QM region. Sander analytic gradients were used with electronic
embedding, QM Ewald/PME, 10 Å cutoffs and an SCC threshold of 10−10.
All 103,153 MM atoms were restrained to the original coordinates at
10 kcal mol−1 Å−2, including MM protein, substrate, solvent and ion atoms;
the 98 QM atoms had no positional restraints. This was a position-restrained
ES relaxation, not an unrestrained ES minimum. The initial 25-cycle run did
not converge (RMS gradient 2.4646 versus a 0.05 kcal mol−1 Å−1 target).
A continuation with a 500-cycle budget retained the same Hamiltonian and
original restraint reference and reached the declared RMS threshold after
107 additional cycles (0.042903 kcal mol−1 Å−1). This is gradient convergence;
no Hessian classification was performed. Optimizer snapshots are coordinate evaluations,
not a molecular-dynamics trajectory. This pilot does not establish TI1
stability, reaction barriers or validated native TI1 force-field parameters.

The subsequent phase-I coordinate-driving pilot retained these MM positional
restraints in addition to two harmonic reaction-coordinate biases. Its
positional reference remained the original preparation restart rather than
the previous reaction window. A separate endpoint-reversion diagnostic removed
both types of restraint (`ntr=0`, `nmropt=0`), with all atomic coordinates active,
no SHAKE and a fixed periodic cell. This run must converge and undergo geometry
review before being promoted to a new ES reference. The user's intended ES
reference protocol is fully unrestrained atomic minimization; the exploratory
restrained model must not be described as satisfying that requirement.
