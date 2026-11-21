# Methods — working draft

## Preliminary NEB-seeded MBAR analysis — executed 2026-11-21

The 63-window 300 K DFTB3/3ob QM/MM campaign used neutral catalytic Thr
N-terminus, ASP272/ASP274, a 98-real-atom QM region of charge −2, and the
original 4,165-atom 12 Å movable selection with fixed exterior. Separate
native Amber umbrellas biased q_acyl = d(C,N_leaving) − d(C,OG1) and
q_H = d(OG1,HG1) − d(N_leaving,HG1), with k = 100 kcal mol−1 Å−2
for each coordinate. Exact native piecewise potentials, including linear
tails and the V = k(q−q0)² convention without an extra factor of one-half,
entered MBAR. No proton-branch censoring was applied.

Accepted complete production blocks were frozen for analysis. Of 23,800
frames, 5,003 statistically thinned samples entered the MBAR reconstruction.
The primary marginal used 25 equal-width q_acyl bins and the density bin
containing the initial E–S umbrella center as reference, with edges
[−1.88410, −1.73231] Å. Fifty within-window 5 ps block-bootstrap MBAR refits
gave pointwise 95% intervals relative to that fixed reference. This procedure
does not establish slow-coordinate convergence or include Hamiltonian error.
Positive-weight direct histogram accumulation was checked against an
independent two-dimensional marginal and a high-dynamic-range regression.
Proton nearest-owner assignments and component distances were retained to
diagnose branch degeneracy. Detailed source hashes and the invalidated initial
histogram supplement are recorded in the
[preliminary analysis record](../docs/neb_us_preliminary_findings_20261121.md).


## ASH272 trajectory comparison — executed exploratory analysis

Two 100 ns WT Nyl12–PA6 L4 replicas with Asp272 protonated in all four
active sites, Asp274 deprotonated and the catalytic Thr233 N terminus neutral
were analyzed at 10 ps/frame after the recorded 0.5 ns equilibration.
Protein Cα and ligand heavy-atom RMSDs were measured after fitting each
replica to its first production frame. We calculated d1 (Thr233 OG1 to
scissile carbonyl C), d2 (Phe113 N to carbonyl O), d3 (Asn186 ND2 to
carbonyl O), the OG1–C–O Bürgi–Dunitz angle, and **signed** projected
Herschlag φ using the existing atom mapping and Tier classifier. Asp272 and
Asp274 O/O separations and side-chain torsions were measured in all four
sites. A geometric water-bridge/H-bond screen sampled every 50 ps with
donor–acceptor distance ≤3.5 Å, H–acceptor distance ≤2.5 Å, and D–H–A
angle ≥120°. Water occupancy is not a proton-transfer rate or a free-energy
estimate. The executed script, exact trajectories, topology hashes, CSVs,
figures, and selected restart/PDB structures are under
`qmmm/pa6_flipped_model30/06_ash272_comparison/04_md_validation_20261113/`.
The matched SI panel preserves the earlier trace/histogram palette and layout,
changing the φ axis to signed values.

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

For noncovalent PA6 L4 sampling, the primary soluble-oligomer state is the
79-atom zwitterion with an explicit N-terminal NH3+ group, a C-terminal COO-
group and net charge zero. This state is appropriate to test for exposed
oligomer termini near neutral pH, but terminal protonation remains dependent on
pH, hydration, dielectric environment and ion pairing. It is not assumed to
describe a terminus buried in bulk or crystalline nylon. Any neutral NH2/COOH
sensitivity model will be parameterized independently and will not share a RESP
fit with the zwitterion.

Four bound-like PA6 L4 conformers were selected from a template-conditioned
model and one snapshot from each of three independent WT trajectories. Their
heavy-atom backbone and amide torsions were retained while all remaining
coordinates were optimized with omegaB97X-D/def2-SVP and SMD water. Def2-SVP
was selected for this geometry-only stage after a def2-TZVP attempt spent more
than three hours in an unnecessary initial `CalcFC` Hessian without beginning
an optimization step. No Hessian is required for RESP geometry sampling. The
subsequent Merz-Kollman ESP convention remains gas-phase HF/6-31G(d), and one
accepted conformer will be compared with def2-TZVP to assess geometry-basis
sensitivity before the force field is certified.

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

## Reaction paths and saddle-point validation

DFTB3/3ob QM/MM will first be assessed for local refinement and exploratory
paths [R9]. After TI1 repair, reverse bond formation and proton transfer will
be used to obtain ES/path guesses. Driven scans or NEB will provide TS starting
geometries, followed by saddle-point refinement. The installed Sander
interface supports DL-FIND P-RFO (`imin=1, ntmin=5`, `optalg='PRFO'`), central-
difference Hessians and Bofill updates [R3,S9]. Newer Sella and pysisyphus
implementations are comparison candidates, not adopted dependencies [R12,R13].

### DL-FIND dimer saddle-point searches

Transition-state refinement was initiated using the native Sander interface to
DL-FIND [R3] and the improved dimer method of Kästner and Sherwood [R43]. The
dimer estimates local curvature from gradients at nearby geometries, rotates
toward a low-curvature direction, and translates uphill along that direction
while relaxing perpendicular directions. It therefore searches for a
first-order saddle point without constructing a full Hessian during the search.
The improved formulation uses limited-memory BFGS optimization to reduce the
number of gradient evaluations [R43].

For the initial flipped-pose Nyl12–PA6 acylation search, image 8 of the final
complete NEB band (cycle 499) supplied the starting coordinates. The initial
dimer direction was the central difference between images 7 and 9, scaled to
a 0.01 Å norm over the declared QM weighting mask. Full-system coordinates
were reconstructed from native active-atom XYZ output using the original E–S
reference for the fixed exterior; the reconstructed image energy was checked
against the corresponding native NEB energy before launch. Calculations
retained the same DFTB3/3OB QM/MM Hamiltonian, topology, protonation state,
98 physical QM atoms, eight link atoms, QM charge of −2, and 4,165 movable
atoms used for the preceding band. The initial search therefore retains
deprotonated Asp272/Asp274 and is distinct from the prospective ASH272 MD
comparison.

Searches used Cartesian coordinates, native `dimer=.true.`, `dimer_mode=1`,
`optalg='LBFGS'`, `dimer_tolrot=2`, and `dimer_maxrot=20`. Declared DL-FIND
controls were `tol=1e-4`, `tole=1e-5`, `maxstep=0.1`, and `maxcycle=1000`.
No reaction-coordinate or positional restraint potentials were applied;
atoms outside the declared movable region remained fixed. These are search
settings, rather than evidence of convergence. Candidate stationary points
require subsequent vibrational characterization showing one imaginary mode
associated with the intended reaction and downhill connectivity checks to
the relevant minima before assignment as a validated transition state.

For the current forward E-S-to-TI1 exploration, a single generalized
coordinate was defined as xi = d(Thr233 OG1-HG1) - d(Thr233 N-HG1) -
d(Thr233 OG1-substrate carbonyl C). A harmonic bias was advanced from
approximately -4.317 to -0.800 Å in 24 windows, with a maximum target increment
of 0.15 Å. No separate target was imposed on either proton transfer or
nucleophilic attack, allowing the two chemical events to advance
asynchronously. Calculations used the accepted 12 Å DL-FIND movable shell and
no positional restraints. Each window retained the component distances and
energy for discontinuity and hysteresis checks. This driven path is treated
only as a source of pathway and umbrella-sampling seeds; it is not interpreted
as a free-energy profile, activation barrier or validated transition state.
An initial two-restraint synchronous trial was terminated after two accepted
windows because it imposed an unsupported timing relationship between the two
chemical events.

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

### Embedded-QM Hamiltonian comparison on the direct proton-transfer path

Five fixed geometries spanning the direct Thr233 OG1-HG1 to N-terminal-N
coordinate were selected from the periodic DFTB3/MM drive. The common QM region
contained 98 real atoms and eight Amber link hydrogens (charge -2; singlet).
Each geometry was evaluated with DFTB3/3OB, GFN2-xTB 6.7.1 and
omegaB97X-D/def2-SVP. For method matching, Amber generated adjusted and imaged
MM point charges within 20 A with QM Ewald and PME disabled; the retained
Gaussian geometry and charge records were converted deterministically to xTB's
legacy point-charge format using a hardness of 999.0. No coordinates were
optimized in this benchmark.

Because the finite 20-A charge list contained 8,088--8,107 atoms depending on
the window, raw finite-cutoff relative energies were not interpreted as a
reaction profile. Method assessment used per-window relative-energy double
differences against omegaB97X-D. DFTB3 and GFN2-xTB gave RMSE values of 1.49
and 3.91 kcal/mol, respectively; their transferred-endpoint errors were -0.57
and +7.93 kcal/mol. The DFTB3-minus-DFT differences were then subtracted from
the original periodic-Ewald DFTB3/MM relative profile as a sparse diagnostic
correction. The resulting non-restraint energies were 0.00, 0.42, 7.27, 12.57
and 15.85 kcal/mol and remained monotonically uphill. These values are neither
a free-energy profile nor proof of stationary-point absence; they do not include
DFT/MM relaxation or Hessian classification.

### Reaction-specific fast-Hamiltonian selection and corrected attack branch

The isolated proton-transfer comparison was not used to select a production
Hamiltonian because it omits nucleophilic attack, carbonyl polarization and
oxyanion formation. A chemically controlled attack branch was instead started
from the exact 24th asynchronous endpoint, which retained all three protons on
the terminal nitrogen and had an OG1--carbonyl-C distance of 2.104 A. Only the
OG1--C distance was advanced, from 2.104 to 1.50 A in increments no larger than
0.075 A. Flat-bottom restraints allowed each terminal N--H distance to move
freely between 0.90 and 1.20 A, with outer walls at 0.80 and 1.35 A, thereby
preventing the untracked H1-to-OG1 exchange observed in the discarded
continuation. The established 12 A DL-FIND movable shell was retained without
global positional restraints. These biased geometries were classified as
benchmark seeds rather than stationary points.

Sparse E--S, transferred-proton, early-attack, late-attack and TI1-like
geometries will be evaluated with DFTB3/3OB, PM6, PM6-DH+, GFN2-xTB and
omegaB97X-D/def2-SVP. Promotion requires preservation of proton identities and
connectivity, reference-consistent state ordering, relative-energy RMSE no
greater than 3 kcal/mol, endpoint errors no greater than 5 kcal/mol, and
agreement in projected reactive-coordinate gradient directions. The production
free-energy path will then use the explicit two-dimensional coordinates
`d(OG1,HG1)-d(N,HG1)` and `d(OG1,C)`, with H1/H2 identity gates, rather than the
chemically degenerate one-dimensional composite coordinate.

## Experimental context and simulation temperature (reported 2026-11-11)

Nyl12 variants were assayed as crude E. coli cell lysates using 5 mg PA6 or
PA66 polymer powder per well of a 96-well PCR plate, in 200 mM Tris buffer
reported as pH 8.0, at 75 °C for 6 h. These conditions were supplied by the
experimental investigator; the reaction volume and temperature of pH adjustment
remain unconfirmed. The model30 simulations used 300 K. The soluble PA6 L4
models and fixed protonation states therefore do not reproduce the complete
polymer-powder/lysate assay. Potential-energy pathways are not assay-temperature
free-energy barriers. Tris temperature dependence further prevents assigning
the incubation pH from the reported value without knowing its measurement
temperature (Brooke et al., 2015, https://pmc.ncbi.nlm.nih.gov/articles/PMC8956001/).
A matched 348.15 K sampling comparison is a proposed follow-up, not completed
work. Detailed provenance: docs/assay_conditions_and_temperature_20261111.md.
## Deposited NylC structure panel and Asp–Asp distances

We examined all nine PDB entries listed in Table 1 of Negoro et al.
(*FEBS J.* 2023, doi:10.1111/febs.16755): 5XYG, 3AXG, 5XYO, 5XYT, 5Y0L,
5Y0M, 7YU0, 7YU1 and 7YU2. Deposited PDB coordinate files were retained
without minimization or protonation changes. For each resolved protein site,
we measured all four Euclidean distances between Asp306 OD1/OD2 and Asp308
OD1/OD2 and used the minimum as an oxygen-label-invariant site descriptor.
Any deposited alternative locations were retained in the complete pair ledger;
the minimum does not assign proton positions or establish hydrogen bonding.
Rigid crystallographic copies preserve these within-site distances and were
not counted as independent observations. An independent ATOM-record inventory
checked target-oxygen completeness in the 15 deposited 3AXG protein chains;
no sixteenth site was assumed or reconstructed.

The scripted extraction ran in tmux with retained command arguments,
environment, source snapshots, input SHA256 hashes and process status.
Authoritative JSON records and a site-level CSV are archived under
`analysis/negoro_2023_protonation/table1_job/panel/`; the runnable scripts are
`audit_table1_structures.py` and `audit_crystal_and_mapping.py` under
`analysis/negoro_2023_protonation/`. Process completion and coordinate coverage
were assessed separately from mechanistic interpretation. No molecular
dynamics, energy calculation or pKa determination was performed for this panel.

## Independent activation/attack grid (started 2026-11-18)

### Proton-occupancy audit and coordinate integration check

For each completed grid point we verified the restart digest and audited all
54 real QM hydrogens against their two nearest QM heavy atoms, with a
full-system nearest-heavy check for boundary escape. The declared geometric
thresholds were 1.25 Å for covalent assignment and 1.35 Å for the second
neighbor shared/ambiguous flag. Assignments are geometric descriptors, not
electronic bond orders or pKa estimates. We also evaluated smooth coordination
`CN(X,H)=sum_H [1+(r_XH/1.25 Å)^6]^-1` for the catalytic N and OG1, relay
water, Asp pair and substrate heteroatoms. The label-independent activation
coordinate is `CN(Thr N,H)−CN(Thr OG1,H)`; both component counts and the
attack distance remain explicit observables. A reusable implementation
provides its analytic Cartesian derivative. Synthetic-coordinate derivative
and hydrogen-permutation tests passed; 12 focused tests covered all 41
executable statements of the coordinate module. No molecular QM forces were
numerically differentiated for these unit tests.

Two one-cycle native DL-FIND evaluations on identical E–S coordinates tested
the proposed Sander NFE/PMD integration, with and without a deliberately
nonzero N_OF_BONDS restraint. The source, QM region, active shell and periodic
Hamiltonian were identical. Native PMD's half-k convention predicts a
0.20 kcal/mol contribution for the prescribed 0.20 coordinate offset and
spring constant 10. The observed energy difference was zero, and PMD/CV
initialization/output was absent. We rejected this integration backend rather
than interpreting process exit zero as evidence of an active bias. The
original 0.40 kcal/mol expectation was retained with an explicit convention
correction; no installation was modified.

A separately authorized unbiased release started from q_004_attack_001,
retaining the accepted source QM/MM input with `nmropt=0`, `ntr=0` and no
reaction restraints. Only the two cycle limits were increased to 5000.
Force/SCC completion and final chemical occupancy were recorded separately.
The initial structure, saved optimization frames (every 25 steps) and final
restart were exported as a 13-frame stripped trajectory retaining WAT1576.
All inputs, source copies, hashes, process records and scripts are archived
under the grid's `proton_occupancy_audit/`, `occupancy_coordinate_preflight/`
and `diagnostic_activation_release/` directories.

### Analytic occupancy-bias callback and bounded pilot

After rejecting the inactive native PMD integration, we added the tested
coordination bias and its analytic Cartesian derivative to native Sander
energy/gradient evaluations through an explicit optimizer callback. Each native
call evaluated one requested geometry with the unchanged periodic DFTB3/MM
Hamiltonian; the retained DL-FIND force-output conversion used the previously
validated source-hash-matched calibration. The native energy was converted
from hartree to kcal/mol before adding the bias; gradients were expressed in
kcal/mol/Å. Native input/output coordinates, single energy/force-frame counts,
SCC status and artifact hashes were checked. An actual two-iteration L-BFGS-B
proof with four native evaluations reproduced a deliberately prescribed
0.20 kcal/mol bias and verified nonzero analytic force addition. The proof's
iteration limit was intentional and was not a stationary-point criterion.

A four-point serial diagnostic started at q_004_attack_001, with
q_occ=CN(N,H)−CN(OG1,H) targeted at 2.10 and independent OG1–C targets
3.1251784, 2.60, 2.10 and 1.60 Å. The bias was
0.5*100*(q_occ−q_target)^2 + 0.5*200*(r−r_target)^2 kcal/mol, with r in Å.
All 54 real QM H atoms contributed; both component counts and real-H
ownership were retained to detect unintended chemistry. The original 102-real-
atom QM region, charge −1, eight links and 4,662-atom movable shell were
unchanged. No additional positional restraints or bounds were introduced.
SciPy L-BFGS-B used scaled variables y=x/(0.1 Å), with corresponding analytic
gradient scaling; this variable scaling does not impose a step cap. Each point
allowed 1,500 iterations and 3,000 function evaluations. A maximum absolute
active biased-gradient component ≤0.12 kcal/mol/Å was required for serial
continuation; energy-based termination alone was insufficient. Accepted
optimizer iterations were distinguished from line-search trial evaluations.
These four diagnostic candidates do not constitute a converged minimum-energy
reaction path or a free-energy calculation.

The initial A minimization's L-BFGS-B energy line search terminated before
the force criterion was met. A retained-trial energy/gradient audit localized
the discontinuity to an ordinary MM pair crossing the unchanged 10 Å
nonbonded cutoff. We therefore resumed from its final accepted geometry using
the installed DL-FIND gradient-based L-BFGS step routine, without an energy
line search. The small ctypes interface records the inspected source and
shared-library hashes and enforces exclusive ownership of native global state.
Inputs to that routine were in bohr and hartree/bohr; proposed displacements
were capped at 0.05 Å per atom. Nonpositive secant curvature caused an explicit,
recorded history reset. The force criterion, physical energy/gradient
evaluation, bias, QM partition, active shell and all electronic settings were
retained. Original failed records were preserved and the continuation appended
native evaluations and accepted snapshots in the same logical point directory.
Runtime and scientific convergence were recorded separately.

After all four points converged, we audited their exact final restart/input
coordinate agreement, retained artifact hashes, native energies, recomputed
bias energies, all real-QM hydrogen owners and original heavy-atom bonds.
Carbonyl pyramidalization was measured as the signed height above the plane
of its three original bonded neighbors, together with their three-angle sum.
We measured both oxyanion-hole donor distances, attack geometry and substrate
C–O/C–N distances; restrained force convergence was not treated as TI1 identity.

A single diagnostic release started at D with both harmonic spring constants
set exactly to zero. Every recorded optimizer state required zero added bias
and equal total/physical energy. Native settings, partition, active shell,
gradient-based optimizer and force tolerance were retained. Saved accepted
snapshots and all native trial inputs/outputs were retained separately. Final
proton ownership, chemistry and force completion were to be evaluated before
any intermediate claim; a Hessian was not included in this release protocol.
In parallel, one exact-D finite-embedding force comparison was initiated with
ωB97X-D/def2-SVP and GFN2-xTB, using the same 110 QM/link atoms and 104,824
adjusted MM charges as the preceding reference protocol. Additional dual
directions shortened substrate O–Phe113 N or O–Asn186 ND2 while holding
proton occupancy and attack distance fixed to first order. The respective
three-coordinate Jacobians were verified full rank and no cut parent moved.

### Saved reference-gradient projections

Without new electronic-structure evaluations, we parsed saved analytic gradients
for E–S, q_004_attack_001 and q_003_attack_005. Gaussian and GFN2-xTB models
contained the same 110 QM/link atoms and 104,824 adjusted MM point charges.
Atom order, coordinates and embedding membership were verified. Local dual
directions obtained from the two-CV gradient Gram matrix changed q_occ or
OG1–C independently while holding the other constant to first order. None
moved a cut parent, so link-position chain-rule terms were zero for these
directions. Dot products yielded local directional derivatives, not relaxed
barriers. Embedded-QM-only Gaussian/xTB derivatives and native total DFTB3/MM
derivatives were reported separately because their energy components differ.
Reference outputs and analysis are retained in
`08_activation_attack_grid/reference_checks/gradient_analysis/`.

### Independent-coordinate grid protocol

To vary proton transfer independently of nucleophilic attack, a 5 × 5 relaxed
PES grid was initiated from the accepted water-in-QM E–S endpoint. The two
restrained coordinates were `q_N = d(OG1,HG1) − d(Thr N,HG1)` and
`r_attack = d(OG1,substrate carbonyl C)`. Five q_N targets span the measured
E–S value (approximately −1.778 Å) to +1.050 Å; five attack targets span
3.099537 to 1.500 Å. Amber harmonic restraint parameters rk2=rk3=100 were
applied to each coordinate. Five one-core row workers run concurrently;
each row first approaches its q_N target at the E–S attack distance and then
sweeps attack at fixed q_N. Short target increments (at most 0.15 Å for q_N
and 0.10 Å for attack) give 135 optimization stages for 25 reported points.
The source topology, 102 real QM atoms, eight links, charge −1, periodic
electrostatics and 4,662-atom movable region were retained. These row sweeps
are deterministic continuations from one seed, not independent MD replicas
or a converged free-energy surface. Measured coordinates, target lag, proton
inventory, process/force acceptance and physical/bias-inclusive energies are
retained separately. No basin or barrier is established while this grid runs.

Sparse ωB97X-D/def2-SVP and GFN2-xTB comparisons use an explicitly separate
finite nonperiodic embedding model with full MM point-charge membership,
audited link treatment and retained analytic-gradient output. The initial
E–S Gaussian force calculation was started with eight CPU cores and 32 GB.
The exported field contains 104,824 adjusted MM charges after excluding the
eight MM link parents. An initial Sander/EXTERN invocation returned exit zero
despite Gaussian exceeding its field-query centre limit; its energy and
zero-atom formatted checkpoint were rejected. The corrected reference driver
runs Gaussian Force directly, retaining the same geometry and charges but
omitting the additional MM field-query coordinates. Acceptance requires normal
Gaussian termination, 110 QM/link atoms and 330 finite Cartesian gradient
components. This gives local embedded-QM gradients; saved-gradient projections
were subsequently analyzed as described above. A complete periodic
DFTB3-versus-DFT total-gradient comparison has not been performed.
Cross-method comparisons require common adjusted charge membership and
consistent point-charge self-energy conventions. Such comparisons diagnose
local electronic-method sensitivity; finite-cluster energies must not be
substituted directly into the periodic-grid PES profile. Exact inputs,
source snapshots and launcher/process records are under
`qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid/`.

## Water-mediated Thr activation diagnostic (executed 2026-11-17)

The accepted ASP272/ASH274 E–S restart with neutral Thr233 N was optimized
using the unchanged 102-real-atom, charge −1 DFTB3/MM partition including
WAT1576 and a 4,662-atom, complete-residue 12 Å movable region. A sequential
25-window scan restrained
`Q = d(OG1,HG1) − d(WO,HG1) + d(WO,H2) − d(Asp272 OD1,H2)`
from −1.692959 to +1.500 Å and independently restrained OG1–carbonyl C
from 3.099537 to 1.650 Å. Both target schedules advanced across the same
windows, although the individual proton transfers could progress unequally.
These are seed-generating restraints; their schedule is not a minimum-energy
path determination. The authoritative force constants and engine settings
are retained in each window's `restraints.RST` and `drive.mdin`.

All calculations used scripted tmux execution with input and output hashes,
force-convergence and chemical-inventory checks, and cumulative VMD/energy/
distance exports. A continuity check initially rejected window 5 because it
included a distant spectator separation. Protocol revision 5 restricted this
check to seven declared reactive distances, preserved the original decision,
validated the five completed windows, and resumed from window 6 without
recalculating them. The chemistry, Hamiltonian, shell, and target schedules
were unchanged. Physical relative energies were computed as
`(EPtot − RESTRAINT)_window − EPtot_ES`; bias-inclusive energies were retained
separately. An unbiased endpoint release was conditional on the endpoint
chemistry and tetrahedral-geometry gates; the final endpoint failed that
gate, so no release was performed. Protocol, scripts, and acceptance records
are under
`qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/`.

## Existing-MD Nalpha-to-leaving-N water screen (executed 2026-11-18)

We screened all frames eligible under the cached Tier A/B and signed-φ
criteria in six existing 100 ns production replicas (10,000 saved frames
each, 10 ps spacing). ASP272/ASH274 was primary; the two ASH272/ASP274
and two ASP272/ASP274 replicas were separately authorized controls. All
used their existing neutral catalytic Thr233 NH2 topologies and source MD
trajectories. Controls were not promoted to the primary microstate.

The existing geometric thresholds were retained: Tier A required d1 ≤3.2 Å,
d2/d3 ≤3.8 Å, Bürgi–Dunitz angle 90–120°, and |φ| 70–110°; Tier B allowed
d1 ≤3.5 Å, d2/d3 ≤4.0 Å, angle 85–125°, and |φ| 60–125°. Only positive
signed φ frames were retained. Here d1 is Thr233 OG1–scissile C, d2 is
Phe113 N–carbonyl O, and d3 is Asn186 ND2–carbonyl O. Native Thr233 was
held out from a Kabsch fit to backbone N/CA/C/O of Phe113, Asn186, Asp272
and Asp274 against `cif/n12c.cif` chain A. The unchanged native limits were
OG1 displacement ≤0.75 Å, Nalpha displacement ≤1.25 Å, |Δψ| ≤35°, and
|Δχ1| ≤20°. The intact-PA6 check required every topology bond within the
four PA6 residues (Amber 1313–1316) to exceed 0.65 Å and remain ≤1.30 Å
when either atom was H, or ≤1.85 Å otherwise. All Tier A/B positive-φ
frames passed this bond-retention check.

An observed one-water bridge required the same topology-owned water O
within 3.5 Å of both Thr233 Nalpha and the scissile leaving amide N.
Productive leaving-leg geometry additionally required an O–H···Nleave
contact with H···N ≤2.5 Å and O–H···N angle ≥120°. An observed Nalpha–H···O
contact used the same donor–acceptor/H–acceptor/angle criteria and defined
a separate strict two-leg class when both legs passed. Absence of that
incoming contact in neutral NH2 E–S did not alone reject a leaving-leg
candidate: the future proton transferred from OG1 is not present on Nalpha
in the screened state. Distances used exact triclinic minimum images,
including the truncated-octahedral cells; fractional rounding alone was not
used as the final distance criterion. Atom names, covalent hydrogen owners,
neutral NH2 inventory and declared Asp state were read from each actual MD
topology. No water coordinates, hydrogens or protonation states were edited.

| Fixed microstate | Replica | Frames | Tier A/B, φ > 0 | All structural gates | One-water bridges |
| --- | ---: | ---: | ---: | ---: | ---: |
| ASP272/ASH274 | 1 | 10,000 | 5,865 | 8 | 0 |
| ASP272/ASH274 | 2 | 10,000 | 30 | 0 | 0 |
| ASP272/ASP274 | 1 | 10,000 | 819 | 0 | 0 |
| ASP272/ASP274 | 2 | 10,000 | 1,296 | 31 | 0 |
| ASH272/ASP274 | 1 | 10,000 | 7,253 | 12 | 0 |
| ASH272/ASP274 | 2 | 10,000 | 3,695 | 0 | 0 |
| Total | — | 60,000 | 18,958 | 51 | 0 |

“All structural gates” means Tier A/B, positive φ, intact PA6 and native
Thr233 together. “One-water bridges” denotes the heavy-atom prerequisite;
its zero count also makes both hydrogen-bond classes zero. Frames are
correlated; these are geometric counts, not independent observations,
protonation populations or transfer rates. This mature-substrate hydrolysis
screen does not address precursor autocleavage.

Per-frame gate tables, all input/trajectory SHA-256 hashes, source snapshots
and separate process/scientific records are retained in
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/00_screen/`.
`primary_attempt1/` and `controls_attempt1/` retain the two screen stages;
`diagnostics_attempt1/` retains eight primary nearpasses and overlapping
native-gate failure distributions. `bounded_screen_acceptance.json` and
`bounded_screen_counts.csv` give the combined counts. The three tmux jobs
completed with exit status zero, while candidate acceptance was negative:
no frame was extracted, no new MD or QM/MM was launched by this screen,
and the primary `selected.json` has null primary/reserve entries.


## Balanced acyl-enzyme endpoint and asynchronous reverse-path test — initiated 2026-11-18

An independent acyl-enzyme/amine coordinate hypothesis was constructed from the
accepted water-in-QM E–S structure. Atom identities, total proton inventory and
QM charge (−1) were retained. Both fragments of the scissile amide were kept;
Thr OG1 formed the ester, the amide C–N bond was separated, and the departing
alkylamine acquired a second proton. One original water proton was assigned to
the departing amine, while original Thr HG1 replaced it on the same water. This
label-balanced construction does not establish the actual proton-transfer
route. The beta-chain N termini remained initially neutral and the initial
acid microstate was ASP272/ASH274.

A bounded coordinate-only search tested 285 rigid leaving-fragment placements
and a 24-value acyl-arm torsion sweep. Internal fragment geometry was preserved
except for explicitly constructed ester/amine coordinates. The selected guess
had no severe nonbonded contacts using thresholds 1.1 Å for heavy–heavy, 1.0 Å
for heavy–H and 0.8 Å for H–H pairs, excluding product 1–2/1–3 neighbors. These
thresholds are gross-overlap preparation checks, not measures of relaxed energy
or hydrogen-bond quality. Maximum heavy/H shifts were 2.172/2.366 Å; protein
and water oxygen coordinates were unchanged during construction. The acute
83.97° starting ester angle was left for energy optimization.

The initiated QM/MM workflow retains 102 real QM atoms, 54 real QM H and eight
link H, with the same periodic DFTB3/3ob Hamiltonian, zero electronic temperature
and complete-residue 12 Å movable region (4,662 atoms). Temporary harmonic
terms target OG1–C=1.4 Å, C–Nleave=3.0 Å and the seed's smooth leaving-N
proton coordination; springs are 200 kcal mol−1 Å−2, 50 kcal mol−1 Å−2 and
100 kcal mol−1 per squared coordination unit, respectively. All terms are
removed for a separate unbiased optimization. No classical covalent-product
parameters are substituted. The maximum active Cartesian gradient criterion
remains 0.12 kcal mol−1 Å−1; endpoint chemical inventory and curvature
verification are assessed separately.

Conditional on force convergence and balanced acyl-enzyme retention, 36 reverse
windows target the scalar ξ = r(OG1,C) − r(C,Nleave) + (1 Å)
[CN(OG1,H) − CN(Nleave,H)] toward the accepted E–S value, using a
100 kcal mol−1 Å−2 harmonic spring. Each CN sums all 54 real QM hydrogens with
the previously defined smooth switching function. Nalpha occupancy and the
individual component coordinates are unrestrained, permitting asynchronous
proton motion, attack and scission. The scalar may still admit chemically
different structures; a final zero-bias release and explicit connectivity audit
are required. This is a potential-energy path search, not a PMF or verified TS.
Native inputs, construction choices, physical/bias energies, all proton owners
and process/scientific outcomes are retained in the campaign manifests.


### Parallel NEB numerical-path diagnostics (initiated 2026-11-19)

Four violet trials were initialized from the retained ASP272/ASH274 NEB
cycle-268 band: uniform 24-image L-BFGS, uniform 48-image L-BFGS, adaptive
24-image L-BFGS and uniform 24-image CG. Each used 13 controls; adaptive
placement doubled the QM arc-length sampling weight near source images 9–10.
Native initialization remained piecewise Cartesian linear interpolation.
The common component step cap was 0.02 bohr. The periodic DFTB3/3ob Hamiltonian,
electronic temperature 0 K, QM charge −1, 102 real QM atoms/54 real H/eight
link atoms, neutral Thr Nalpha, frozen endpoints/exterior and original
12 Å/4,662-atom movable selection were unchanged; no additional restraints
were introduced. L-BFGS remained primary. CG was limited to 1,750 evaluations
(about 79 full 24-image cycles) as a curvature-history diagnostic.

Installed-source inspection identified a hardcoded >90° kink repair applied
after step capping. Comparisons therefore record native Max/RMS gradients and
target ratios, reset counts/timing, coherent image/chemical gaps, elapsed time
and actual energy-evaluation counts at 50/100/250 cycles where feasible, rather
than rank transient profile maxima. Frozen-endpoint 24/48-image bands require
22/46 interior evaluations per full cycle. One shared reacting H20318 in the
48-image seed was accepted explicitly for that unconverged seed only, with its
strict retention failure and SHA-pinned review retained; spectator gates were
unchanged. Separate startup reviews accepted native mapping at 1e-6 Å while
preserving original 1e-7 Å diagnostic failures. Neither review constitutes
force/chemical convergence. Configurations, reviews, coherent movies and
source-hashed comparison histories are retained under
`09_water_relay_acylation/04_neb/spacing_tests/` in the ASH274 campaign.

## Exploratory pKa screens: mature/precursor NylC and mature Nyl12

PROPKA 3.5.1 was applied to protein-only coordinates as a static sensitivity
screen. For mature Nyl12, three Model30 Tier A MD frames represented geometries
sampled under ASP272/ASP274, ASH272/ASP274 and ASP272/ASH274 assignments.
All eight mature chains and explicit chain boundaries were retained. Protein
hydrogens and the nonstandard PA6 ligand were omitted, and ASH heavy atoms
were normalized to ASP for the common titratable-group model. Consequently,
these predictions compare conformations; they do not estimate the free energies
or populations of the original fixed-charge MD microstates. All four beta
sites were reported, including the catalytic Thr233 free amino terminus.

For NylC, we compared mature 5XYG with uncleaved 7YU0 and 7YU1 using the
deposited sequences, without back-mutation to WT. The precursor structures
contain H130Y/N266A/T267A and D122G/H130Y/T267C substitutions, respectively.
Thr267 is internal in those precursor chains; its free amino terminus is
represented only after cleavage in mature 5XYG. A malformed TER record in the
initial mature input prevented terminal recognition. Correcting only the
terminal record restored the Thr267 N-terminal group without changing atom
coordinates. Only the corrected attempt3 mature results are retained for
interpretation. Precursor results are from attempt2; its earlier mature
values are superseded. Default and coupled-alternative crystal predictions
were retained separately. The original Nyl12 screen did not inspect coupled
alternatives. The pKa subset is distinct from our nine-entry structural
Asp–Asp distance audit, which itself performed no pKa calculation.

An exploratory MCCE4-Alpha/DelPhi calculation completed for mature NylC, but
failed structural/scientific acceptance because internal coordinate gaps
created 13 spurious terminal groups in the alpha chains. Its extreme titration
bounds are diagnostic failures, not accepted pKas. A subsequent complete
Nyl12 protein-only step-1 preflight retained eight chains but recognized zero
NTR and zero CTR groups and emitted 2,840 error diagnostics. No full Nyl12
MCCE titration was launched; no accepted precursor MCCE result was obtained.
Original inputs, commands, hashes and explicit rejection receipts are
preserved. Detailed results and source paths are indexed in
`docs/propka_asp_pair_pilot_20261116.md` and the student handoff's
`evidence/pka/source_mapping.json`. No calculation was rerun for this
documentation consolidation.


## NEB-seeded direct-acylation umbrella sampling (2026-11-19)

The force-converged, regular 24-image E–S-to-acylenzyme NEB for the
neutral-Thr233-N/ASP272/ASP274 microstate was used to seed a separate
DFTB3/3ob QM/MM umbrella campaign (98 real QM atoms, charge −2; electronic
0 K). The NEB maximum, 35.80895 kcal mol⁻¹ relative to its reactant endpoint,
was not assumed to be a verified transition state or a free-energy barrier.
The Hamiltonian and exact original 4,165-atom movable region were retained;
100,757 exterior atoms were frozen. The resulting PMF is conditional on this
boundary and fixed proton inventory. Legacy Amber GROUP cards encoded the
movable selection; a native probe verified its identity, zero frozen-atom
displacement and thermostat degrees of freedom (12,495).

Separate native Amber potentials biased q_acyl=d(C,N_leaving)−d(C,OG1) and
q_H=d(OG1,HG1)−d(N_leaving,HG1), allowing asynchronous fluctuations. Each
Amber RK value was 100 kcal mol⁻¹ Å⁻², with central potential k(q−q0)² and
native linear tails outside ±10 Å. No NEB spring or Cartesian image potential
was retained. All 24 image centers were kept; additional centers gave 63
windows with adjacent two-coordinate Euclidean spacing ≤0.10 Å. Each used
the nearest validated original full-system restart, without coordinate
interpolation. This path-centered sampling does not establish complete
coverage of the two-dimensional chemical space.

At 300 K, each window was allocated 1 ps heating, 5 ps equilibration and
six 5 ps production blocks, with a 0.5 fs step, Langevin friction 2 ps⁻¹,
independent recorded stage seeds, fixed volume and no bond constraints.
Production allocation was at least 30 ps/window; completion alone was not
sampling convergence. Full-precision coordinate-derived CVs were used to
cross-evaluate both native biases in every window. CPU-only PyMBAR analysis
retained valid samples from competing proton branches, assessed correlation,
effective samples and overlap, and withheld global surfaces for disconnected
overlap graphs. In addition to the HG1-labeled coordinate, a smooth
permutation-invariant N/O–H coordination difference monitored proton transfer
over the nine exchangeable real-QM hydrogens; Thr-N coordination and labeled
bond distances were retained. Link H atoms were excluded from this pool.
Unsampled bins were not interpolated. Interpretation and sampling extensions
require block/replica stability and adequate orthogonal chemical coverage.
Protocol, code/input hashes and live acceptance records are retained in the
parent NEB directory's `umbrella_sampling/README.md` and manifests.


Distributed umbrella execution retained byte-identical scientific inputs
across violet and azalea, using separately hashed native executable builds.
A 100-step compatibility probe reproduced the reference initial potential
energy at printed precision (acceptance tolerance 0.001 kcal mol⁻¹). A static
manifest assigned 16 windows to violet and 47 to azalea with eight and 24
single-threaded workers, respectively. Ongoing violet windows were adopted
without restarting MD. Completed remote observables and cross-bias matrices
were checksum-validated and published on the analysis host with acceptance
metadata last, avoiding partially transferred blocks in MBAR. Full native
trajectories remained on their execution hosts; compact structural views
and analysis products were mirrored separately. Neither the distribution
nor the control-process library isolation changed the QM/MM Hamiltonian,
thermostat, seed schedule, frozen exterior or chemical bias definitions.
