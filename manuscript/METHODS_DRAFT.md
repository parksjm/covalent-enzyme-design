# Methods — working draft

Status: 2026-10-21. Sections explicitly distinguish completed preparation from
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
beta-chain NH2 group and Thr OG1H. A neutral relay water is included in the
working mechanism. These state definitions are audited hypotheses; predictor
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

## Parameter-model quantum calculations — executed; RESP pending

### Capped anionic TI1 bonded reference — prepared; optimization running

A local bonded-parameter surrogate was constructed as
ACE-THR(-O-C*(O−)(CH3)(NHCH3))-NME, retaining the tetrahedral carbon and intact
leaving-group C*–N bond. Its explicit-hydrogen graph contains 37 atoms
(C10H20N3O4−), charge -1 and singlet multiplicity; Thr CA/CB configurations
are S/R. Methyl groups truncate the reacting substrate amide. The surrogate
has neutral peptide caps and does not represent the native catalytic NH3+
group, full substrate termini or relay water. Its charges cannot be transferred
directly to the native TI1 model.

The user-selected protein reference `cif/n12c.cif` was parsed with an
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

### Exploratory local relay-water refinement — in progress

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

Validated catalytic atomic motifs will guide focused, symmetry-tied sequence
redesign with protected maturation/catalytic residues. Local RFD3 backbone
diversification is a later option [R10,R11]. Independent noncovalent ES
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
