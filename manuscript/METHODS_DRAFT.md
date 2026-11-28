# Methods — working draft

## Prespecified TS-ensemble redesign and validation workflow

For the initial matched sequence-validation cohort, the two highest-ranked
original-backbone candidates were selected using local-pocket LigandMPNN
scores across both original TS-like parents. Candidates 06 and 03 contained
92 and 99 substitutions, respectively, and differed at 41 of 328 precursor
positions while retaining all 101 protected identities and four-copy sequence
ties. Scores served only as sequence/structure diagnostics. WT and both
designs entered a fresh matched two-seed, template-free Boltz comparison on
Violet, with four unguided PA6 L4 molecules. Each design used its own exact
328-residue precursor MSA search and query-aligned mature-chain pairing;
the unchanged validated WT MSA was reused. The queued workflow applies identical unrestrained Amber relaxation after
canonicalization and severe-clash screening. Positive-phi Tier A/B selection
considers all 48 site–ligand–amide combinations before one-to-one assignment,
with all four sites reported before and after relaxation. No comparative outcome
is yet claimed for this cohort. These are structural tests, not estimates of
catalytic efficiency.

The bounded redesign protocol selects chemically coherent, decorrelated
representatives near supported relay free-energy bottlenecks, retaining source
window, block, frame, proton inventory and reaction-atom mapping. Templates
selected before free-energy and reaction-coordinate convergence are designated
TS-like structures; they are not claimed to be verified first-order PES saddle
points or committor-defined transition states. Original backbones provide a
baseline. The current design workflow excludes RFD3. Its candidate pool comprises the
12 fixed-backbone LigandMPNN designs from two chemically coherent MD-derived
TS-like parents, with the 101-position protection mask and identical sequences
across all four copies. Thermal conformational asymmetry is allowed; exact D2
coordinates are not imposed on MD snapshots. Parent acceptance independently
requires reaction mapping, proton inventory, ligand graph and catalytic geometry;
multiple chain conformations alone do not establish TS compatibility. Historical
RFD3 remodels and their 12 sequences remain diagnostic artifacts outside this
current pool; no replacement generation or RFD3 qualification is required.

Fixed-backbone LigandMPNN design followed by independent template-free,
paired-MSA Boltz prediction and matched Amber relaxation can assess
sequence-driven repacking and assembly; it does not promise improved activity.
Only clear evidence of a backbone limitation would motivate a future controlled
RFD3 comparison evaluating coupled perturbations of A/E, B/F,
C/G and D/H, with independent physical-movement, clash, assembly, fixed-context
and chemical checks. It is not implemented or qualified for production, and
exact coordinate symmetry is not a universal acceptance requirement.

LigandMPNN sampling ties equivalent alpha and beta positions across the four
copies. The current design mask fixes 101 of 328 precursor positions: 24 core
policy positions and 77 additional recurring inter-repeat contact positions.
The structural contact rule requires a heavy-atom distance of at most 4.0 Å
in at least three of four equivalent copies in each of three original TS-like
parents, excluding the intrinsic A/E, B/F, C/G and D/H alpha/beta pairs.
This is an operational assembly constraint, not evidence that every fixed
position is experimentally essential. The full mask, WT identities, overlapping
reasons, source coordinates and reproducible visualization are provided in
the Supporting Information. Sequence protection is distinct from backbone
coordinate restraint.

The historical pilot generated 24 raw samples, six on each of two original and
two remodeled parents accepted under the historical graph/coordinate criteria.
Only the 12 original-parent sequences enter the current design pool; these
receipts did not qualify an RFD3 production workflow. Complete sequences are deduplicated while
retaining every lineage, and WT/design cross-scores use identical model weights,
masks and decoding orders across TS-like and matched E–S structures. Likelihood
rank is used for computational triage, not as an activation-energy estimate.
At most two initial design representatives and WT enter a matched two-seed,
template-free, precursor-paired-MSA Boltz comparison, with four unguided PA6 L4
molecules per assembly. AF3 is an optional separately calibrated comparison.

Raw predictions are canonicalized and checked for chain identity, substrate
connectivity and severe interchain overlaps before matched unrestrained Amber
relaxation. The existing OpenMM adapter is an alternative backend; it is not
an additional required stage and cannot replace only the design arm without
matched control qualification. Canonical models are measured before and after
relaxation over all 48 site–ligand–amide combinations. Positive signed Herschlag
phi and the declared Tier A/B criteria are applied before pose assignment;
native-Thr similarity remains diagnostic. A small, diverse group of finalists
can proceed to matched independent MD replicas, with MM/QM/MM optimization
and reaction-path tests limited further to persistent, chemically compatible
candidates. US/MBAR is reserved for the strongest cases with a matched WT
reference. All stages retain input/source hashes, maps, environments, process
receipts and separate scientific acceptance. Detailed stage gates and reused
entry points are recorded in the [protocol](../docs/NYL12_TS_ENSEMBLE_REDESIGN_PROTOCOL.md).

## Supplementary interface-confidence analysis

Directional ipSAE scoring with ipSAE2 [R64] was technically qualified on matched
WT and design-06 Boltz predictions. Version 0.1.11 at commit
`fe560f4c2461eeb91aebded61ac8c1e83b264f60` was installed in a separate,
uv-locked Python environment. Each unmodified native CIF was bound to its
matching PAE, confidence, token-pLDDT and processed-structure artifacts by
hashes. Native Boltz token records, atom order and pLDDT-to-CIF correspondence
verified PAE-axis identity before canonicalization; dimensions alone were
insufficient for identical protein copies and ligands. Protein-only structural
correspondence established mature alpha/beta pairing independently for each
prediction. No PAE transposition was applied.

Original chain pairs, virtual mature protomers and the whole-protein group
were scored separately, retaining both directions and all four ligands. For
protein groups compared with a ligand, the explicitly protein-aligned score
was retained because the grouped scalar maximum can select the opposite
direction. Scores annotate native predictions through explicit identity maps;
they are not newly estimated confidence values for Amber-relaxed coordinates.
No confidence threshold substitutes for positive-phi Tier A/B or physical
checks. PA6-specific calibration, binding affinity and catalytic-activity
prediction are not claimed. Technical checks and first scores are reported
in [Supporting Information](SUPPORTING_INFORMATION_DRAFT.md#supplementary-ipsae2-qualification--28-november-2026).

## Precursor-derived paired MSAs and template-free WT controls

A local ColabFold/MMseqs search used the exact 328-residue WT precursor
against the installed UniRef30 and environmental databases. Homolog sequences
were split at the query-aligned 232|233 boundary to provide paired alignments
for the mature alpha 232 and beta 96 chains. Query gaps and insertions were
handled in alignment space; integer pairing keys recorded the common source
row. Duplicate segments were filtered symmetrically to prevent independent
native chain deduplication from breaking pairs. Of 2,795 nonquery rows,
1,387 homolog pairs were retained. The older independently generated
alpha/beta alignments were not paired by synthetic header matches.

Two matched WT controls (seeds 20261126/20261127) retained eight mature protein
chains and four unguided zwitterionic PA6 L4 molecules. Stock Boltz-2 inputs
contained no structural templates or substrate contact restraints. Actual
native dataset checks verified coherent pairing across all eight protein
chains, empty template records and exact agreement with stock dummy-template
features before inference. Canonical export and matched unrestrained Amber
relaxation were handled separately from prediction completion, with all 48
site–ligand–amide combinations retained for geometric evaluation. Source
hashes, database/settings provenance and the completed local search are
preserved; the working setup is described in
[the local MSA workflow](../docs/LOCAL_MSA_WORKFLOW.md).

The current hard minimum protects canonical 113, 156, 186, 221–235 and 272–276
in all copies, with initiator M1 and structural/interface constraints added.
Lys156 is retained as a precaution informed by NylC Lys189 mutagenesis, without
claiming direct evidence of Nyl12 essentiality. Leu235 is retained within the
mature beta N-terminal 1–3 context; its individual catalytic role is unproven.
Historical protection of 275–276 likewise does not establish direct
participation in maturation chemistry. An exploratory 100-position pilot
allowed Lys156 substitutions; all 24 generated sequences changed this residue
and were excluded from advancement when the 101-position mask was restored.
Its immutable outputs are retained as diagnostics, without silent sequence
reversion or relabeling as current-policy designs.

The conditional release of the corrected 24-sequence pilot required at least
one independently accepted WT control with conserved identity/connectivity,
force-converged unrestrained relaxation and passing positive-phi geometry.
A successful seed was not used to relabel another seed as accepted. The
second paired-MSA control failed because of severe predicted interchain
overlaps; its unsuccessful minimization and native artifacts were retained,
and its endpoint was excluded from accepted relaxed-model comparisons.

## Finite substrate models, binding register and protein-template omissions

Zwitterionic PA6 L4 was used as a tractable local model of a potentially much
longer polyamide strand. The original contact-guided L4→L2+L2 models were
motivated by observed oligomeric products, including L2 and L4, rather than
by a measured L4-specific cleavage reaction. Register-independent validation
examines all three internal amides and records the selected bond separately
from fragment lengths: C54–N59 joins units 2–3 and would give L2+L2; C55–N57
joins units 3–4 and would give L3+L1 for this finite model. The latter was
selected by the recent protein-template WT controls. A finite-model product
split does not establish endo/exo selectivity for long polymers. The ongoing
relay reaction-path/PMF model retains its original central-amide mapping;
different-register results are not pooled or treated as the same pathway.

The commercial-polyamides manuscript identifies structurally variable loops
95–103 and 289–300 around the occupied site, including contributions from
neighboring subunits. Existing loop/tail-omitted crystal CIFs were audited
before reuse. Their masks differ and some retain waters/PEG; filename alone
does not establish which residues are omitted. For the bounded new design
screen, an explicit mask was derived from the validated protein-only mature
assembly: observed alpha-subunit coordinates 95–103 and 223–232 were omitted
in all four copies, preserving every other observed coordinate. All full-length
alpha 232/beta 96 protein queries and homolog MSA rows were retained. The mask
removed 486 of 9,426 protein heavy atoms; beta289–300 remained available as
template context. This frees template coordinates without deleting sequence
residues or introducing loop-length changes. No PA6 coordinates/contact
restraints or designed TS-backbone template were supplied to the predictor.

Template qualification subsequently checked the tensors supplied to the
native predictor, rather than relying only on coordinate omission and stored
alignment ranges. The installed Boltz featurizer used an offset without
clipping to the nominal alignment boundaries. Original WT and design template
tensors were identical, but all retained template features at alpha query
positions95–103, which were intended to be untemplated. Changing only the WT
alignment-record start did not change these tensors. The original six outputs
were retained as exploratory results. An isolated adapter maps both observed
alpha-subunit residues 8–94 and 104–222 and mature beta-subunit residues 1–93, leaving the
installed predictor unchanged. The actual native dataset path was checked for
all six cases against features calculated independently from source-author
residue indices. All template tensors agreed exactly with this reference and
across the six inputs, with 1,196 retained protein residues per assembly and
no template features at omitted loop/tail or ligand positions. Original
processed structures and chemistry retained their source hashes. The mapping
adapter passed 14 regression tests with 100% statement coverage. Corrected
WT and design predictions retain the same sequences, seeds, MSAs, template
coordinates and relaxation protocol; downstream outcomes are recorded
separately from input acceptance.

Legacy capped-substrate models and GAFF payloads are retained as provenance.
The available historical payload does not yet certify its exact charge-fitting
method. It is not labeled as the high-quality RESP/GAFF2 parameter release,
and capped and uncapped force-field calculations are not pooled quantitatively.

## N-terminal-relay path construction and sampling qualification

The relay hypothesis was investigated by coordinate-driven optimization from
the common force-converged window033 structure toward E–S and acylenzyme in
opposite directions. The ordered branches supplied a shared 39-node reference;
they were not refined by NEB. The TI1-like starting structure's full-minimum
qualification remains unresolved. We initialized 31 activation and 26 proton-
delivery umbrellas using two distance-difference biases per window. All windows
retained the same DFTB3/3ob-3-1 QM/MM Hamiltonian, proton inventory, movable
region and 300 K protocol. Continuous accepted segments were analyzed by MBAR
after one initial 5 ps production discard per independent history.

Projected-coordinate pileup, importance-weight concentration, distance drift
and initialization dependence were assessed separately from statistical
interval width. Corrective sampling comprised a 33-center physical-coordinate
activation grid, reciprocal-start and product-side checks, and six additional
late-activation/transition-region initialization controls. All existing data
were retained and process completion, chemistry and sampling acceptance were
recorded separately. The exact atom mapping, source construction, native bias
convention, equilibration, reweighting, ESS definitions and prospective final
uncertainty checks are detailed in [SI S1–S8](SUPPORTING_INFORMATION_DRAFT.md).
No convergence claim follows from the planned sampling duration alone.

Supplemental trajectories were incorporated into the same global MBAR fit
as soon as their native bias-energy parity, identity and production-length
checks passed. Independent starts at the same bias remained separate time
series for equilibration and uncertainty analysis. The first eligible pooled
fit used 57 main 50 ps histories and four supplemental 10 ps histories, retaining
51,700 frames after the same one-time 5 ps discard per history. Pooled estimates
were labeled point-only when a newly eligible history retained only one 5 ps
block; uncertainty bands from the original 57-history fit were not reassigned
to the pooled curve. An automatic source-signature monitor updates progress
fits and verified movie pointers; final sampling/uncertainty qualification
remains a separate analysis.

## Exploratory backbone and broad-sequence qualification — launched 2026-11-26

**Prospective design-selection policy, updated 27 November:** Boltz is
prioritized for the next matched validation screen. New designs require a
PI-approved positive-signed-phi amide under the current coherent atom mapping
and Tier A/B geometry after unrestrained relaxation. All amide/site candidates
are evaluated before selecting a passing pose; original WT control tiers are
retained separately. A WT flipped-pose pass is not required to test designs.
Loop insertion/deletion reconstruction is proposed as a separately qualified
arm, with identical edits in all four copies and explicit WT/design/Amber
maps. No loop-indel calculations have yet been performed in this study.

Matched unrestrained relaxation of the two protein-template WT controls was
subsequently performed with the previously validated four-substrate Amber
protocol. Each system retained four zwitterionic PA6 L4 molecules with released
modular RESP/GAFF2 parameters, ff19SB protein parameters, neutral Thr233
alpha-amino groups and anionic Asp272/Asp274 in all four mature beta chains.
TIP3P solvent was added in a truncated-octahedral box with 12 Å padding,
with 150 mM NaCl plus neutralizing ions. Both used a single all-atom
unrestrained periodic minimization, maximum 10,000 cycles with 1,000 initial
steepest-descent cycles, 10 Å nonbonded cutoff and RMS-gradient target
0.05 kcal mol⁻¹ Å⁻¹. Positional/contact restraints and bond constraints were
absent throughout (`ntr=0,nmropt=0,ntc=1,ntf=1`). The same source-pinned
Violet pmemd.cuda executable was used sequentially. Matching dry topology,
restart and conversion metadata restore canonical identity before paired
all-amide/four-site geometric analysis. Native completion, force convergence
and geometric acceptance are recorded separately; no MD equilibration or
QM/MM stationary-point acceptance is implied by this minimization.

A sequential two-target RFDiffusion3 qualification job was launched from an observed relay
barrier-region parent, using partial-noise amplitudes of 2 and 5 Å. These are noise
settings, not molecular-dynamics times or final structural-displacement limits.
The full eight-chain mature assembly and one physical substrate were retained.
Twenty-four pocket/entrance-loop positions in the occupied repeat were mobile;
other context and catalytic/maturation anchors were fixed. Exact four-copy
sequence ties were specified for subsequent inverse folding. This historical
local run did not impose backbone symmetry; it remains a historical local
remodeling diagnostic, outside the current baseline.

The original PA6 component IDs collided with unrelated installed CCD entries.
A reversible serialization-only relabeling retained 33 ligand heavy atoms and all
32 typed heavy-atom bonds, verified through the native feature pipeline. The source
E–S graph is declared ligand context for the actual sampled TS-like coordinates;
no stable Thr-adduct graph or fractional electronic structure is implied.
Native RFD3 removes eight protein-terminal OXT atoms. Raw output is retained;
only these source OXT coordinates may be restored to a derived export after
proper alignment and verification that their terminal anchors remain fixed.
Any unexpected missing atom fails qualification. Geometry, stereochemistry,
loop closure, clashes, sequence and ligand correspondence are checked
separately from process exit. A subsequent source-validated export adapter
also tested the native centered-BF16 serialization signature before restoring
fixed source coordinates and typed connectivity. It did not adjust mobile
atoms. The 5 Å parent01 export passed the unchanged audit; the 2 Å alternative
retained genuine packing clashes and was rejected. A second TS parent subsequently passed the historical graph/coordinate audit;
both remodeled parents remain outside the current 12-design pool.

The conditional broad pilot permits 229 of 328 precursor positions, protecting 99
including M1 and the established catalytic/maturation motif. Recurring interface
core means a heavy-atom contact within 4 Å in at least three of four copies in each of
three original parents; this is an operational protection, not proof that every
protected position is experimentally essential. Cross-scoring uses an 18-position
reactive-core neighborhood separately from global sequence compatibility.
Original and remodeled arms share the same broad mask and sampling policy.
Unguided AF3/Boltz validation requires independent WT protocol calibration,
with fresh noncovalent E–S ligand chemistry, four substrates and one diffusion
sample per seed. The first controls used query-only MSAs and empty templates;
separate follow-ups test homologous protein MSAs and protein-only crystal
templates without substrate contacts or designed-TS template guidance.
The qualified prediction policy must be matched across WT and variants.
Every internal amide is evaluated with its
own carbonyl C/O/N and acyl-side adjacent C; four ligands × four sites × three amides
yield 48 local evaluations. Legacy geometric tiers remain descriptive, and
native-Thr similarity is not a hard filter. These are exploratory templates,
not committor-qualified transition states or validated improved enzymes.

The historical broad batch comprised 24 raw sequences, six per original or
accepted remodeled template for two TS parents; it is not a 24-candidate
current nomination pool. Existing canonicalization
and modular-ligand tools generate biological and mature-local exports and
sequential Amber heavy inputs with explicit atom/residue maps. Coordinate
and catalytic-metric roundtrips are checked before publication. Matching
topology/restart/conversion metadata are used to recanonicalize Amber outputs
before subsequent RFD3/MPNN use; an engine-written chain label alone is not
an identity map. Source raw files remain immutable.

## Planned complementary amine-release assay — 2026-11-26

**Prospective protocol; no new assay data have been generated.** Purification
of WT Nyl12 plus two or three selected designs is considered feasible by the
PI. A TNBS primary-amine-release assay will first be qualified on purified WT,
using NylC precedent [R61] and the published 96-well nylon-hydrolase method
[R62, section 1.4]. Its derivatization conditions are separate from the
enzymatic hydrolysis temperature. An amine-free reaction buffer such as
phosphate will be evaluated against the existing Tris condition using paired
calibrated product-MS measurements; pH measurement temperature and enzyme
stability will be recorded. Tris, accessible protein amines and continuing
autoprocessing/proteolysis can contribute background.

Qualification will include reagent, enzyme-only, substrate-only and time-zero
controls; mature-enzyme quantification; calibration against 6-aminohexanoate
and available L2/L3/L4 standards; and matrix-spike recovery across any cleanup.
Detection/quantification limits and a usable time/concentration range must be
established before comparing designs. Independent enzyme preparations and
technical replicates will be distinguished. Background-corrected soluble
amine equivalents will be reported alongside individual calibrated L2/L3/L4
products. Neither measurement alone establishes total bond cleavage or polymer
conversion; higher oligomers lack current quantitative standards. The full
assay-development plan is retained in
`../docs/design_ranking_and_iteration_20261022.md`.

## Matched WT protonation-state comparison — queued 2026-11-21

After all twelve matched WT replicas complete, four independent 100 ns replicas
per state will be compared for NH2/ASP272/ASP274, NH2/ASP272/ASH274 and
NH3+/ASP272/ASP274. The occupied AE site will be analyzed at every saved
production frame (10 ps), using the same periodic geometry, signed Herschlag
phi, Tier thresholds and crystal-relative native-Thr diagnostics. Tier A, exclusive
Tier B, inclusive A-or-B and native intersections will be reported separately.
Equal-replica means, sample SD/range and 10 ns block stability will describe
replica spread and temporal behavior; correlated frames and the four protein
sites will not be treated as independent replicas.

The four remaining trajectories will receive the same 50 ps directed relay
screen used previously. Eight completed matched screens will be reused only
after their trajectory/topology hashes match newly analyzed inputs. Unique-frame
relay counts will be joined to all-frame native geometry, retaining direct-route
controls, assisted networks, conditional proton chemistry and failures explicitly.
Historical unmatched cohorts are excluded. This section records the queued
protocol; it does not claim that final comparison results are available.

## Leaving-group proton-relay screening — executed 2026-11-21

Sixteen completed 100 ns WT trajectories (1.6 μs) from the retained original-
and flipped-pose cohorts were screened at 50 ps spacing (32,000 coarse frames).
Completed trajectories and topology-specific proton inventories were verified
before analysis; incomplete files were excluded. Actual bonded donor hydrogens
and periodic minimum images defined paths containing up to two waters and one
declared internal residue. Each directional leg required D–A ≤3.5 Å,
H–A ≤2.5 Å and D–H–A ≥120°. Starts included Thr OG1, Nalpha, protonated
Asp oxygens and Lys156 NZ in their prepared states. Neutral-Nalpha donation
was conditional on activation, and amide-N acceptance represented a hypothesis
for reaction preparation. No proton or coordinate was changed.

Assisted routes were separated from zero-intermediate direct controls before
ranking. Catalytic d1/d2/d3, attack angle and signed Herschlag phi were measured
using the existing definitions. Twenty diverse shortlisted source frames were
audited against the recorded crystal-relative Thr displacement/torsion criteria.
The ±4 saved-frame neighborhood was rescanned: native saved spacing was 10 ps
for fourteen trajectories and 50 ps for the two original-pose trajectories.
Exact restarts, source hashes, water identities, conditional chemistry and failed
gates were retained. Coarse/refined observations of the same source frame were
deduplicated for candidate counts. These structural screens establish neither
reactive proton affinities nor a free-energy barrier.

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

## Experimental nomination and optional mechanistic follow-up

The PI may nominate candidates from sequence, assembly, structural and dynamic
evidence for expression, maturation and activity testing against WT. QM/MM PMFs
and verified PES transition states are not prerequisites for experimental
nomination. Resource-permitting chemistry may precede testing or follow measured
improvement to investigate mechanism. Computed barrier-improvement claims still
require matched, supported WT/design calculations; post hoc protocols must not
be tuned to reproduce assay outcomes.

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
D2-like, but thermal conformational differences do not require exact coordinate
symmetry. Any future optional RFD3 experiment must independently qualify
physical movement, assembly and retained catalytic context; this is not a
production-qualified arm. Symmetric sequence design does not require
identical instantaneous loop conformations or ligand occupancy. Variable-length variants are deferred
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

### PMF-first selection of saddle and ensemble design templates — proposed 2026-11-25

The current priority is a converged WT acylation free-energy profile with
connected sampling, uncertainty estimates, chemical retention/proton-inventory
analysis and checks for omitted slow coordinates. A favorable barrier value
alone will not establish mechanistic validity. A stable discrete TI1 is not a
prerequisite for the reaction or this design strategy. The following template
comparison is proposed and has not been executed.
If PES-saddle certification remains unresolved after the PMF/ensemble passes
its qualification, the single-saddle control may be deferred without calling
the ensemble templates verified PES stationary points.

Verified PES first-order saddles will provide chemically interpretable reference
geometries and single-template controls. Finite-temperature transition-region
templates will be actual cluster representatives from equilibrated, decorrelated
and reweighted frames near each relevant PMF bottleneck. We will initially retain
three to five representatives when supported by distinct conformational/contact
clusters, preserving substrate register, atom/proton inventory and source
window/block/frame. Coordinates will not be averaged. The selected PMF band is
a conditional structural ensemble; it is not automatically a dynamical TSE.
Bias-free QM/MM shooting with independent thermal momenta, declared adjacent
basins, binomial uncertainty and explicit uncommitted outcomes can test whether
representative configurations have committors near one half [R53]. Targeted PES
saddle optimization and dynamical qualification answer different questions;
unminimized ensemble parents will be retained to avoid losing thermal diversity.

A bounded, matched-budget LigandMPNN comparison will test one verified saddle
against a small ensemble arm, with identical masks, protected identities,
sampling settings and explicit alpha/beta symmetry ties. Whole candidate
sequences will be evaluated across all templates and E–S/acylenzyme/access
controls, rather than requiring identical sequences to recur independently.
The multiple-PDB option is batching, not built-in joint multistate design [R10].
LigandMPNN context/likelihoods do not encode a computed reaction barrier or
substitute for QM charge redistribution. A computed barrier-improvement claim
will require matched mutant/WT activation free energies referenced to each enzyme's own
reactant basin, plus checks for impaired later chemistry or access. Lauko et al.
combined QM-derived reactive motifs with multistate PLACER analysis [R27];
Zoi, Antoniou and Schwartz analyzed 15 uncorrelated, committor-defined TS
configurations per enzyme in a theoretical engineering study [R53]. These
precedents and related multistate work [R29] motivate our proposed comparison;
they do not demonstrate this LigandMPNN workflow or validate Nyl12 activity.
No new redesign or shooting campaign was launched by this documentation update.

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


## Complementary pKa preparation and curve audit (2026-11-20)

A formatting-only repair enabled a complete Nyl12 MCCE4-Alpha/DelPhi pilot
on the same HIS-normalized heavy-atom coordinates used for the new PROPKA
default/coupled screen. All 9,928 original heavy atoms and sixteen physical
terminal groups passed preparation checks. The native pilot used level-1
conformers, dielectric 4/80, 0.15 M salt, 298.15 K and a pH 0–14 grid.
Native completion, structural acceptance, curve/fit consistency and scientific
protonation assignment were assessed separately. We summed conformer
occupancies using native charges, distinguishing protonated neutral acids from
protonated cationic bases, checked occupancy normalization, and inspected
six-run occupancy variability. Out-of-range or inconsistent fits were not
treated as resolved microscopic pKas.

Matched-coordinate comparisons and a 6 Å titratable-residue neighborhood
screen are retained in `docs/pka_trust_and_active_site_audit_20261120.md`
and its source-hashed analysis directory. No MD/QM protonation was changed.
The earlier rejected mature-NylC result and unrun precursor inputs remain
separate records. PDBFixer repair of mature 5XYG is a new explicitly modeled
preparation; its missing-residue reconstruction and cleavage/terminal gates
must pass before downstream predictions.

## Matched sequence-mask and backbone pilot (2026-11-22)

Two fully minimized WT tetramer conformations with four PA6 L4 molecules were
used for a matched LigandMPNN comparison. The original mask had 17 alpha-chain
and two beta-chain positions; the broader nested mask added ten alpha-chain
positions. New positions required relative solvent accessibility of at least
0.25 in all four equivalent chains of both structures and no protein heavy-atom
contact to another chain within 4 Å. Exposure was calculated from protein
coordinates with RDKit-bundled FreeSASA using Lee–Richards/ProtOr radii and a
1.4 Å probe, normalized by the theoretical allowed-region maxima of
[Tien et al. (2013)](https://doi.org/10.1371/journal.pone.0080635).
These conservative design-mask choices were not fitted to activity labels.
In full precursor numbering (alpha 1–232, beta 233–328), the focused
eligible positions were 50, 51, 52, 73, 76, 101, 111, 112, 114, 115,
116, 151, 152, 153, 185, 187, 188, 271 and 273. The broader mask added
28, 37, 38, 92, 96, 97, 160, 161, 162 and 164. Beta local positions
39 and 41 therefore correspond to precursor 271 and 273. The additional
positions were drawn from 6–12 and 12–16 Å substrate shells; the exposure
and interface criteria above applied to additions, not retroactively to the
original pocket mask. Full criteria and source records are in
`docs/ligandmpnn_mask_site_rationale_20261123.md`.
Catalytic/maturation protections and chain lengths were unchanged; explicit
A–D and E–H sequence ties used equal weights. Two temperatures (0.1/0.2), two
backbones, two masks and six samples per job yielded 48 sequences; four
parent/temperature RNG seeds were paired across masks. Every FASTA was checked
for chain identity, protected and fixed positions, and sequence symmetry.
All samples and duplicate lineage were retained.

The broader arm produced 24 unique sequences from 24 samples (5–10 changes,
median 8); the focused arm produced 15 unique sequences from 24 samples
(3–8 changes, median 5.5). Sequence uniqueness does not establish structural
heterogeneity or catalytic improvement. Up to four representatives per arm
were reserved for a deterministic, stratum-balanced sequence-diversity
selection, with matched WT, K59D and V111C controls and five prediction seeds.
The same template/contact conditioning and unrestrained Amber protocol were
retained. WT predictions were reused only under scientific-input/runtime
identity; mismatched host builds were regenerated on Violet. Protein formal
charge was derived residue by residue, including the declared neutral beta
N-termini, and counterions were recalculated for each sequence rather than
assuming WT charge. No new TS geometry, activity-based threshold or change in
protonation state was introduced in this comparison.

For the follow-up, two independent 100 ns replicas were planned for each of
mask_07 seed 04, mask_03 seed 03, V111C seed 02 and WT seed 04, with four
ligands and matched neutral Thr233/ASP272/ASP274 conditions. All four assigned
sites of each starting structure have negative projected phi under consistent
scissile-amide mapping. This PI-authorized flipped-pose selection tests a possible
Thr N-terminal relay to the leaving amide nitrogen, without assuming closer
nitrogen proximity: current distances require rearrangement or possible water
mediation for transfer. Neither phi sign nor a shorter distance establishes
relay directionality or exclusive binding. Production
has no orientation-holding restraint. Initial structural differences and the
WT native-Thr deviations remain structural diagnostics, not hard rejection
criteria, following the calibrated 2026-11-24 policy. Execution
status is recorded in `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/queues/`.

## Conditional N-terminal-relay acylation paths (prepared 2026-11-23)

Follow-up paths were conditioned on qualified full-active-region numerical
minimum verification of the force-converged window 033 tetrahedral-like model.
On November 24 the PI explicitly authorized exploratory pathway generation
while full-minimum qualification remained unresolved. The original failed
minimum gate was retained; a separate source/input-hash- and force-acceptance-
checked dispatcher records this limited authorization. The source is described
as force-converged and tetrahedral-like, not a verified minimum or accepted TI1.
Two branches were prepared from its exact full-system coordinates and fixed
exterior: reverse activation toward E–S and forward N-terminal proton delivery
to the substrate leaving amide nitrogen with C–N cleavage. Composite coordinates
were d(Nalpha,H)−d(OG1,H)+d(C54,OG1) and
d(Nalpha,H)−d(N59,H)+d(C54,N59), respectively. Only the composite was restrained,
so component changes need not be synchronous. All N-terminal protons were
tracked; no proton-identity or positional restraint was added. The unchanged
98-real-atom DFTB3 QM region contained no water. Reaction biases were removed
at candidate endpoints before independent force/connectivity acceptance.
Eligible interior physical-PES maxima seed native DL-FIND dimer searches;
force convergence alone is not first-order saddle certification. Full input
hashes, numerical budgets, runtime and per-window inspection records are in
`qmmm/pa6_flipped_model30/12_nterminal_relay/`. These preparations and conditional
queues do not constitute completed reaction paths or validated kinetic barriers.


## Bounded acylenzyme hydration reconnaissance (2026-11-24)

To assess solvent access before deacylation modeling, we prepared an exploratory
QM/MM hydration protocol from the retained-product acylenzyme endpoint of the
previous direct-acylation calculation. Independent native force/electronic
checks accepted this endpoint; the associated saddle remains unverified. The
full 104,922-atom system, original 98-real-atom DFTB3 QM partition (charge −2 and
eight link atoms), and4165-atom movable 12 Å region were retained. The original
fixed exterior, composition and proton inventory were unchanged; the leaving
amine fragment remained present. No water was inserted or manually translated.

The bounded protocol comprises 1 ps NVT heating from 25 to 300 K followed by 5 ps
NVT at 300 K, using a 0.5 fs timestep, Langevin thermostat, flexible water, no
QM SHAKE and no positional or reaction-coordinate restraint. Group cards
encode the exact movable selection without mask-string truncation. Native
completion, SCC history, finite thermodynamics and unchanged periodic box/
frozen exterior are checked separately from acyl/proton connectivity. Loss
of acyl connectivity prevents automatic continuation, with all data retained.
Water remains MM during this initial hydration screen: no hydrolysis chemistry
is inferred. A participating water requires explicit QM promotion and subsequent
optimization before reaction-path calculations. This short pilot does not
establish equilibration, a free-energy barrier or the rate-limiting step.
Source/input hashes and separate process/scientific receipts are stored under
`qmmm/pa6_flipped_model30/13_deacylation_pilot/`; actual launch/completion status
is established by the native receipts, not by this prospective protocol text.

Because native belly mode omits fixed–fixed energy terms, its total-energy
offset and virial-based printed Ewald diagnostic differ from non-belly
optimization. A same-coordinate comparison at the starting acylenzyme gave
identical analytic gradients for all 4,165 movable atoms at retained precision
(maximum/RMS difference 0.0 Eh/Å; threshold 1e-9 Eh/Å), and identical printed QM
energies. We therefore do not interpret the native total-energy offset as a
chemical reaction energy. This two-evaluation diagnostic establishes force
parity at the starting geometry, independently of subsequent sampling quality.

The acyl-Thr linkage, scissile amide atoms, transferred proton and newly formed
amine headgroup were within the QM region. Original topology connectivity
therefore did not impose a classical scissile C–N or OG1–HG1 bond: Sander
removes all-QM classical bonded terms. Unchanged protein/substrate tails
retained ff19SB and RESP/GAFF2 descriptions, respectively. QM–MM Lennard-Jones
parameters and mixed terms at the eight link boundaries were retained from
the starting topology; ester/amine-specific LJ retyping was not introduced.
The leaving N and its first two alkyl carbons were QM. This avoids a complete
classical acylenzyme/amine parameter fit but retains boundary/nonbonded
approximations requiring sensitivity assessment before quantitative claims.

The retained QM–MM Lennard-Jones types are held constant across chemical states
for a continuous Hamiltonian; they are not a claim that a quantum amine remains
electronically an amide. Retyping only the product endpoint would change the
Hamiltonian and invalidate direct profile comparisons. Any nonbonded/partition
sensitivity test must use the same alternative setup for all compared states.
Fully classical MD of the acylenzyme or released amine would require suitable
validated classical parameters, which this pilot does not supply.


## Native-Thr diagnostic policy and matched velocity replicas (2026-11-24)

Catalytic geometry eligibility uses the original Tier A/B criteria. Native
Thr233 resemblance to the supplied n12c crystal-derived reference is reported
separately and is not used to reject or preferentially rank candidates. The
original OG1/Nalpha displacement and psi/chi1 deviation measurements and
reference limits are retained; native_pass and joint_pass are descriptive
legacy fields. Chemical integrity, explicit protonation, atom mapping, native
process completion and force criteria appropriate to stationary-point work
remain unchanged. This policy does not establish catalytic activity.

Both independent 100 ns velocity replicas of WT, focused mask_03, broader
mask_07 and V111C receive the same fixed-assignment analysis of all four PA6
active sites. The source-bound analysis samples every tenth 10 ps production
frame (100 ps spacing; 1,000 frames per site per replica). Per-site all-frame
denominators, Tier A/exclusive B/A-or-B occupancies, signed Herschlag phi,
d1/d2/d3, attack geometry, native-similarity scalars and first/second-half
summaries are retained without native-similarity censoring. Complete native
production steps, accepted stage/process receipts, topology/restart identity
and actual trajectory frame/time coverage must pass before analysis. Analyses
are queued on the execution host and do not interrupt MD or transfer raw
trajectories. Independent replicas are compared explicitly; nested protein
sites and correlated frames do not increase the independent replicate count.


## Water approach analysis and bounded hydration continuation (2026-11-24)

All saved configurations from the retained-product acylenzyme hydration pilot
were screened against all water oxygens using triclinic minimum-image distances.
Measured coordinates were Ow–acyl-C distance, Ow–C–Ocarbonyl angle, and contacts
with Thr233 Nalpha/OG1, the leaving amine, and both chemically equivalent
carboxylate oxygens of Asp272/Asp274. Hydrogens were retained to assess explicit
O–H···acceptor geometry. Near-attack (≤3.5 Å; 95–125°) and hydrogen-bond contact
(heavy-atom distance ≤3.5 Å; donor–H–acceptor angle ≥135°) labels are exploratory
annotations, not validated activity or reaction criteria. Continuous metrics,
water identities and closest-approach frames/times are retained. Heating and
300 K denominators are separate.

The approved additional 20 ps protocol comprises four sequential 5 ps NVT
segments at 300 K and 0.5 fs, using distinct recorded Langevin seeds and accepted
predecessor velocities without reheating. Original DFTB3 settings, zero
electronic smearing, 98 real QM atoms/eight links/charge −2, the 4,165-atom 12 Å
movable region, fixed exterior, solvent and complete amine product are retained.
No water insertion, reaction bias, positional restraint or QM-water promotion
is introduced. Each block requires native process, SCC, finite-state/temperature,
unchanged box/frozen exterior and acyl-retention checks before continuation.
Failed blocks stop explicitly. Accepted-block water analysis and matching stripped
VMD movies are refreshed through a separate source-pinned postprocess hook.
Process and scientific/sampling acceptance remain separate.

The four-block hydration continuation was launched on Violet on November 24,
2026, using a retained native/source-pinned tmux runner. Water analysis and the
portable matching-topology VMD trajectory are regenerated after accepted blocks;
these updates do not establish chemical or sampling convergence. Original pilot
outputs and failed preparation/export attempts remain preserved separately.


## Conditional hydration follow-on and assay-temperature preparation (2026-11-24)

A bounded conditional queue was submitted to wait for the accepted 26 ps
hydration endpoint and complete water analysis. If no saved water met the
exploratory near-attack annotation, it would continue from that accepted restart
for five further 5 ps blocks, reaching 51 ps total. The same annotation is
checked after every block and triggers a hold for review. Original 300 K
Hamiltonian, retained amine product, QM/MM partition, charges and frozen
exterior are unchanged. Native failures or invalid provenance stop the queue.
At submission this described scheduled sampling. By 26 November, all segments
had completed and passed their recorded native and chemical gates, yielding
51 ps total: 1 ps heating and 50 ps at 300 K. No deacylation PMF was calculated.

A further bounded continuation to 101 ps total was approved on 26 November:
ten sequential 5 ps segments on one Violet CPU, starting from the accepted
51 ps restart with its velocities and without reheating. Native Hamiltonian,
partition, bath, mobile region and retained product are unchanged. Water
geometry and the stripped trajectory are updated after accepted blocks;
a near-attack candidate or a failed process/SCC/chemical/boundary gate stops
automatic continuation for review. All waters remain MM in this reconnaissance.
A candidate water must be included explicitly in a revised QM region and the
hydrated acylenzyme reoptimized before modeling water attack/proton transfer.
The 51 ps endpoint has neutral NH2 groups at both the Thr233 N-terminus and
the amine product; HG1 is attached to the latter. An earlier derived annotation
calling the N-terminus protonated was identified as stale and is corrected in
the continuation analysis. These hydration simulations do not substitute for
a converged deacylation pathway or free-energy profile.

Two additional stochastic hydration branches were initiated from accepted
26 ps and 41 ps checkpoints, retaining their coordinates and velocities and
using distinct recorded Langevin seeds. Each branch is bounded at 25 ps of
new sampling in five 5 ps NVT segments. Hamiltonian, proton inventory, QM/MM
partition, movable region and retained product match the main continuation.
Historical trajectory prefixes are referenced as read-only chronological
context; they are not counted again as independent replica data. New sampling
statistics are identified by branch-local production segments. All three runs
are exploratory water-access calculations, with the same candidate/failure
stop criteria; they are not independent equilibrated estimates of water-entry
kinetics or deacylation free energies.

Separately, both completed WT noncovalent E–S 100 ns endpoints were staged for
a prepared-only 348.15 K branch: 200 ps NVT temperature ramp from 300 K, 1 ns
NPT density adaptation and 2 ns unrestrained NPT, followed by scientific review.
Full-system velocities, box, ff19SB/TIP3P/PA6 parameter release and neutral
Thr-NH2/ASP272/ASP274 protonation are preserved. The 130,815-atom four-substrate
WT MD system differs from the 104,922-atom single-substrate legacy QM/MM system;
these preparations cannot share atom indices or direct energy references.
No thermal MD has been launched and no 348.15 K equilibration is claimed.
Subsequent mechanistic simulations require state-specific full-system thermal
preparation and explicit QM mapping. Covalent acylenzyme remains QM-treated.

Temperature-specific transition-state-theory illustrations use
k=(kBT/h)exp(−ΔG‡/RT), with κ=1. At 348.15 K, 25 and 28 kcal/mol imply
1.47e−3 and 1.92e−5 s−1, respectively; at 300 K the same hypothetical 25 kcal/mol
barrier gives 3.84e−6 s−1. These do not establish how ΔG‡ changes with temperature.
The crude-lysate assay used 75 °C for 6 h in nominal-pH-8 Tris, whereas the
purified commercial-polyamide assay used 65 °C for 24 h in phosphate. The
Tris pH-adjustment temperature is unknown; its temperature dependence prevents
assuming identical assay and room-temperature proton activities. Product
conversion is not converted into kcat or a unique mechanistic barrier.

For new relay PMFs, both proton-transfer and acyl-bond progress must be resolved,
with continuous proton inventories and forward/reverse path checks. A composite
coordinate that hides another slow coordinate can produce misleading projected
barriers (Rosta et al., 2009). Existing umbrella/MBAR infrastructure remains the
planned method. PathCV/OPES is a documented methodological option if coordinate
hysteresis persists (Rivas-Fernández et al., 2026), not an installed or executed
replacement campaign. Protonation-microstate population penalties and sampling
uncertainty must accompany comparisons of barriers from different states.


## Relay connectivity and prospective PMF qualification (2026-11-24)

Exploratory reverse activation and forward N-terminal-relay acylation share the
same force-converged window 033 restart, atom/proton inventory, Hamiltonian and
frozen exterior. The reverse branch must be read from its released E–S endpoint
back to window 033 before joining the forward branch. A dedicated postprocessing
audit was commissioned to verify accepted-window lineage, topology/mapping,
chemical retention, force acceptance, tracked and unlabeled proton ownership,
component-distance continuity, and unbiased endpoint identities. Initial
window 033 minimum qualification remains unresolved. A stitched geometric path
or force-converged endpoint does not certify a minimum-energy path, TS or PMF.

The candidate acyl coordinate is r(C,Nleave)−r(C,OG1). Proton-transfer and
heavy-atom progress will be resolved independently when choosing umbrella
coordinates. Existing all-QM-N/O hydrogen-pool coordination observables diagnose
label exchanges and N-terminal proton residence without adding proton-identity
restraints. These nonlinear observables are not assumed to be implemented as
bias forces. Separate local free-energy surfaces, if used, require a connected
overlap region and a common statistical normalization; independently zeroed
one-dimensional profiles cannot be added to obtain an overall barrier.

A new PMF campaign is conditional on chemical connectivity and a supported,
validated reaction-coordinate implementation. The preregistered sampling gates
in `docs/qmmm_umbrella_sampling_and_pmf_protocol_20261109.md` require connected
overlap, correlation-aware effective sample sizes, block/independent-start
stability, uncertainty and estimator checks, and stationary orthogonal
observables. A staged unbiased committor test and sparse matched embedded
DFT/MM checks are planned for a candidate bottleneck. Neither the expected
biochemical barrier range nor a smooth reconstructed PMF is a success criterion.
Reactant-state normalization must include conformational preorganization:
a path from an unusually relay-ready E–S pose yields a conditional barrier
unless equilibration connects that pose to the dominant bound ensemble.
Nalpha–leaving-N geometry and substrate register are candidate slow-coordinate
diagnostics; neither a selected snapshot nor Tier occupancy alone determines
a preorganization free-energy correction.
These qualification steps are prospective, not completed results.


## Planned matched comparison of acylation mechanisms (2026-11-25)

The direct Thr233-OG1-to-substrate-N proton-transfer mechanism (A) and the
Thr233 alpha-amino-mediated proton relay (B) will be evaluated with matched
chemical/proton inventories, substrate binding register/direction, force field,
QM/MM partition and Hamiltonian, temperature, movable/frozen environment,
equilibration, bias conventions and sampling-convergence criteria. Actual
identity hashes, rather than atom counts alone, establish comparability.
Reaction coordinates/window centers may differ to cover the channels, but
shared component-distance/proton-coordination observables, a common free-energy
measure and one E–S basin reference are required for quantitative comparison.
Exact sampling biases will be removed with MBAR; WHAM provides an estimator
cross-check where appropriate. Correlation-aware block/independent-start tests
and bootstrap uncertainty, including reference uncertainty, are required.

The completed direct-seeded dataset will first be reanalyzed for reusable
support. Its preliminary `q_acyl` marginal mixed proton-transfer channels and
was nonconverged; it is not a pure-direct barrier and cannot be subtracted
from a separately normalized relay path-arclength profile. Channel switching
will be retained and reported rather than suppressed with proton-identity
walls. Where a common multidimensional reconstruction is supported, channel
bottlenecks will be compared in that common normalization. Conformational
preorganization or additional sampling-restraint costs must be included or
explicitly identified as limitations of conditional barriers. This section is
a prospective protocol; no new matched direct/relay result is claimed.

## Provenance and nomenclature of window 033 (2026-11-25)

The window 033 relay/addition-like seed was drawn from
`window_033/03_production_block_05_5ps/trajectory.nc`, frame 46 (one-based),
recorded time 28.299999 ps, SHA-256
`98ea1bac883503741de3236126786560233094cad73d15a20682e10923287845`.
The initial extracted-coordinate hash is distinct from the subsequently
optimized restart hash. Source identities, partition/mobility hashes, native
SCC history, initial/final geometries and independent gradient checks were
retained. The umbrella bias was removed before minimization. Candidate
retention, force convergence, chemical connectivity, full-Hessian minimum
qualification and thermal lifetime are separate outcomes; the directory name
`ti1_window033` is a historical identifier, not chemical validation.
The full same-protocol comparison and limitations are documented in
`docs/neb_us_nterminal_relay_release_20261123.md`.


## Live SCC guard interpretation during relay follow-up (2026-11-25)

Intermediate SCC nonconvergence/reset messages were distinguished from final
unresolved electronic evaluations. An initial conservative guard incorrectly
aborted at a recoverable70-cycle warning; the accepted window 033 history itself
contained70→96-cycle recoveries. Original aborted attempts are retained. The
corrected isolated guard permits bounded recovery through500actualSCCiterations,
rejects fatal/unresolved final evaluations, and leaves electronic temperature,
SCF/density tolerance and final force criteria unchanged. The repaired downhill
attempt starts from its original independently validated displacement. Guard
completion or a recovered SCC trial does not establish chemical retention or
transition-state index.
## Pointwise uncertainty and sampling support for relay PMFs

For cumulative relay PMF checkpoints, the cross-window MBAR analysis,
three-coordinate path projection, histogram edges, original equilibration
discards and statistical thinning were retained unchanged. PyMBAR 4.0.3
histogram covariance estimates provided asymptotic standard errors for
free-energy differences relative to a fixed sampled E–S reference bin.
The displayed pointwise 95% intervals were calculated from 200 bootstrap
replicates. In each replicate, completed 5 ps production blocks were resampled
with replacement independently within each umbrella window, retaining each
block's original selected configurations; MBAR was then refitted to the
resampled cross-evaluated bias energies. The histogram grid and reference bin
were fixed across replicates. Thus uncertainty in the MBAR normalization was
re-estimated rather than holding frame weights fixed. Replicate bin/reference
support was recorded explicitly; percentile intervals were displayed only for
bins supported by every replicate. These are pointwise percentile intervals,
not simultaneous bands or uncertainty estimates for an optimized maximum.

Support diagnostics included the number of selected configurations, the number
of contributing umbrella windows, and the bin-local weight effective sample
size, `(sum_i w_i)^2 / sum_i w_i^2`, with sums restricted to configurations in
the bin. The latter diagnoses weight concentration and is distinct from the
number of statistically thinned configurations. Statistical intervals remain
conditional on the sampled configurations and estimated correlation times;
four completed blocks per window cannot exclude slower structural or chemical
relaxation. Sampling-convergence and stationary-point acceptance were assessed
separately. The implementation follows the
[PyMBAR histogram-FES uncertainty interface](https://pymbar.readthedocs.io/en/stable/fes_with_pymbar.html)
and the MBAR framework of
[Shirts and Chodera (2008)](https://doi.org/10.1063/1.2978177).

## Continuous-window equilibration and correlation analysis

### Path-coordinate and configurational-support diagnostics

Nearest-segment projection onto the three-distance-difference reference path
was tested for probability accumulation at its vertices. We quantified both
the raw frame fraction and MBAR-weighted fraction assigned exactly to each
vertex, and calculated Kish weight concentration within every PMF bin. This
quantity measures importance-weight concentration, not temporal independence.
A differentiable path-coordinate sensitivity and physical two-dimensional
histograms used the same global MBAR weights. The activation analysis retained
both q_act and q_acyl because a marginal in q_act alone merges configurations
from different chemical stages. No kernel smoothing or removal of inconvenient
proton configurations was used to obtain these diagnostics. Confidence
intervals remain conditional on sampled modes and cannot validate a coordinate
projection or compensate for unsampled configurations.

**Corrective sampling protocol specified 27 November; execution/results are
recorded separately.** Supplementary activation umbrellas span q_acyl = −1.25,
−1.50 and −1.75 Å and q_act = −1.0 to +1.0 Å in 0.2 Å increments (33 states).
The two native biases retain k = 25 kcal mol⁻¹ Å⁻² with E = k(ξ−ξ₀)²;
no additional q_del or proton-identity restraint is introduced. Accepted
same-Hamiltonian coordinates initialize 5 ps equilibration and 30 ps production
at 300 K in 5 ps segments. The existing reciprocal-start and product-tail
diagnostics are continued from their 10 ps checkpoints to 30 ps without
reinitializing velocities. Every independent lineage retains its own initial
5 ps production discard and bootstrap resampling groups; independent starts
are not concatenated into an artificial continuous trajectory. All umbrella
biases must be cross-evaluated and checked against native restraint energies
before pooling. At 10 ps production only one retained 5 ps block exists, so the
specified 5 ps block-bootstrap interval is unresolved; eligibility begins at
15 ps production, with limited precision from only two retained blocks.

### Revised primary analysis: fixed burn-in and correlated-block uncertainty

The following protocol supersedes the latest-of-thirteen truncation rule for
new relay PMF publications; the earlier estimates remain archived below as
analysis-sensitivity results. Statistical truncation procedures apply to QM/MM
and classical trajectories alike, but their reliability depends on sampled
relaxation times and transitions rather than on the Hamiltonian or trajectory
duration alone [R55–R58]. Short trajectories can provide noisy autocorrelation
estimates and cannot reveal unvisited states. No universal equilibration or
production duration was assumed.

Each umbrella window had a separate 1 ps heating and 5 ps equilibration stage.
These stages were excluded from all PMFs. Uninterrupted production segments
were concatenated in timestamp order. We provisionally discarded an additional
5 ps of production once per window, specified before examining the new PMF,
and retained every subsequent configuration through the latest common complete
sampling horizon. This transparent analysis choice was motivated by the
time-series audit and is not a demonstrated equilibration time for every
window. No late excursion or undesired chemical configuration was selectively
removed. Fixed production discards of 0, 2.5 and 10 ps and early/late halves of
the post-discard interval were evaluated as sensitivity checks.

All retained frames contributed to the MBAR point estimate using their actual
window counts. Consecutive frames were not treated as independent for error
estimation: 200 bootstrap data sets resampled nonoverlapping contiguous 5 ps
chunks with replacement separately within each window, retained the full
chunk chronology and observables, and refitted MBAR for every data set.
Resampling 2.5 ps chunks assessed block-length sensitivity. The Hamiltonian,
cross-evaluated umbrella potentials, path projection, histogram edges and
fixed E–S reference bin were identical across these analyses. Pointwise 95%
percentile intervals were reported only for bins supported in every bootstrap
refit; missing reference support or a disconnected overlap graph invalidated
that draw rather than being silently omitted. No IID analytical confidence
interval was assigned to the correlated all-frame estimate. Weight-based ESS
was labeled as weight concentration, without implying temporal independence.

Correlation estimates for the three reaction coordinates and ten additional
chemical diagnostics were retained for diagnosis only, alongside 5 ps means,
temporal PMF comparisons, overlap and per-bin support. A few resampled chunks
cannot establish that the block length exceeds all relevant relaxation times.
These intervals therefore describe conditional statistical precision, not
systematic uncertainty from initial conformations, hidden slow coordinates,
QM/MM model error or incomplete equilibration [R56–R58]. Persistent temporal
drift calls for further equilibration/sampling or independent initialization,
rather than increasingly late automatic discards. Statistical criteria were
not adjusted to favor a lower activation barrier.

Recent QM/MM precedents illustrate complementary checks rather than a fixed
sampling prescription. Jäckering et al. combined an adaptive string, explicit
substrate-orientation coordinates, umbrella integration and a TS-region
confidence criterion [R59]. Zheng et al. analyzed the last 15 ps of 30 ps/window
QM/MM sampling, assessed profile convergence and performed committor analysis
[R60]. Those choices do not establish adequate sampling for Nyl12/PA6.

### Matched analysis of direct-NEB-seeded and N-terminal-relay sampling

The fixed-discard and contiguous-block MBAR procedure was applied to both
campaigns through a shared statistical implementation. All 63 direct-seeded
windows had 30 ps production. A matched 25 ps analysis retained 400 frames
per window after the single additional 5 ps discard; a separate full-data
30 ps analysis retained 500 per window. Native heating and equilibration
were excluded independently. Actual q_acyl/q_H umbrella potentials were
retained for the direct data; the relay's different coordinate combinations
and force constants were not substituted during reweighting.

The archived direct PMF was a marginal along q_acyl, whereas the relay PMF
used a three-dimensional reaction-coordinate path. To apply consistent density
conventions, distances were recovered from the native direct trajectories for
five reactive atoms, using periodic minimum-image coordinates. Source hashes,
times, all four previously recorded distances, q_H = q_act + q_del, and
archived cross-evaluated bias energies were verified for every production
segment. No frames were removed by proton ownership or catalytic geometry.

Both data sets were projected in (q_acyl, q_act, q_del), with equal Angstrom
weighting and the same extended-terminal-tangent arclength algorithm. Each
retained its physical geometrical reference path; the direct reference used
the original 24 NEB images, verified against their recorded restart hashes.
Normalized arclength therefore parameterizes different reference curves and
does not identify identical molecular configurations between pathways. The
histogram grid width was 0.025, with a fixed reference at s = 0.0125. The direct
grid was extended by one terminal bin to preserve fluctuations beyond the
relay plotting range; interior bins and reference were unchanged. Pointwise
95% intervals, block-length/discard sensitivities and temporal-half profiles
were computed identically. Direct-seeded proton visits to the N-terminus were
reported explicitly rather than censored to enforce a mechanistic label.
These projected densities remain conditional on sampled configurations and
are not, by themselves, a comparison of mechanism-specific kinetic rates.

### Initial structures and initialization-dependence controls

The direct campaign was initialized from the original 24-image NEB path.
For the N-terminal relay, separate coordinate-driving branches were initiated
from the same relaxed window033 configuration: one toward E–S and one toward
acylenzyme. The E–S-directed branch was reversed in order and joined to the
acylenzyme-directed branch at the identical shared restart, included once in
the 39-node geometrical reference. Intermediate CV centers were interpolated;
initial atomic coordinates were selected from accepted restart structures.
The native umbrellas restrained distance-difference combinations, rather than
RMSD to an NEB image or a geometrical path. No new NEB refinement was required
as an initialization prerequisite.

All umbrella distributions were analyzed together with one global MBAR
normalization. Independently zeroed partial PMFs were not spliced or added.
Different procedures for generating initial configurations do not invalidate
an equilibrium comparison, provided that each window samples its intended
biased equilibrium distribution with adequate overlap. A shared geometrical
join alone does not establish equilibration of protein, solvent or chemical
degrees of freedom. The retained reference curves also define route-specific
plotted projections, which were documented separately from native biases.

**Planned initialization-sensitivity control, not yet completed:** after the
extended-sampling analysis identifies persistent drift or weak overlap, a
small number of windows near the relay junction, barrier or product basin will
be initialized independently from configurations approached from opposite
sides of the reaction. For each tested window, both initial conditions will
use the exact same umbrella center and potential, topology, proton inventory,
QM/MM Hamiltonian, movable region, bath and sampling/analysis protocol.
Chemical-distance/proton-owner distributions, environmental hydrogen bonding,
water occupancy and reweighted estimates will be compared after the same
initial-discard procedure. Persistent dependence on the starting side will be
treated as unresolved sampling, even if marginal CV histograms overlap or
bootstrap intervals are narrow. Such controls test initialization hysteresis;
they are not additional NEB calculations or reverse-rate measurements.

Both existing campaigns are being extended toward 60 ps production per window
in retained 5 ps segments, without reheating or changing the native biases.
Direct sampling is assigned to Violet and relay sampling remains on Azalea.
Equal nominal duration is a matched-checkpoint diagnostic, not a declaration
of convergence. The fixed-discard, all-frame estimator and correlation-aware
uncertainty/sensitivity criteria are applied consistently to both campaigns;
convergence may require unequal final sampling lengths.

### Archived automatic-truncation sensitivity analysis

For the revised relay analysis, accepted production records were concatenated
in physical-time order for each umbrella window, verifying a continuous output
schedule across segment boundaries. The same three reaction coordinates and
ten additional distance/coordination diagnostics were used as in the initial
analysis. A candidate truncation point was estimated for each nonconstant
observable using PyMBAR 4.0.3's implementation of Chodera's automated
truncation-point selection procedure [R55]. Candidate discard times were
ranked by the estimated effective sample count, approximately
`N_retained(t0) / g(t0)`, where `g` is statistical inefficiency estimated from
the observable's time autocorrelation. The exact installed implementation
and its index convention were retained in the source-bound audit.

We used a **study-specific aggregation rule**: the largest individual discard
index among the thirteen observables defined one shared suffix per window.
Statistical inefficiencies were then re-estimated over that suffix, and their
largest value set the conservative integer subsampling interval, `ceil(g)`.
The latest-of-thirteen rule and common-suffix maximum-g selection are
additional choices in this study, not procedures prescribed or validated by
Chodera's single-observable method. The resulting one discard and correlation
estimate per window replaced independent treatment of each 5 ps segment;
no production segments were reclassified as new simulations.

The unbiased histogram PMF was refitted with identical physical Hamiltonian,
umbrella potentials, path projection, histogram edges and fixed sampled E–S
reference. Matched analytical PyMBAR pointwise intervals were used to compare
selection policies. Supplementary uncertainty was assessed with 200
window-stratified resamples of the selected, whole-window-decorrelated
configurations, refitting MBAR in every replicate. Those intervals assume
effective independence after thinning and condition on the estimated burn-in
and correlation time; they exclude uncertainty in those estimates themselves.
They are distinguished from the archived complete-five-ps block bootstrap.
Automatic truncation was treated as a bias–variance heuristic, not as a test
that equilibrium or adequate configurational sampling had been achieved.
Clark, Cole and Michel [R56] demonstrate that correlation-estimation noise
can produce late, variable truncation points and distinguish this problem
from detecting insufficient sampling. Accordingly, the time-series audit and
burn-in sensitivity were recorded separately; temporal stability remained a
required acceptance check. Statistical intervals remain conditional on the sampled
configurations and analysis choices; broader sampling-quality limitations
were interpreted following Grossfield et al. [R57].

Production continuations were prepared from the existing 30 ps checkpoints
to a bounded total of 60 ps/window in 5 ps segments. Restart coordinates,
velocities, thermostat settings, restraints, QM/MM model and frozen environment
were retained without reheating. Stage-specific Langevin seeds were distinct
and recorded. Original accepted-stage identities remained tied to the original
execution plan; the new stage identities were tied to a separately reviewed,
hashed continuation plan.

## Diagnostic audit of automatic truncation choices

For the eight windows with the largest automatic production discards, all
thirteen observables, the driving observable's raw and 0.5 ps display-averaged
time series, and the effective-sample objective as a function of discard were
retained. Objective curves were checked against the installed PyMBAR 4.0.3
API. Both the default approximate lag summation and the full lag summation
were evaluated. Fixed initial-discard sensitivity was assessed diagnostically
using the same maximum-g conservative subsampling rule; these counts did not
replace the source PMFs or establish an accepted equilibration time. This
distinction between truncation selection and a demonstration of sufficient
equilibrium sampling follows Chodera [R55], the methodological reassessment
by Clark, Cole and Michel [R56], and the sampling-quality guidance of
Grossfield et al. [R57]. Full bibliographic entries and DOIs are retained in
`manuscript/REFERENCES.md` and `manuscript/recent_design_and_relay_references.bib`.


### Diagnostic checks on projected PMF coordinates (26 November checkpoint)

In addition to common bin grids and bootstrap analysis, we audited the scalar
projection of the same retained three-CV samples onto the reference polyline.
We counted exact interior-vertex projections, reconstructed the saved histogram,
and examined raw counts, Kish weight ESS, bootstrap support and temporal-half
estimates per endpoint bin. Kish weight ESS measures weight concentration and
is not a count of decorrelated frames. Interior nearest-segment clipping can
map finite regions of CV space to individual path vertices; reference-bin
subtraction also pins that bin's uncertainty to zero in each bootstrap refit.
Consequently, narrow conditional intervals do not exclude coordinate artifacts
or unsampled conformations. Terminal tangent-extrapolation bins were identified
explicitly rather than interpreted as new reference minima. The bounded binning-sensitivity calculation retained the same physical
biases, frames and MBAR normalization; the original profile was reproduced
exactly before reporting half-bin offsets, doubled bin widths and weighted
vertex contributions. All histograms used the same original reference-bin
density rather than independently shifting each curve to its own minimum.
No smoothing or alternate projection silently replaces the primary comparison.

For PMF reporting, each bin is additionally screened for raw count, Kish weight
ESS, largest normalized frame weight and bootstrap support; bins beyond either
reference endpoint are marked, and their support is assessed separately from
the reference-basin density. Projection pileup is quantified as reweighted
corner mass per bin, not inferred from raw frame counts alone. Time-block and
independent-start/replica stability, adjacent-window overlap, and declared
bin-edge and projection sensitivities are reported alongside the profile.
Pointwise bootstrap intervals remain conditional on the sampled states,
coordinate/reference, and resampling units and are not total uncertainty or
evidence by themselves for equilibration or endpoint support.



### Matched finite-bin WHAM diagnostic (28 November 2026)

We compared the displayed 82,000-frame, 100-history pooled MBAR point estimate
with finite-bin WHAM using the same post-discard observations, native bias
functions, final s bins and E–S reference. Because the umbrella biases depend
on different pairs drawn from three physical distance-difference coordinates,
a one-dimensional harmonic approximation in projected s was not used. We
formed uniform three-dimensional physical-coordinate cells of width 0.1,
0.05, 0.025 and 0.0125 Å and evaluated all native bias functions at occupied
cell centers. Repeating each cell-center bias column according to its sample
multiplicity solves the finite-bin WHAM equations using the existing PyMBAR
numerical engine [R65,R66]. This is a binning-sensitivity test, not an
independent software implementation. Original per-history counts were retained;
no correlation correction was introduced into only one estimator.

The primary comparison distributes each WHAM cell mass among its original
frames and uses those frames' original s values, preserving within-cell
empirical observable proportions. A separate cell-center projection tests the
additional projection discretization. Final histogram widths and energy zero
are unchanged. A half-cell origin shift tests the finest grid. Ten analytic
tests covered known probability masses, state-energy gauge shifts, duplicate
bias splitting and probability conservation; core coverage was 98%. A 0.1
kcal/mol finest-grid/origin comparison was declared as a numerical diagnostic,
not a sampling-convergence gate. These point comparisons carry no new
statistical confidence intervals. All source/input hashes and failed gates are
retained; no dynamics or molecular Hamiltonian was changed.

## Relay supporting-information record (27 November 2026)

The [supporting-information draft](SUPPORTING_INFORMATION_DRAFT.md) consolidates
the actual coordinate-driven relay seed construction, atom mapping, native
two-CV biases, fixed-discard/global-MBAR analysis, ESS definitions and the
45 ps support/drift diagnostics with source and figure links. It distinguishes
the running 33-state activation supplement and four continuation checks from
validated convergence; the original dated analyses above remain historical
evidence. Initial TI1 minimum qualification and TS validation remain unresolved.
