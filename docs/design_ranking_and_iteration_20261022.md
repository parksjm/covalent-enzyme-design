# Design ranking, iteration and symmetry — 2026-10-22

Status: living evaluation protocol; dated sections distinguish completed pilots
from proposed follow-up. This extends the October 16 focused-redesign protocol in
`/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/NYL12_FOCUSED_REDESIGN_WORKFLOW_20261016.md`
and the current backbone/access strategy. Freeze numerical settings before
ranking new designs; calibrate with matched WT and available mutant controls.

## User observations and interpretation

The user reports lower prediction confidence and increased effective strain
for many recent single mutants relative to WT. Retain this as a project
observation, with its underlying datasets/protocols, rather than assume every
mutation is deleterious or every departure from WT is undesirable. Coordinated
sequence redesign can accommodate interactions that isolated substitutions
disrupt. That is a design hypothesis to test for Nyl12, not an established
benefit for every LigandMPNN output.

LigandMPNN conditions sequence generation on a supplied protein structure and
local atomic context; it does not require an MSA as a design input [R10].
Published experimental binding successes support testing this approach, but
do not guarantee folding, affinity or catalysis for our focused redesigns.
Folding confidence, thermodynamic folding stability and binding affinity are
different quantities. LigandMPNN sequence likelihood is also distinct from
AF3 confidence. A favorable likelihood is not a predicted dissociation constant.

Independent single-sequence/no-MSA prediction is a useful planned test of
sequence–structure consistency. Use the same prediction settings for WT,
single-mutant controls and coordinated designs. If an MSA-assisted arm is used,
declare its construction, depth, templates and pairing separately; do not rank
WT-with-MSA directly against designs-without-MSA as a controlled comparison.
Record templates, seeds, model version, sampling budget and guidance in each
arm. Single-sequence AF3 performance for these designs remains untested here.

## Broader sequence and backbone redesign — revised plan, 2026-11-26

**Current PI steering supersedes the restrictive prospective masks, contact
guidance and single-register selection below.** The completed five-site,
24-sample pilot remains a technical baseline with four retained sequences;
it is not the limit of the next design search. The PI prioritizes general
PA6 hydrolysis, accepts any chemically consistent amide register and requests
consideration of TS-conditioned RFDiffusion3 partial diffusion before broad
LigandMPNN design. SolubleMPNN remains a separate later surface-design layer.
This section is a revised plan, not a claim that a new diffusion, AF3 or MD
campaign has launched. Keep the current relay US and targeted checks running.

### Required pose and bounded loop-length follow-up — PI steering, 2026-11-27

**Correction:** “flipped = negative phi” was an assistant interpretation error.
The PI explicitly approves the positive-phi Boltz pose. Numerical measurements
and historical cohorts are not relabeled. New pose gates require explicit
`target_phi_sign=+1` against the retained reference coordinates; negative
functions remain only as legacy diagnostics.

Prioritize Boltz for the next design-model screen on the evidence from the
matched protein-template controls and unrestrained Amber relaxation. AF3 is
an optional diagnostic, not a required parallel predictor for every design.
Keep the selected protein-MSA/template policy and unrestrained relaxation
matched across WT and designs; PA6 placement remains unguided. This choice
is specific to the tested workflow and does not establish universal predictor
superiority.

**Advancement of new designs requires the flipped local substrate pose:**
positive signed Herschlag phi with the current ordered, graph-consistent
scissile-amide definition, together with Tier A/B catalytic geometry after
unrestrained relaxation. A negative-phi Tier A site does not satisfy this
combined requirement. Retain the positive-phi WT controls and their original
legacy tiers as evidence; WT need not meet the design-specific pose target
to serve as a control. Do not postpone all design predictions until a flipped
WT control passes. Current fixed-length samples may fail the combined gate.

Apply pose eligibility to **all 48** site–ligand–amide combinations before
selecting the best eligible pair or one-to-one assignment. Filtering only
the legacy best assignment could conceal a flipped A/B alternative. Require
at least one eligible site for a model to advance, report all four sites and
prediction seeds, and prefer reproducible occupancy over an isolated favorable
snapshot. Refresh sign checks against the mapped relay E–S reference whenever
atom ordering or ligand serialization changes. No contact restraints or tier
threshold changes are introduced to obtain passing models. LigandMPNN TS
compatibility and later MD/QM/MM validation remain separate requirements.
Before release of prediction-based design selections, qualify the combined
pose/tier gate with explicit negative-only rejection and hidden-positive-A/B
fixtures. The historical control-analysis tiers themselves remain unchanged.

**Design identifiers and structural parents:** `relay_NN` labels a generated
LigandMPNN sequence, not a reaction-path window or frame. The 24-design pilot
sampled six sequences per backbone: two selected relay-US representatives and
their two RFD3-remodeled counterparts. `relay_23` is sample5 on remodeled
parent03 (106 WT substitutions); its source representative is
delivery_full_007/block05/frame95, s = 0.671934. `relay_02` is sample2 on
original parent01 (97 substitutions). Prediction seeds20261126/20261127 are
another sampling level. Stable job IDs are preserved; the source-checked
`design_identity_catalog_20261127.json` records all24 sequence lineages.

**Prospective refinement, 27 November:** use RFD3 as an optional, bounded
backbone sampler. Compare unchanged TS-like parents with a small number of
gentle whole-assembly partial-diffusion outputs, including the active-site
loops. Preserve a minimal catalytic/substrate heavy-atom motif, required
maturation and assembly anchors, chain lengths and explicit four-copy sequence
ties. Broad sequence freedom remains allowed. Input noise is not final RMSD
or evidence of improved catalysis; measure resulting backbone displacement,
motif retention, peptide closure, clashes and interfaces. The qualified local
24-residue diffusion masks do not establish readiness for a globally mobile
tetramer. Global partial diffusion requires its own native input/output audit.
Keep a no-RFD3 baseline and the existing bounded sample budget; do not add
indel or de novo scaffold branches before the current screen is assessed.

RFD3 has experimental de novo enzyme-design evidence in Butcher et al.,
Science (2026), DOI10.1126/science.aeh2982 [R31]. That evidence supports its
use as a candidate generator, not a claim that this specific Nyl12 remodeling
protocol improves activity. The study does not replace our controls or
activation-free-energy comparisons. Preferential TS stabilization relative
to E–S is the objective; stronger substrate binding alone is insufficient.

Distinguish structural file roles. RFD3 partial diffusion still needs a
starting structure from the mechanistic ensemble; a CIF used for that purpose
is not an external crystal template for Boltz/AF3. The preferred next E–S
validation screen uses MSAs without structural templates or substrate contact
restraints. Qualify WT first, because prior MSA-only controls failed pocket
geometry. Do not feed the designed TS backbone back as a prediction template
and call its recovery independent validation. E–S predictions need not retain
every TS distortion; score fold/assembly and catalytic approach separately
from TS compatibility. Retain template-assisted predictions as a labeled
sensitivity/control condition where useful, without mixing protocols in a
design ranking. The currently running corrected six-model campaign continues.

A focused next MSA test is a full-precursor homolog alignment split at the
query-aligned cleavage boundary into mature alpha and beta alignments, retaining
same-homolog row provenance and explicit pairing keys. This is an untested
hypothesis for improving interface prediction. Keep the actual prediction as
eight mature chains; a precursor MSA is not a change to the cleavage state.
Do not pair independently generated alpha/beta rows by matching
`homolog_NNNNN` labels. Use fresh output directories and verify the native
record has no structural templates and template features are absent; do not
reuse processed records from a templated run. Cache searches by exact query,
database version and settings, and reuse equivalent chains. New local search
setup/downloads require a separate dependency/resource decision.

**Loop insertions/deletions are a prospective, separate qualification arm.**
The existing fixed-length pilot continues unchanged. Length changes can
reshape the substrate channel and donor/acceptor geometry beyond side-chain
substitution, but added flexibility can impair catalytic preorganization.
The PI specifies WT **95–103 first**, then **289–300** if another loop is
needed. In the occupied A/E pocket these are neighboring D-alpha and F-beta
loops, respectively; assess the full assembly. Use structure, existing homolog
alignments, packing and recurrent failure geometry to qualify each edit.
Before 30 November, first reuse audited loop-coordinate-omitted templates
with full-length protein queries; do not start true indels before reviewing
the current bounded model screen. Protect the
cleavage junction, free catalytic Thr terminus, catalytic/oxyanion-hole anchors,
conserved structural elements and required inter-repeat contacts. Avoid
terminal-ion capture or obstruction of longer polyamide strands.

**Input-role audit, 27 November:** both qualified RFD3 parent CIFs retain
complete alpha chains 1–232 and beta chains 233–328 in canonical exports.
Both loop regions are present. Their fixed-length partial-diffusion inputs
move 24 other pocket residues; alpha 95–103 and beta 289–300 (mature-local
57–68) are fixed in every copy. Thus the current two remodeled parents do
not test backbone remodeling of these loops. For a subsequent fixed-length
95–103 experiment, retain the full starting CIF and explicitly free the
appropriate loop backbone atoms in the native fixed/mobile mask, while
preserving the TS motif and required anchors. Audit all four copies and
numbering maps before qualification. Removing loop coordinates alone is not
this operation. Protein-prediction templates separately omit alpha 95–103
and tail 223–232 while keeping full query sequences; true sequence deletions
require a separately qualified RFD3 contig/inpainting specification and
regenerated maps. Existing qualified CIFs and submitted jobs are not edited
in place.

The installed RFD3 documents explicit source-anchored contigs with generated
segments and omitted source residues, supporting loop-length changes. This
is distinct from the currently qualified fixed-length partial-diffusion path:
the installed parser rejects `length` with partial diffusion. A loop-indel
input/export/graph qualification is required before production use. See the
[official contig specification](https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/input.md#contig-strings)
and installed `models/rfd3/docs/input.md` lines191–211 and
`models/rfd3/src/rfd3/inference/input_parsing.py` lines279–282,628–652.

For a first proposed arm, compare loop lengths L−1, L and L+1, with at most
two backbones per length and two sequences per accepted backbone (six
backbones/twelve raw sequences maximum). These are proposed additional
budgets, not an expansion of the frozen 24-sample pilot or a launch request.
Choose the exact loop length per trial and apply identical edits to all four
corresponding copies; do not independently sample their lengths. Preserve
four-copy sequence symmetry and actual assembly geometry. Qualify closure,
stereochemistry, clashes, TS-motif retention, sequence/graph mapping and the
complete Amber/predictor roundtrip before extending the arm.

Maintain three explicit identities: WT reference residue, sequential design
chain/residue, and sequential Amber residue/atom index. Give inserted residues
stable metadata identities and mark deletions explicitly; do not use PDB
insertion codes in simulation inputs. The WT reference label Thr233 remains
an annotation if its numerical design position changes. Rebuild precursor/
mature-chain cleavage maps, protection masks, symmetry groups, terminal charge
edits, predictor sequences, topologies, QM masks and analysis selections for
each length. Hardcoded 232/96,328,1328,229/99 and Amber929 assumptions are
valid only for the current fixed-length campaign. Standard inserted amino
acids use the existing protein force field, and the unchanged PA6 release
remains reusable; each altered protein requires a fresh topology. No indel
calculation has been launched under this proposal.

### Broaden sequence freedom and test backbone remodeling separately

Allow **all precursor positions except an explicit protected set** to be
considered for the broad mask; do not reduce it again to five nearest sites
or require every mutable residue to be solvent exposed. Freeze known catalytic
and maturation identities, the Asp272–Gly273–Asp274 motif in the first test,
disulfides and required inter-repeat interface contacts. Audit conserved
structural interactions before changing them. Publish the actual canonical
mutable/fixed lists and counts before generation; no new exact count is claimed
until that structure-based census is complete. Treat previously protected
positions as conservative research choices where essentiality is unproven.
Use a declared, common native-identity preference at distal sites if needed
to control gratuitous mutation burden; broad permission need not force every
permitted site to change. Do not apply the old global RSA filter to the pocket,
second shell or buried packing positions needed for coordinated redesign.

Use the **same broad sequence mask** on original TS parents and remodeled
parents. This separates the benefit of sequence freedom from that of backbone
changes. Fixed-backbone broad design is a useful matched control, not a reason
to defer the backbone test. The objective is a well-organized reactive ensemble
that is accessible from E–S. Excessive flexibility can reduce preorganization
and increase the entropy cost of reaching that ensemble; structural diversity
alone is not evidence of improved dynamics or a lower barrier.

For backbone remodeling, initially retain WT mature lengths and connectivity:
four alpha chains of 232 residues and four beta chains of 96 residues, with
the authentic processed termini. Keep the TS reactive heavy-atom motif and
necessary catalytic anchors fixed while allowing selected pocket/entrance
loops and their structurally coupled flanks to move. Include the whole
tetramer as context and preserve inter-repeat interfaces and maturation
architecture. Sequence freedom can extend well beyond the segments whose
backbone coordinates are remodeled. Do not introduce indels or redesign the
whole fold in this first qualification.

### RFD3 qualification and tetramer symmetry

The installed Foundry code supports partial diffusion with atom-specific
coordinate masks; coordinate and sequence masks are independent. `partial_t`
is a noise amplitude, not a timestep count or a promised final RMSD. Its
default fixed-atom mask is true, so an explicit frozen complement is essential
to get actual remodeling. See the
[official input documentation](https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/input.md).
Local inspection identified checkout `c97ad246877e2c7cb3b50e60aab7d76d02e8d393`
and installed package `0.1.11.dev34+gcee116dc5`; pin both. A CPU diagnostic of
the imported sampler confirmed that partial noise reaches mobile coordinates
and leaves fixed ones unchanged. Subsequent full GPU runs and source-validated
exports qualified the two 5 Å TS-parent remodels on 27 November; see the dated
native acceptance manifests. This does not qualify length-changing inputs.

Exact four-copy **sequence symmetry is mandatory**. Identical sequences need
not have identical conformations, particularly with one occupied active site.
The directly inspected nonsymmetric partial-diffusion path retains the eight
physical chains and can remodel the occupied repeat within a fixed tetramer
context; subsequent LigandMPNN ties alpha copies separately from beta copies.
Do not fuse the mature chains or independently generate/average four loops.
Strict geometric C/D symmetry is an optional qualified route: use measured
assembly operators, not an assumption that tetramer means C4. The prior
scaffold audit was D2-like. A single TS ligand must not become four invented
occupied sites to satisfy a generator's assumptions. Check fixed-motif
retention after inverse folding, relaxation and independent prediction.

The PI recalled a prior RFD3 tetramer setup in `redesigner-agent`. The inspected
checkout contains demonstrated LigandMPNN `--homo_oligomer 1` usage on four
merged 328-residue input chains and a folding split of 232/96. Targeted
source/configuration and ignored-file searches found no actual RFD3 input,
command or execution artifact there. This supports reuse of earlier sequence
symmetry ideas, not a claim that RFD3 partial diffusion on today's eight mature
chains is already validated. Preserve these different representations/maps.

Before GPU use, audit custom P6N/P6I/P6C import: the one L4 molecule spans four
residues on chain I. The explicit multiresidue-ligand option renumbers residues;
all 33 ligand heavy atoms, inter-residue bonds and their source correspondence
must survive. Preserve the TS geometry as atom coordinates with an explicit
reactive-state annotation, not an invented stable covalent FF model. Do not
add terminal capture/burial conditioning. Retain proton inventories separately.

### Finite new exploration and prediction budget

**Execution refinement, 27 November:** the first broad pilot is capped at
24 raw sequences allocated equally across four templates: original TS
parents01/03 and their independently accepted5 Å remodels, six sequences each.
Both parent adapters passed their native/independent audits. The source-hashed
queue generated all 24 unique designs with the same 229-position mutable mask,
99 protections and four-copy ties, and completed native TS/E–S cross-scoring.
This replaces the provisional three-per-four-remodel allocation below for
this first batch; additional RFD3 expansion is reserved, not automatically
queued. Deduplicate and cross-score complete sequences against all four TS
templates and the matched E–S reference. Do not equate generation with
experimental nomination or launch variant predictions on a failed WT baseline.

The original query-only WT controls and the two subsequent homologous-MSA
controls failed to recover Tier A/B sites. Therefore the prediction-calibration
branch tests one protein-only crystal-template control per failing predictor,
with PA6 still unguided. Keep that calibration separate from broad sequence
generation. Original control slots and explicitly authorized diagnostic
follow-ups are separately counted; do not rerun the old eight-job controller.

Every published predictor result includes biological numbering (A–D alpha
1–232; E–H beta233–328), mature-local numbering (beta1–96), and an Amber
heavy-input export using the established modular converter. Four-ligand
controls have 1,328 sequential dry residues; the one-ligand TS parents have 1,316.
After pmemd/sander, use the matching topology, restart and conversion map to
restore canonical identity before RFD3/MPNN. Reapply neither a guessed offset
nor a graph inferred from close contacts in a reacting geometry. Retain the
native model and exact bidirectional atom/residue maps. Heavy-input exports
do not establish hydrogen placement, protonation or MD topology acceptance.

1. **Qualify two RFD3 outputs**, using one accepted TS parent at two modest
   noise settings, provisionally `partial_t=2` and `5` Angstrom. Require actual
   backbone movement, fixed catalytic motif retention, intact PA6 atom/graph
   mapping, correct chain breaks/lengths, loop closure and acceptable clashes.
   Noise amplitude is not an output-displacement tolerance. Failed import or
   retention stops expansion; do not substitute a different ligand silently.
2. If qualified, make **at most eight backbone outputs total**, including
   those two: two source TS parents × two noise settings × two random seeds.
   Select at most one valid diverse backbone per parent/noise combination
   (four total). Keep all rejected attempts and reasons.
3. Sample **24 new raw sequences maximum**: twelve on the two original
   backbones (six each), and twelve on up to four remodeled backbones (three
   each). Use the same broad mask, protections and sampling policy in both
   arms. Fewer valid backbones reduce the budget; no quota-filling retries.
   Keep the existing narrow A187S sequence as a baseline without regenerating
   the old pilot. Deduplicate whole sequences and cross-score all parents,
   matched E–S and available acyl references; retain worst-parent behavior,
   all-copy compatibility, mutation burden and diversity.
4. For independent E–S prediction, use **no substrate contact restraints**
   and do not supply the designed TS backbone as the prediction template.
   AF3 is a reasonable candidate, but first compare WT and one preselected
   representative broad design, two matched seeds each, with AF3 and Boltz:
   **eight predictions**. Use the same declared MSA/native-template policy
   across sequences; keep designed-backbone templates out of this independent
   check. Choose the primary predictor using WT assembly/chemistry/pose quality,
   runtime and memory, not whichever favors a particular design. Current AF3
   timing evidence covers a smaller one-site system, not the full tetramer;
   full-system hardware compatibility/cost must be measured. Boltz's validated
   four-L4 tetramer logs took about 7.6 minutes/model. Do not extrapolate the
   small AF3 timing to the larger assembly; see its
   [official performance guidance](https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md).
5. Use the selected unguided predictor for **up to six designs + WT, two seeds
   each** (14 predictions, reusing the four matching pilot predictions).
   Include representatives of broad-original and broad-remodeled arms plus
   the narrow baseline where useful. Then add three seeds each for top3+WT
   (12 more). Including the four pilot predictions from the other predictor,
   the entire new prediction ceiling is **30**, not 30 per arm. All generation
   failures remain in the denominator. Relax with the same unrestrained
   noncovalent protocol and no TS bonds/FF assumptions.
6. Retain the existing downstream ceiling: **two designs + matched WT**, two
   100 ns MD replicas each; at most six QM/MM optimizations and three matched
   coordinate-driving paths; US/MBAR only for the strongest design and WT.
   Use the declared common occupancy/protonation/Hamiltonian protocol. Broader
   sequence changes increase the importance of maturation/fold checks; they
   do not justify more expensive campaigns automatically. SolubleMPNN follows
   selection of a catalytic parent and remains separately controlled.

The immediate deliverable before November 30 is a functioning, qualified
backbone/sequence/prediction funnel and a defensible shortlist. Do not promise
converged new mutant PMFs or experimental catalytic enhancement by that date.
Current relay-PMF convergence remains the scientific CPU priority. Only Violet
and Azalea are in scope for this new work. No new dependency is approved by
this planning update; request approval if qualification requires one.

### Register-independent analysis still requires exact chemical identity

Enumerate **every internal amide** from the actual PA6 graph at every occupied
site. For each candidate bond, use its own carbonyl C, O and amide N together
for OG1 attack, oxyanion contacts, Bürgi–Dunitz angle and N-terminal relay
geometry. Exclude terminal carboxyl groups from the amide list. Accept any
coherent register/direction that supports the intended chemistry; do not
require the legacy C54/O49/N59 labels or enforce one phi sign before evaluating
the local geometry. Retain signed Herschlag phi, exact atom map, selected amide,
orientation and all other tested registers in the report. A single ligand at
two possible sites is not two independent recoveries; report denominators.

During MD, report both any-amide catalytic occupancy and the persistence of
each particular register. Do not silently switch which C/O/N defines a single
reaction coordinate from frame to frame. Before QM/MM/RCD/US, freeze the
selected amide identity and rebuild/audit its QM boundaries, charge and product
inventory. Different registers can change terminal-charge proximity and
which atoms belong in QM. Candidate PA6 directionality/register diversity
remains a result to report, not an undefined atom-selection shortcut.

### Assay constraints change what success means

The PI clarified that the current mass-spec assay quantifies only products
with standards, currently **PA6 L2, L3 and L4**, and that testing will probably
remain in crude lysates. Larger oligomers are not quantitatively covered and
may remain undetected by the current targeted method until further hydrolysis.
Do not infer absence of cleavage from low L2–L4 signal, assume an assay for
monomer without confirming its standard, or interpret summed calibrated
L2–L4 yield as total polymer conversion. Earlier papers/figures can have
different analyte panels; preserve assay-version provenance.

The practical experimental objective is improved **calibrated L2–L4 product
yield under matched lysate conditions**, reported product by product, with
expression/maturation/stability contributions acknowledged. It is distinct
from intrinsic kcat/Km and from total amide-bond cleavage. Favor substrate
access and product release as well as chemical compatibility; cleavage only
to larger soluble fragments may be useful but requires a complementary assay.
Report uncertain assay-negative variants as having no detected calibrated
short-product improvement, not as chemically inactive enzymes. K59E remains
an unresolved noisy measurement, not a clean label.

**PI confirmation, 26 November:** purification of WT plus a few promising
designs should be feasible. The planned validation panel is **WT plus two or
three selected designs**, using the existing validated construct/purification
approach. The PI also favors developing a primary-amine-release assay even
if it requires changing buffer. This is an agreed experimental direction,
not a claim that the assay is validated or that new experiments have started.
Compare mature fraction, soluble amount, assembly/stability and
time-dependent product formation. Purification permits enzyme-normalized
comparisons and cleaner orthogonal assays; it does not automatically establish
active-site concentration or Michaelis–Menten kinetics on a heterogeneous solid.
If purification is unavailable, retain lysate screening with matched culture/
lysis/loading, biological replicates, WT on each plate, empty-vector and time-zero
controls, standard-spike recovery and soluble processed-target quantification
where feasible. Total lysate protein alone does not normalize target expression.

Cell-free synthesis is a future expression/screening option, not a remedy for
the current analytical uncertainty. A useful primary precedent is
[Landwehr et al., Nature Communications 2025](https://doi.org/10.1038/s41467-024-55399-0),
which used cell-free engineering and subsequently purified the best enzyme for
comparison. Nyl12 would need a WT proof of soluble folding, tetramer assembly,
autocatalytic maturation and activity before a cell-free library. Cell-free
reagents/lysates add amino acids and other background; purification/cleanup
can still be required. Do not assume that easier synthesis means equal active
enzyme concentration or correct processing.

Useful complementary assays, in order of relevance to the present limitations:

- **Retain calibrated targeted LC–MS:** the nylon-specific
  [Bell et al., Nature Communications 2024 study](https://doi.org/10.1038/s41467-024-45523-5)
  and [authors' UHPLC–MS/MS protocol](https://doi.org/10.17504/protocols.io.kxygx331dg8j/v1)
  are practical starting points. A targeted method still requires analyte
  standards; changing instruments alone does not remove that limit.
- **Broader LC–HRMS/MS profiling of the soluble fraction:** seek higher oligomer
  series and use exact mass/MS/MS for assignments; without authentic standards,
  label their abundance as relative/semiquantitative, not absolute conversion.
  First check whether existing raw data already contain full-scan spectra;
  a targeted-only acquisition cannot be retrospectively expanded into them.
- **Soluble primary-amine release:** TNBS has direct NylC precedent in
  [Negoro et al., JBC 2012](https://doi.org/10.1074/jbc.M111.321992) and
  [Shiraishi et al., PLOS ONE 2025](https://doi.org/10.1371/journal.pone.0318641).
  OPA/fluorescamine are alternatives requiring their own validation. Current
  200 mM Tris contains a primary amine and lysates contain many reactive amino
  groups; do not apply these methods directly to that matrix. Validate an
  amine-free assay/cleanup protocol, blanks, recovery and detection range.
  Amine equivalents in the soluble fraction do not by themselves count every
  bond broken throughout the original insoluble polymer.
- **SEC/GPC on recovered polymer:** a molecular-weight distribution shift
  complements short-product release and can reveal chain scission. It needs
  appropriate polyamide dissolution, calibration, sampling and enough signal.
  Gravimetric loss or surface imaging alone is weak evidence at low conversion.
- **Analytical post-hydrolysis, proposed rather than validated:** isolate the
  soluble reaction fraction after removing original polymer, then convert its
  oligomers exhaustively to calibrated 6-aminohexanoate for LC–MS quantification.
  This could estimate soluble PA6 repeat-unit recovery with one standard and
  detect larger-fragment release. Establish quantitative digestion and recovery;
  do not assume NylB alone converts every chain length, or digest the original
  polymer during the analytical step. Include processed lysate blanks and
  authentic oligomer spike controls before interpreting totals. Crucially,
  [6-aminohexanoic acid](https://pubchem.ncbi.nlm.nih.gov/compound/564) and
  [isoleucine](https://pubchem.ncbi.nlm.nih.gov/compound/6306) (also leucine)
  share formula C6H13NO2. Lysate protein hydrolysis would release the latter
  abundantly; exact precursor mass alone cannot distinguish them. Demonstrate
  chromatographic/MS/MS discrimination with authentic standards and matrix
  controls. This makes post-hydrolysis a method-development option, not the
  immediate replacement for the calibrated short-oligomer assay.

Defined oligomers/fluorogenic analogues can diagnose bond chemistry but do not
replace the actual polymer assay. Keep polymer lot, particle size, pretreatment,
temperature and mixing matched; apparent gains can arise from accessibility.
No new wet-lab assay is claimed validated by this plan.

### Primary-amine assay development — agreed next experimental step

Start with **purified WT**, then apply the accepted method to the small finalist
panel. TNBS is the first literature-supported candidate; OPA/fluorescamine remain
alternatives if their sensitivity and background are better after validation.
Retain calibrated L2/L3/L4 MS as the product-specific channel and, where feasible,
split aliquots from the same reaction for both readouts. Do not replace one
unqualified endpoint with another.

1. **Select a compatible reaction buffer.** Evaluate an amine-free option such
   as phosphate, checking WT activity and stability at the intended assay
   temperature. Distinguish the hydrolysis buffer from the derivatization
   reagent/buffer. Record pH and its measurement temperature; a nominal pH
   match between Tris and phosphate does not ensure the same pH at 75 °C.
   Run paired MS measurements in the original and candidate reaction buffers
   to distinguish buffer effects on hydrolysis from changes in detection.
2. **Control enzyme and reagent backgrounds.** Include reagent/buffer blanks,
   substrate-only/no-enzyme samples, time-zero samples, and enzyme-only
   time courses under the same heating conditions. Proteins themselves carry
   reactive amines, and continuing enzyme maturation/proteolysis can produce
   new termini even without PA6. Check the mature fraction before/after
   incubation; purification alone does not remove this confound. Use a
   validated protein-removal procedure if needed before derivatization.
3. **Validate calibration and sample recovery.** Establish detection and
   quantification limits, usable range and precision at the low product
   concentrations actually expected. Calibrate amino equivalents with
   6-aminohexanoate and check available L2/L3/L4 standards rather than assume
   identical responses. Spike before and after any cleanup to distinguish
   recovery loss from matrix effects. Use validated non-nylon consumables;
   report filtration/size cutoffs and any loss of larger oligomers.
4. **Measure time-dependent activity under matched conditions.** Use the same
   polymer lot, particle preparation, loading, temperature, mixing and
   quantified mature-enzyme loading for WT/designs. Report independent
   expression/purification batches separately from technical replicate wells.
   Include active WT on each plate and seek a time interval supporting an
   apparent initial-rate comparison. A prolonged endpoint also reflects
   thermal inactivation and accessibility; mature fraction does not establish
   the concentration of catalytically active sites.
   Retain the existing lysate screen for throughput, with expression/processing
   controls; do not pool it with purified-enzyme measurements as one assay.
5. **Keep the two inferences distinct.** Report background-corrected soluble
   primary-amine equivalents and calibrated individual L2/L3/L4 yields.
   More amine release with little change in calibrated short products may
   support release of larger fragments, but requires chemical corroboration.
   Cleanup or solubility limits can exclude products. Neither readout alone
   establishes total polymer conversion, all bond-cleavage events or kcat/Km.

Advance from WT assay qualification to finalist comparisons only when blanks
are controlled, expected signals are quantifiable, spikes recover acceptably
and paired MS confirms activity in the selected buffer. Freeze numerical
acceptance thresholds with the experimental collaborators before testing the
designs; no laboratory-specific values are invented here. Cell-free synthesis
remains optional later, and broader MS/SEC/post-hydrolysis options remain
separately qualified follow-ups rather than new requirements for this pilot.

## Main and possible secondary bottlenecks — clarification, 2026-11-27

The current main projected maximum remains s ≈ 0.69. The s ≈ 0.4 shoulder
and s ≈ 0.54 hump are not yet additional established transition states.
The shared TI1-like reference maps to s = 0.537776; the subsequent shallow
free-energy trough is strongly coordinate sensitive (SI S9). Retain a small
chemically distinct early addition/proton-transfer template set from the
s ≈ 0.4–0.6 search band alongside the main barrier templates, if the structural
audit supports distinct reactive configurations. This is a bounded secondary
arm/control within the agreed total budget, not authorization to multiply
the sequence or prediction campaign. Do not uniformly sample every s value.

Whole candidate sequences must be cross-evaluated across both reactive sets
and E–S/acylenzyme compatibility controls. LigandMPNN compatibility is not a
free-energy difference. Prefer balanced candidates without a severe failure
at one relevant state; do not simply reward tighter binding to TI1, which
could make its subsequent escape harder. Reducing the main barrier can expose
an earlier or later bottleneck; the final choice requires matched chemical
barrier comparisons for the few shortlisted designs. Existing backbone and
sampling gates remain in force; no new design run is launched by this note.

## Bounded relay-ensemble redesign campaign — historical first pilot, 2026-11-26

The PI authorized a small exploratory LigandMPNN ensemble based on the actual
relay-US frames near **s = 0.7**, the current high-free-energy region. Execution
and immutable inputs belong in
`design/nyl12_relay_ts_ensemble_20261126/`. This dated plan supersedes the older
24-sample split that required an available optimized saddle. It authorizes
**24 raw samples: eight per each of three actual parent frames**; it does not
authorize an unbounded search until a preferred number of unique designs appears.
The running relay US and its targeted initialization/endpoint checks take CPU
priority. The initially idle Violet GPU 0 can handle the bounded sequence pilot.

### Hypothesis, parents and sequence constraints

The hypothesis is that new local packing and hydrogen-bond interactions can
favor the relay's reactive heavy-atom geometry relative to the corresponding
E–S ensemble while retaining maturation, assembly, substrate entry and product
release. Tighter E–S binding alone need not lower the activation free energy.

Select three **observed** medoids from retained s in [0.65, 0.75], using reactive
heavy-atom and pocket geometry, conditional regional MBAR weights and distinct
source windows where possible. Preserve the exact parent coordinates, window,
native frame/time, bias, proton inventory, atom map and input hashes. Do not
average coordinates or minimize every parent into one conformation. Current
regional weight ESS is about 35.7 and the largest frame contributes about 14.8%
of regional weight; the very small globally normalized probability of visiting
the barrier region is not itself a reason to reject its conditional ensemble.
Three templates remain a bounded structural sample, not a convergence claim or
three independent molecular-dynamics replicas. Keep source coordinates with
hydrogens and a separately mapped heavy-atom LigandMPNN input.

Each source contains the complete mature tetramer and **one physical substrate**.
Do not invent three additional TS ligands. Canonical protein chains A–D contain
alpha residues 1–232 and E–H contain beta residues 233–328; paired mature repeats
are AE/BF/CG/DH. Design one precursor sequence with exact four-copy alpha and
four-copy beta sequence ties. For the one occupied active site, sampling uses
symmetry weights **1,0,0,0** at each designed position, so the other three apo
sites do not outvote the TS context. Audit the whole assembly afterward; these
weights are a sampling choice, not a joint multistate energetic objective.

Start from the previously tested focused pocket mask. Intersect it with the
fresh templates' local reactive-core contacts, and remove disulfides, protected
catalytic/maturation positions and inter-repeat interface contacts. Do not
exclude every alpha–beta contact within a mature repeat: that would remove
much of the catalytic pocket itself. Do not add positions whose only contact
is with an artificial short-oligomer terminal group. Freeze Thr233, Asp272,
Asp274, Phe113, Asn186, Lys156 and the established maturation protections;
the exact audited mask and per-position exclusion reasons are authoritative
in the campaign manifest. Use fixed-side-chain context and verify that the
parsed catalytic atoms are present. All mutation reports use canonical
precursor numbering, with mature-chain maps retained separately.

Fresh preflight narrowed the first pilot to **V111, I112, D114, A187 and V188**
in all four alpha copies. Beta chains are fixed. The initially eligible G273
is conservatively frozen to retain the Asp272–Gly273–Asp274 loop in this first
test; that is a scope decision, not a mutational finding. Charge changes at
D114 must be reflected in later topology/counterion/QM-region audits.

### A finite funnel for complete sequences

| Stage | Maximum new work | Advancement criterion |
|---|---|---|
| LigandMPNN | 3 parents × 8 = **24 raw samples** | Exact symmetry/protection checks; deduplicate complete precursor sequences and retain every sampling lineage, including WT duplicates |
| Cross-template assessment | All unique samples plus WT on all 3 parents | Mean and worst-template local compatibility, all-copy context costs, mutation burden and sequence diversity; retain complete sequences, never combine position-wise votes |
| Initial Boltz-2/relaxation | **Up to 8 designs + WT**, 2 matched seeds each = **18 models** | Assembly/fold integrity; chemically mapped, unrestrained post-relaxation geometry and ensemble consistency; sequence likelihood is supporting evidence only |
| Focused Boltz confirmation | **Top 3 designs + WT**, 3 additional seeds each = **12 models** | Two additional guided seeds and one guidance-free sensitivity seed per sequence; total ceiling **30 predictions**, with guided/unguided arms reported separately |
| Classical MD | **Up to 2 designs + WT**, 2 independent replicas × 100 ns = **6 trajectories** | Persistent binding/relay and oxyanion geometry, stability across replicas, no maturation/interface penalty; use a matched 20 ns checkpoint before committing the remaining time |
| QM/MM optimization | Same up-to-3 sequences, up to 2 matched snapshots each = **6 optimizations** | Correct proton/charge inventory; comparable unrestrained local relaxation and chemically viable reaction geometry |
| Coordinate driving | **Up to 2 designs + matched WT**, one identical-protocol path each = **3 paths** | Candidate-relative E–S-to-barrier changes, chemistry/connectivity and comparable sampling of environment; a PES improvement is provisional |
| US/MBAR | **Only the strongest design + matched WT** | Sufficient evidence from earlier stages, a defined common protocol and available resources; converged uncertainty/overlap/hysteresis checks, not a fixed trajectory length alone |

These are ceilings, not quotas. Zero candidates advancing is an acceptable
outcome. Review the retained Pareto shortlist after each stage; do not choose
eight nearly identical sequences to fill slots. Cross-score every complete
sequence on all parents; record its worst-parent behavior and mutations instead
of presenting a sequence designed on one frame as an ensemble optimum. An
E–S and acyl-state compatibility comparison is required before final chemical
nomination, but an unavailable exactly mapped state must be marked pending,
not silently replaced. LigandMPNN structural compatibility across states is
not a TS-minus-E–S free-energy calculation. Different mutant absolute energies
must never be compared as barriers.

Only the sequence pilot is launched by the present instruction. Structure/MD
and chemistry stages remain bounded follow-ups selected from actual results.
Use the existing WT, K59D and V111C evidence to challenge whether a structural
gate is useful; earlier known-active controls sometimes failed those gates.
K59E's discordant lysate/MS measurements remain unresolved and cannot serve
as a clean negative or quantitative activity label. Native-Thr similarity is
diagnostic, not a hard rejection filter. Four sites in one model/trajectory
are nested observations, not four independent seeds/replicas. Current noisy
lysate product assays do not provide intrinsic turnover constants or reliable
quantitative training labels.

### Boltz inputs, contact guidance and structural assessment

Use the mature tetramer and **four noncovalent PA6 L4 molecules**, one per site,
for the model-comparison stage, consistent with the PI's existing instruction.
This differs explicitly from the single occupied TS source. Boltz predicts
E–S-compatible structures; it is not required to retain the reactive substrate
bond geometry. Use the validated PA6 graph/atom map at each ligand instance,
the same crystal template/MSA policy and matched random seeds for WT/designs.

Retain minimal broad contact guidance as a **model-building aid**: the previous
two contacts per site can be reused only after their atom/residue mappings and
protected anchor identities are audited. Do not set the exact distances later
being used as success criteria. Remove all generation contacts for the validated
Amber/OpenMM relaxation; report their initial satisfaction separately. Include
the one guidance-free seed per finalist and WT, without pooling that small
sensitivity arm with the guided arm as if all models had the same protocol.
Whether a guided pose survives unconstrained minimization and MD matters more
than whether the predictor satisfied a supplied contact. Do not reject a known
active control solely because an unguided predictor failed to place its ligand.

Evaluate all four sites with a verified graph and canonical numbering: signed
Herschlag phi/register, OG1–carbonyl-C distance, Bürgi–Dunitz angle, both oxyanion
contacts, N-terminal-N–leaving-N separation and potential relay contacts. Report
Thr chi1/psi, crystal-aligned N/OG1 displacement and CB–OG1–carbonyl-C angle as
diagnostics. Add local clashes, unsatisfied polar contacts, fold/assembly
confidence, inter-repeat interfaces, terminal-only interactions and access to
the acyl/deacylation geometry. Show distributions and denominators across seeds
and sites rather than a single best-scoring picture. E–S Tier A/B criteria do
not apply unchanged to tetrahedral or TS geometries.

Amber residue 929 is an identity in the declared topology ordering, not a
preferred-site rule. Rank eligible site–ligand–amide combinations across all
four sites before selecting the mechanistic system. Record the physical
alpha/beta chains, ligand chain/residues and coherent carbonyl-C/O/leaving-N
mapping. Derive downstream QM masks, driving coordinates, restraints and
viewer selections from the selected conversion map; the best pose need not
involve residue 929 or the historical substrate C54/N59 pair. Contact-restraint
policy does not alter this requirement. Preserve every site's measurements
and the one-to-one assignment as well as the final selection rationale.

For the expensive mechanistic MD/QM/MM stage, standardize to **one physical
PA6 L4 molecule per tetramer** in every design and WT, matching the present
mechanism campaign. The four-ligand Boltz screen above tests pose compatibility
at all equivalent sites; it does not establish the physical occupancy. Retain
one mapped site before identical full-system preparation, resolvation/ion
accounting, minimization and equilibration; do not silently delete ligands
from an existing simulation or compare different occupancies. Fix parameters,
mature termini, protonation, ionic conditions, 300 K and preparation across
sequences. Select a representative valid model per sequence, then use two
independent velocity seeds; do not launch many MD starts and retain only the
ones that bind. Reuse old WT sampling only after a protocol audit establishes
the same occupancy and Hamiltonian; otherwise prepare a fresh matched WT.
In QM/MM, retain the remaining physical system as MM context. Use the validated QM subsystem,
charge accounting and 12 Å movable-region protocol consistently; preserve the
declared frozen exterior. No new classical parameters for a TS/covalent state
are introduced by this pilot. Evaluate reaction barriers relative to each
sequence's own matched E–S state. Replicate variability and protocol limitations
must accompany any apparent improvement.

### Timing and stopping

Archived matched logs give approximately 7.6 min/Boltz model, so the complete
30-model ceiling is about **4 GPU-hours** of prediction, plus preparation and
relaxation. A 100 ns MD replica previously took about **10 hours** of production
at 224–247 ns/day; six replicas cost about **60 GPU-hours**, roughly 30 hours
on two GPUs before equilibration/queue overhead. These are measured historical
estimates, not live performance guarantees. Use one idle Violet GPU initially
and preserve scientific CPU jobs. Do not recruit shared hosts without renewed
authorization. A sequence shortlist and initial structural screen are realistic
before November 30; rigorous new mutant PMFs are not promised by that date.

Stop and diagnose the preparation if mapping, protected residues or symmetry
fail. If a structural gate rejects WT/known-active controls, investigate its
calibration instead of making it stricter or silently changing thresholds.
If the sequence pool is repetitive or collapses on all parents, report it
before expanding masks. Advance at most two diverse designs to MD; no automatic
second large prediction/MD campaign. Reserve expensive PMFs for evidence that
justifies the next comparison. Keep independent native jobs running throughout.

### Solubility and product selectivity are separate design variables

Defer SolubleMPNN to a second, matched surface-design layer on a selected
catalytic parent. Protect the reactive network, maturation region, conserved
functional/structural sites and inter-repeat interfaces; solvent exposure and
distance alone do not establish dispensability. Keep four-copy sequence ties.
A later bounded parent-versus-surface-variant comparison can assess fold and
hydrophobic patches, but expression/solubility and correct maturation still
need experiments. Do not confound the first mechanistic pilot with extensive
simultaneous surface redesign.

The first objective is **improved PA6 cleavage**, with measured product profiles,
not enforced carboxyl-end/dimer-selective cleavage. Nyl12 predominantly produces
linear trimers in the Drufva study; L4→L2 is our tractable modeling hypothesis,
not an established isolated-substrate reaction. A central cut of L4 cannot
distinguish endo activity from release of a dimer at either end. The corresponding
long-polymer chain ends, directionality and attacked register remain unresolved.

Do not deliberately reward short-L4 terminal ion pairing in the first mask.
Carboxyl-end recognition could later be a useful specificity objective, but
requires a separate end-binding/register hypothesis, longer-chain or capped
controls, product-exit checks and product-resolved assays. An end-only enzyme
may have fewer accessible attack sites on a long polymer; stronger end/product
binding could also cause product inhibition. These are testable concerns, not
established Nyl12 limitations. Do not begin another parameterization campaign
for capped polymers merely to enlarge this pilot. Distinguish the mature
RESP/GAFF2 zwitterionic-L4 release from provisional AM1-BCC capped controls.

## Gates followed by a diverse Pareto shortlist

### Reference states and design readiness — updated 2026-11-25

The current priority is a statistically stable, chemically interpretable WT
US/MBAR profile. Barrier convergence, coordinate adequacy and method limitations
must be assessed independently of whether the number looks favorable. This
supersedes the older requirement to obtain a stable TI1 before TS-guided design.
Prepare matched noncovalent E–S, each relevant acylation bottleneck, and
acylenzyme references with explicit atom/proton inventories. Include a relay
or tetrahedral intermediate only if supported; a discrete TI1 is optional.
LigandMPNN can use exploratory structural templates before full mechanism
certification, but such designs must retain their provisional provenance.
Claims of a true PES TS require force convergence, one chemically relevant
negative mode in the declared optimization space and downhill connectivity;
an umbrella frame or a partial-QM Hessian alone does not establish that claim.

Use the existing apo/WT assembly as a structural control. An acyl-enzyme and a
product-bound/deacylation-water compatibility model should be included before
final candidate selection; full free-energy pathways for all later steps are
not prerequisites for exploratory LigandMPNN sampling. Preserve symmetry,
residue/atom maps, maturation residues and substrate register across states.
Distinguish early template-based sampling from validated TS-guided design.

### Clarification: exploratory templates before final PMF convergence — 2026-11-26

Neither a globally converged PMF nor an optimized first-order PES saddle is a
LigandMPNN input requirement. A small, explicitly exploratory pilot may use
three to five chemically consistent, distinct barrier-region frames while the
full PMF is still being qualified. Select across plausible bottleneck locations
and conformations when the peak location remains uncertain; preserve the
actual parent coordinates, atom/proton inventory, bias and source provenance.
This does not relax the convergence requirements for quantitative barriers or
justify claiming an experimentally improved catalyst from sequence likelihood.
The priority remains finishing the WT relay analysis before a large campaign.

A verified PES saddle is a useful chemical anchor and single-template control,
but need not replace finite-temperature snapshots. Do not optimize every MD
parent to a saddle or unrestrained minimum before ensemble design: that can
remove the diversity being tested. A committor-qualified transition-state
ensemble is distinct from frames near a projected PMF maximum; unbiased
shooting toward defined adjacent basins tests that dynamical claim
([McGibbon et al., JCP 2017](https://doi.org/10.1063/1.4974306)).

The installed LigandMPNN parser masks ligand H atoms (`Y_t != 1`), and its
atomic context encodes elements and geometry rather than electronic charge
redistribution or activation energies. Fixed catalytic side chains require
`--ligand_mpnn_use_side_chain_context 1` with the correct fixed-residue mapping;
audit the parsed context so Thr OG1 and protected catalytic partners are visible.
Keep protonation/bond-order metadata for subsequent QM/MM checks even when
those attributes are not model features. This follows the
[primary LigandMPNN model description](https://doi.org/10.1038/s41592-025-02626-1)
and the installed `data_utils.py`/README, not a claim that LigandMPNN scores
proton transfer. Evaluate complete proposed sequences across E–S, candidate
reactive templates and acylenzyme controls; preserve assembly sequence ties
and maturation/catalytic protections. No new design run is authorized by this
clarification.

### PES saddles and transition-region ensembles — proposed 2026-11-25

Use a hybrid reference set rather than choose one instantaneous frame or one
optimized geometry as universally best. This is a proposed comparison, not a
completed design campaign or a claim that multiple snapshots guarantee activity.
An unfinished PES-saddle certification must not indefinitely block a credible,
qualified ensemble-based pilot: defer the unavailable single-saddle control and
retain the narrower transition-region design claim until that control is ready.

| Structural input | Purpose | Qualification and limit |
|---|---|---|
| Verified PES first-order saddle | Reference reactive geometry and a single-template control | Method/model-specific stationary point; not itself the finite-temperature transition-state ensemble |
| Three to five actual representative frames per relevant PMF bottleneck | Preserve sampled protein/substrate/solvent arrangements for an ensemble design arm | Select from equilibrated, decorrelated, reweighted data; call these transition-region templates until dynamical qualification |
| Matched E–S, acylenzyme and available apo/open/product-compatible models | Check selectivity, turnover compatibility and access | Strong E–S or intermediate stabilization alone is not evidence of a lower activation free energy |

1. **Establish the WT profile first.** Apply the eight sampling gates in
   `qmmm_umbrella_sampling_and_pmf_protocol_20261109.md`. Examine both chart
   families, proton coordination/ownership and orthogonal relaxation. A PMF
   maximum along a selected coordinate need not be the dynamical dividing
   surface, and a stable profile can still describe an incorrect or unfavorable
   modeled mechanism. Do not tune the protocol to obtain a target barrier.
2. **Select a small ensemble.** Identify every kinetically relevant barrier,
   not merely the highest single sampled energy. Cluster decorrelated frames
   near each barrier using mapped reactive geometry, backbone/side-chain
   arrangements and hydrogen-bond/solvent networks; retain three to five actual
   medoids initially, fewer if no distinct clusters are supported. Record source
   window/block/frame, sampling weight, chemical inventory and cluster coverage.
   Raw umbrella frame counts are biased and correlated. Reweighted cluster
   populations describe the selected coordinate band, not automatically reactive
   flux or committor-defined TSE weights. Do not average atomic coordinates.
3. **Qualify chemistry without quenching away diversity.** Use the existing
   verified PES saddle, or a bounded targeted search from selected representatives,
   as an anchor. Retain unminimized ensemble parents separately: TS optimization
   can collapse distinct thermal conformations to the same stationary geometry.
   Where transition-state-ensemble claims or coordinate adequacy matter, test
   representative configurations with bias-free QM/MM shooting at the same
   physical temperature/model. Draw independent thermal momenta and count
   commitment to explicitly defined adjacent basins. A committor near one half
   supports a dynamical TS; report binomial uncertainty and uncommitted shots.
   A small shooting pilot is diagnostic, not proof of a precise committor.
4. **Compare templates under a matched design budget.** Once the profile is
   accepted, propose a bounded pilot of 24 raw sequences per chosen bottleneck:
   12 from one verified saddle, when available, and 12 split among three representative ensemble
   templates. If five templates are justified, fix an equal total allocation
   before sampling. Use the same mutation mask, protected catalytic/maturation
   identities, sampling temperature and matched WT controls. This is a proposal,
   not authorization to launch a new campaign now. Treat all four sites as
   nested observations and tie equivalent alpha and beta sequence positions
   separately; instantaneous conformations need not be identical.
5. **Evaluate whole sequences across states.** Pool and deduplicate candidates,
   then score/evaluate every candidate on all retained templates and controls.
   Compare balanced-state average performance and poor-state failures, retaining
   a diverse Pareto shortlist. LigandMPNN `--pdb_path_multi` batches independent
   inputs; it does not implement a joint multistate objective. Do not require
   identical sampled sequences across templates or splice position-wise consensus
   residues without evaluating the resulting whole sequence. Its coordinates/
   element context and likelihoods do not supply QM charge redistribution or
   an activation energy. Verify that reactive atoms and fixed catalytic context
   survive parsing; keep proton/bond inventories in explicit metadata [R10].
6. **Test catalytic improvement physically.** Independent prediction/relaxation
   and matched replicated MD test structural robustness, not barrier lowering.
   Reserve matched mutant-versus-WT QM/MM PMFs for a small shortlist, comparing
   each sequence with its own reactant basin. The target is
   `DeltaDeltaG‡ = DeltaG‡(mutant) - DeltaG‡(WT)`, not an absolute energy or
   LigandMPNN-score difference between sequences. Check other relevant barriers,
   deacylation compatibility and substrate/product access so one improved step
   does not create a downstream trap. Native-Thr crystal similarity stays a
   diagnostic, not a hard rejection gate.

**Precedent versus proposal.** Lauko et al. used a QM-derived reactive motif
and PLACER analysis across catalytic states [R27]; Choi et al. report PLACER
ensembles for TS1/acyl/TS2/apo states [R29]. These learned structural ensembles
are not Boltzmann samples or verified QM/MM transition-state ensembles.
Zoi, Antoniou and Schwartz explicitly analyzed 15 uncorrelated, committor-defined
TS structures per enzyme in a theoretical enzyme-engineering study [R53]. This
is direct precedent for multiple TS configurations, but not a demonstration of
our proposed LigandMPNN workflow or experimental nylon-hydrolysis improvement.
Use primary links recorded in `../manuscript/REFERENCES.md` for these claims.

Do not compress the following into an arbitrary weighted activity score.
An attractive secondary score cannot rescue a failed chemistry gate. Missing
calculations remain incomplete; scientific failures remain in the denominator.

| Priority | Evidence | Planned use |
|---|---|---|
| Required | Sequence, maturation, assembly and molecular graph | Protected identities; separate alpha/beta ties; correct termini, charges, stereochemistry, covalent state and explicit water H; interpretable complete topology; no severe unresolved clashes |
| Primary | Independent ES recovery and catalytic geometry | Without covalent enzyme–substrate bonds or imposed catalytic contacts, assess d1/d2/d3, a_BD and projected Herschlag phi after the declared preparation/relaxation; also report Thr chi1/psi, crystal-aligned Thr N/OG1 displacements, and CB–OG1–substrate-C as orientation diagnostics; retain full distributions and raw pass/attempt counts across seeds |
| Primary as references become available | State compatibility and access | Check supported intermediates and acyl-enzyme/deacylation-water states as well as apo/ES; use each state's chemistry and geometry; screen substrate entry and product exit rather than optimize only a closed intermediate pocket |
| Supporting | Local AF3 confidence and consistency | Ligand/reacting-group confidence, catalytic-region confidence and ligand–protein/assembly PAE or interface metrics with explicit chain selections; report seed dispersion and pose recovery, not global pLDDT alone |
| Supporting | Physical quality and structural deformation | Local clashes, packing, unsatisfied polar groups, preserved interfaces; matched WT/self-repeat effective-strain baseline and regional maps; useful deformation is not automatically a failure |
| Optional, secondary | Noncovalent interface_delta | Include only if a matched WT/mutant pilot shows useful information beyond other metrics; use identical chemistry/scoring conventions and do not call it binding free energy |
| Diagnostic | Sequence likelihood and iteration convergence | Track position-wise frequencies, diversity, repeated sequences and cycles; use to monitor sampling, not as evidence of activity |
| Final selection | Diversity and experimental observables | Choose complementary passing sequences from tradeoffs rather than many copies of one motif; assay expression/maturation and detected product identities as well as activity |

Use the existing ES Tier A/B implementation only with the verified substrate
atom map and geometry convention. Near-attack ES cutoffs do not apply to a
tetrahedral carbon, a covalent bond or a TS. No universal confidence cutoff is
declared for the very-low-confidence substrate. Initial pilot results and
matched controls must establish which confidence quantities are informative.

Here, reproducible catalytic geometry primarily means independent **noncovalent
ES** near-attack recovery, without a covalent bond or imposed catalytic contacts.
TI1 is a separate conditional test: retaining its imposed covalent linkage is
insufficient. Check tetrahedral stereochemistry/bond geometry, oxyanion contacts,
protonation, any explicitly included relay water and local physical quality with
QM/MM; a classical TI1 test additionally requires validated state-specific
parameters. Classical TI1 parameters are not a prerequisite for QM/MM sampling.
A true PES TS1 has a separate saddle/Hessian/endpoint-connectivity requirement; prediction
agreement alone cannot establish a transition state. These are distinct gates,
not interchangeable measurements of independent binding.

Assess every ligand against all four equivalent sites. Report all site results,
with one model-level any-site recovery count for a single ligand; a ligand at
BF instead of AE is not a failure. Avoid double attribution with multiple
ligands. Sites within a model are nested observations, not independent seeds.
Report sensitivity to occupancy and guidance in separately labeled arms.
Neither independent pose recovery nor its seed fraction is a binding affinity,
occupancy probability or catalytic rate.

## Energetic comparisons and effective strain

The inspected legacy helper
`/data/jpc/code/redesigner/scripts/amber_relaxation/pyrosetta_interface_delta.py`
calculates E(complex) − E(extracted protein) − E(extracted ligand), without
separate reoptimization. Its values are Rosetta score units, not experimentally
calibrated binding free energies. It handles one residue named LIG and silently
falls back to the last residue if LIG is absent. Do not reuse that fallback;
a future adapter must reject missing/ambiguous ligands and validate parameters,
assembly, ligand selection and scoring terms. Report restraint contributions
separately; do not reward a design for satisfying imposed contacts.

Use conventional separated-interface scoring only for noncovalent ES (or
separately defined noncovalent product complexes). For TI1 or acyl enzyme,
extracting the ligand breaks a covalent chemical graph and is not the same
energy comparison. Evaluate these states with validated parameters and local
geometry/interaction analyses; do not compare absolute total energies across
different sequences or proton/bond inventories as reaction barriers. Later
barrier calculations must compare each sequence's corresponding reactant basin
and TS with balanced atom/proton inventories and the same protocol.

Effective strain denotes structural deformation relative to a reference; it is
not force-field strain energy. Keep those metrics separate. Establish the
same-sequence prediction variability baseline before interpreting WT–design
differences. A coordinated redesign may make a useful local rearrangement.
Backbone correspondence and regional confidence must be explicit, especially
before any future indels. Low deformation alone does not establish good binding.

## Iteration and stopping

Keep parent backbone, starting sequence, mutation mask, designed sequence,
round, random seed, state, guidance, MSA/template settings and validation
outcomes in a per-attempt ledger. Deduplicate sequence identities while retaining
their sampling counts and lineage. Do not treat repeated proposals as
independent successes. Track whether geometry/robustness and the Pareto set
improve, alongside sequence changes and diversity. Sequence convergence may
reflect sampling settings, protected positions or a model's preferences.

Use a bounded pilot and an independent final evaluation with seeds/settings
reserved from iterative selection. Freeze a stopping budget and plateau rule
after the pilot; stop when the budget is reached or independent accepted
metrics plateau, not merely when two rounds emit the same sequence. No numeric
budget or plateau threshold has been selected yet. Retain diverse alternative
lineages and evaluate experiment rather than indefinitely optimizing predictor
self-consistency.

## Verified LigandMPNN symmetry behavior

The official README and run.py were checked on 2026-10-22 [R10; links below].
`--homo_oligomer 1` derives ties from matching residue identifiers across all
parsed chains with equal weights, overriding explicit tie lists. This is
appropriate for genuinely identical chain types, not the eight mature Nyl12
chains of unequal alpha/beta lengths. For our (alpha-beta)4 representation use:

```text
--homo_oligomer 0
--symmetry_residues "A113,B113,C113,D113|E40,F40,G40,H40|..."
--symmetry_weights "0.25,0.25,0.25,0.25|0.25,0.25,0.25,0.25|..."
```

This is a schematic excerpt, not a runnable complete mutation mask. Generate
all intended tie groups from the checked correspondence map, retain fixed
identities and verify final alpha copies and beta copies independently. The
recent October 16 protected preparation already uses this explicit approach.
The older precursor-chain homooligomer workflow represents different chain
boundaries and must not be copied unchanged.

These options enforce sequence equality; they do not request a spatial point
group. The retained scaffold's independent repeat audit supports D2-like
geometry, with three near-180-degree rotations. Experimental biological
assembly operators still need checking before a production RFD3 symmetry
choice; a crystal space group is not that assignment. No historical explicit
point-group declaration has yet been recovered in the bounded source search.

Primary sources:

- [LigandMPNN README](https://github.com/dauparas/LigandMPNN/blob/main/README.md)
- [LigandMPNN implementation](https://github.com/dauparas/LigandMPNN/blob/main/run.py)
- [LigandMPNN paper](https://doi.org/10.1038/s41592-025-02626-1)

Installed implementation audited at commit
`26ec57ac976ade5379920dbd43c7f97a91cf82de`. Direct publisher/author-PDF
retrieval failed during this review; publisher-indexed primary-paper abstract
and result excerpts support the limited method/validation statements above.
No new design, AF3, Gaussian, force-field-fitting or refinement job was
launched to prepare this document. The ranker and matched evaluation arms
remain to be implemented.
