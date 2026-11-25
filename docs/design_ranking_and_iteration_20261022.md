# Design ranking, iteration and symmetry — 2026-10-22

Status: proposed evaluation protocol, not an implemented ranker or an executed
design campaign. This extends the October 16 focused-redesign protocol in
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
