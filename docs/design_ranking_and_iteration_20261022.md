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

### Reference states required for the first design round — clarified 2026-10-25

LigandMPNN can technically condition on a supplied ES or covalent geometry;
it does not require a complete computed reaction cycle. For the intended
TS-guided first round, prepare a consistent WT reference set: (1) noncovalent
ES with its declared protonation; (2) TI1 with tetrahedral connectivity,
oxyanion and the intended terminal-NH3+ state; and (3) the connecting TS(s).
If a protonated-Asp/attack intermediate is a separate minimum, retain it and
the separate saddles rather than force the route into a single TS1 label.
Bias-stabilized endpoints are templates, not certified minima. Saddle acceptance
requires the intended unstable mode and downhill connection to adjacent states.

Use the existing apo/WT assembly as a structural control. An acyl-enzyme and a
product-bound/deacylation-water compatibility model should be included before
final candidate selection; full free-energy pathways for all later steps are
not prerequisites for exploratory LigandMPNN sampling. Preserve symmetry,
residue/atom maps, maturation residues and substrate register across states.
Distinguish early template-based sampling from validated TS-guided design.

Do not compress the following into an arbitrary weighted activity score.
An attractive secondary score cannot rescue a failed chemistry gate. Missing
calculations remain incomplete; scientific failures remain in the denominator.

| Priority | Evidence | Planned use |
|---|---|---|
| Required | Sequence, maturation, assembly and molecular graph | Protected identities; separate alpha/beta ties; correct termini, charges, stereochemistry, covalent state and explicit water H; interpretable complete topology; no severe unresolved clashes |
| Primary | Independent ES recovery and catalytic geometry | Without covalent enzyme–substrate bonds or imposed catalytic contacts, assess d1/d2/d3, a_BD and projected Herschlag phi after the declared preparation/relaxation; also report Thr chi1/psi, crystal-aligned Thr N/OG1 displacements, and CB–OG1–substrate-C as orientation diagnostics; retain full distributions and raw pass/attempt counts across seeds |
| Primary as references become available | State compatibility and access | Check TI1 and acyl-enzyme/deacylation-water states as well as apo/ES; use each state's chemistry and geometry; screen substrate entry and product exit rather than optimize only a closed TI1 pocket |
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
protonation, explicit relay-water geometry and local physical quality with the
validated TI1 force field, including sensitivity to preparation restraints.
TS1 has a separate saddle/Hessian/endpoint-connectivity requirement; prediction
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
