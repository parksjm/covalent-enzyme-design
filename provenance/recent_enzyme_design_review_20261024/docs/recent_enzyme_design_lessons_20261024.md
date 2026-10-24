# Recent enzyme design and the Nyl12 workflow — review 2026-10-24

## Decision and limits

Do not make exhaustive resolution of every proton-transfer hypothesis the
prerequisite for exploratory Nyl12 redesign. Establish a chemically plausible,
protein-compatible local catalytic model, retain explicitly labelled competing
mechanisms where needed, and screen compatibility across several states.
Reserve expensive saddle-point/barrier validation for reference mechanisms and
a shortlist. This is a proposed staged strategy, not a claim that we already
have a validated TI1, TS1, force field or improved enzyme.

The pseudo-ring is excluded only within the tested rigid AF3-derived pocket,
contact windows and steric criteria. It is NOT experimentally ruled out, nor
is a transient solvent-mediated proton relay ruled out. The finite front and
opposite searches do not identify the operating reaction mechanism. See
water_relay_mechanism_evidence_20261024.md for primary sources and
front_water_feasibility_20261024.md for numerical evidence.

## Verified 2025/2026 precedents

### Lauko et al., serine hydrolases — Science 2025, DOI10.1126/science.adu2454

The final article uses RFdiffusion, sequence/structure refinement and PLACER
ensembles across catalytic states. Early acylating designs often stalled at
the acyl-enzyme; later selections included apo/acyl-enzyme preorganization.
The useful lesson is multistate screening and measuring sustained turnover,
not simply obtaining a stable covalent pose. This is ester hydrolysis by a
Ser/His/acid system, not Nyl12 amide hydrolysis. The 2024 ChemNet preprint
(PMC11384011) differs from the final PLACER article (PMC12288761); do not mix
numerical activities or method names between versions. DFT basis details
from the preprint are not asserted as verified final-paper specifications.
Sources: https://pubmed.ncbi.nlm.nih.gov/39946508/ and
https://pmc.ncbi.nlm.nih.gov/articles/PMC12288761/ .

### Kim et al., metallohydrolases — Nature 2026 issue, online 2025-12-03

DOI10.1038/s41586-025-09746-w. DFT-generated catalytic configurations were
scaffolded with RFdiffusion2; LigandMPNN/Rosetta iteration and PLACER/structure
prediction assessed organization. Multiple catalytic arrangements were explored,
and a later campaign explicitly included the catalytic base. This shows the
value of testing alternative local chemistries and revising the catalytic motif;
it does not show that diffusion discovers an unknown mechanism. Small activated
ester/metal chemistry is not a benchmark for nylon cleavage. The 2026 issue
date must not be confused with a new 2026 online study.
Source: https://www.nature.com/articles/s41586-025-09746-w .

### Choi et al., cysteine proteases — Nature, online 2026-10-14

DOI10.1038/s41586-026-11143-w. The workflow combines backbone/motif refinement,
sequence design and AF3 assessment. PLACER compares apo, TS1, acyl-enzyme and
TS2 models. This is a strong multistate precedent; the inspected material does
not demonstrate that every TS-labelled input was an independently optimized
QM saddle. Primary extended data and official pipeline support these claims.
Sources: https://www.nature.com/articles/s41586-026-11143-w and
https://github.com/hojaec/RFD2_MI_cysteine_protease_pipeline .

### Chen et al., metalloproteases — Nature, online 2026-10-14

DOI10.1038/s41586-026-11167-2. A relevant hypothesis-building precedent is
motif derivation from a natural enzyme bound to a transition-state analogue,
rather than a newly located QM saddle. That detailed provenance was verified
in the authors' preprint (aminopeptidase N, PDB4QHP), not final SI. Its peptide
recognition problem is closer to polyamide register selection than ester
hydrolysis, but it does not validate nylon cleavage or polymer accessibility.
Sources: https://www.nature.com/articles/s41586-026-11167-2 and
https://www.biorxiv.org/content/10.1101/2025.11.20.689622v3.full .

### Butcher et al., RFdiffusion3 — Science, online 2026-10-15

DOI10.1126/science.aeh2982; primary metadata:
https://pubmed.ncbi.nlm.nih.gov/42879533/ . Final abstract/metadata were
verified, but full final Methods/SI were inaccessible; do not substitute
2025 preprint experimental results. Official enzyme instructions describe
selective atomic motif constraints, donor/acceptor and exposure conditioning,
and removal of hydrogen atoms during preprocessing. Thus RFD3 can support
backbone exploration around a supplied chemistry hypothesis; it does not
resolve explicit proton relays or compute electronic barriers. Documentation:
https://rosettacommons.github.io/foundry/models/rfd3/tutorials/advanced_enzyme_design_tutorial.html .

### Other groups and polymer redesign

Listov et al., Nature 2025, DOI10.1038/s41586-025-09136-2, combine catalytic
motifs with natural-fold design and FuncLib stabilization, including active-site
and core evaluation. For Nyl12, the transferable principle is improving the
packing and preorganization of an existing architecture rather than replacing
it with an unrelated fold. Their Kemp chemistry does not establish amidase
performance. Source: https://www.nature.com/articles/s41586-025-09136-2 .

Wang et al., Communications Biology 2025, DOI10.1038/s42003-025-08364-6,
redesigned PET-hydrolase loops with PRODA and replicated MD geometry indicators,
then tested polymer conversion. Their tetrahedral substrate model used geometry
from prior QM/MM work, illustrating reuse of a chemical template for redesign
rather than a new enzyme-specific saddle search for every variant. PET is a
polyester; its chemistry, material properties and assay cannot be treated as
nylon validation. Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC12174350/ .

## Proposed Nyl12 adaptation: separate chemistry from search breadth

1. **Reference chemistry, a bounded task.** Compare the native n12c structure,
   the predicted pocket and published NylC catalytic networks. Map residue
   identities before assigning Asp272 a general-base role. Specify the
   protonation/charge and missing environment explicitly. Test a small number
   of justified mechanisms, not an unrestricted inventory of proton paths.
   An accepted geometry is required before fitting; raw continuum Hessians or
   solvated charges do not automatically become ff19SB/GAFF2 parameters.
2. **Two uses of TS models.** For exploratory structural screening, a validated
   catalytic intermediate or approximate transition-state-like template can
   define an ensemble of local motifs. Label it as such. For quantitative
   barrier comparisons and a claim of a true TS1, retain the intended Hessian,
   reaction-mode and endpoint-connectivity validation. Neither a covalent
   prediction nor PLACER/RFD3 supplies these electronic-structure tests.
3. **Keep the existing fold first.** Begin with focused LigandMPNN designs on
   WT-length native backbones, preserving N-terminal maturation, catalytic
   topology and alpha/beta interfaces. Tie equivalent alpha and beta copies
   separately across the tetramer. Optional RFD3 is a local, WT-length pocket/
   entrance-loop arm when the native framework cannot support the desired
   interactions; it is not a compulsory whole-protein redesign step.
4. **Multistate ranking.** Assess independent noncovalent ES, TI1, acyl enzyme
   and a deacylation-water/TI2-compatible model. Include apo/open-compatible
   and product states. For early screens these can be geometrical models;
   they are not all asserted optimized minima. Reject severe incompatibility
   with a later state. Compare design ensembles against WT and matched single
   mutants, reporting distributions/frequencies rather than one best seed.
5. **Do not freeze a polymer into a small-molecule pocket.** Use a local amide
   fragment for economical chemistry, then restore full assay-relevant oligomer
   graphs for specificity/access screening. Sample independent conformations
   on both sides of the reaction centre, alternative cleavage registers,
   termini and lengths. Evaluate core reaction-atom geometry separately from
   tail mobility; low distal confidence/global ligand RMSD is not by itself
   catalytic failure. Retain a subset of productively bound conformers without
   demanding a single rigid global substrate conformation. PLACER samples are
   structural predictions, not Boltzmann populations or binding free energies.
6. **Match what the assay detects.** PA6 L4-to-two-L2 conversion is a project
   hypothesis to test, not established by detecting L4 and L2 products. Compare
   single-terminus recognition and soluble, doubly charged oligomer scenarios
   separately. Include actual capped MAXA-related chemistry and off-target
   cleavage registers. Track soluble product identities/yields and turnover;
   oligomer activity is not proof of insoluble-polymer depolymerization.
7. **Preserve access and release.** Seek a locally organized catalytic core
   while retaining open-compatible entrance conformations. Do not reward
   rigidity of every loop/tail, strongest substrate binding, or sequence
   convergence as substitutes for catalysis. Use physical sampling only after
   the relevant state parameters pass validation; small-molecule design papers
   do not demonstrate polymer entry/exit or the dynamics of our symmetric enzyme.

## What can proceed now

Document native-to-literature residue correspondence and a small set of local
chemical hypotheses; prepare the multistate/oligomer/register evaluation panel
and its WT controls. These tasks need not wait for every mechanistic question.
Production relaxation, force-field fitting and claims of TS-guided improvement
remain gated by accepted chemistry and parameter validation. No new design
software, dependency, simulation campaign or Gaussian run was launched by this
review. Recent-paper practices inform proposals; they do not retroactively
validate previous failed clusters or replace the TS1 acceptance requirements.
