# Calibrate the filters before enlarging the sequence search

Prepared 2026-11-20 for the 2026-11-21 planning checkpoint. Documentation and
source review only: no model installation, weight download, inference, API
charge, topology change, or new scientific calculation was performed.

## Decision

If experimentally active single mutants fail the current Tier A model/MD
filters, the first question is whether the filters recover known activity.
More sequence generation does not resolve an unvalidated rejection rule.
Retain hard chemical validity gates; treat specific bound geometry and Tier A
occupancy as testable surrogates until retrospective controls establish their
usefulness. Do not optimize Tier A occupancy as the primary design objective
or widen its cutoffs merely to rescue preferred variants.

RESP/GAFF2 work, nonstandard beta-N templates, explicit protonation inventories
and pKa analysis remain valuable prerequisites. Their existence does not prove
the modeled substrate pose, protonation populations, dynamics, or relationship
between a geometric threshold and experimental turnover. Classical MD samples
the specified potential and protonation state; ordinary fixed-topology MD does
not itself simulate hydrolysis or establish a catalytic free-energy barrier.

The existing manuscript's kinetic caveat applies: polymer conversion in crude
lysate is not a direct measurement of kcat. Expression, soluble active enzyme,
maturation, thermal survival, adsorption/access to polymer, binding, chemistry,
and product release can all contribute. The documented assay temperature
(75 degrees C) also differs from the approximately 300 K historical MD. An
active label is nevertheless a necessary positive control for the complete
prediction/filter pipeline. A failed positive control challenges that pipeline;
it does not, by itself, identify which modeled component is wrong.

There is also a substrate mismatch to test explicitly. Activity on PA6 polymer
powder does not establish activity on isolated zwitterionic L4, cleavage into
two L2 molecules, or preference for our modeled register. The predominant
observed soluble product is L3. Consequently, a known-active mutant is a
positive control for the overall polymer-hydrolysis objective, not proof that
it must maintain one chosen L4 pose. Match assay labels to their actual
substrate/product endpoint before treating computational disagreements as
force-field or protonation errors. A defined-oligomer assay would distinguish
this uncertainty if it becomes experimentally practical; none is assumed here.

The PI clarified that the reacting chain is likely much longer than L4; short
models were chosen because long, flexible polyamides are difficult to model
with Boltz-2. Neither binding direction nor the distance of a physical polymer
terminus from the scissile amide is known. Treat L4 as a tractable local polymer
proxy, not evidence for terminal recognition or isolated-L4 turnover.

For a later bounded sensitivity comparison, two neutral continuation caps
represent an interior segment, while one free end plus an opposite cap tests
proximity to the N or C terminus separately. Keep the uncapped L4 model as the
historical reference. Chemically cap the chain rather than silently changing
the protonation of its existing ends: capping changes connectivity and can
change atom counts, total charge and parameter scope. Define cap structures,
central-amide mapping, register, microstate and charge/counterions explicitly.
The current validated PA6 release does not automatically certify these caps.
These are future comparisons, not new campaigns or changes to active jobs.

## Bounded retrospective calibration

Freeze the sequence numbering, measured products, assay labels and current
filter version before looking at alternative thresholds. Use every existing
matched experimental control, rather than inventing inactive labels for
unmeasured sequences. First rescore retained structures/trajectories; launch
new calculations only for a documented missing comparison.

| Control group | Experimental record | Matched computational comparison | What it tests |
| --- | --- | --- | --- |
| WT | Conversion, uncertainty, assay conditions and active-enzyme normalization if measured | Existing WT structures/MD and the same declared state/pose/protocol used for mutants | Baseline acceptance and reproducibility |
| Each confirmed active single mutant | Exact sequence and assay substrate/product; retain numerical measurements | Same mature assembly, ligand graph/register, protonation inventory, parameter release, relaxation protocol and analysis denominator as WT | False-negative rate among known positives |
| Each experimentally inactive/poor single mutant | Detection limit, expression/maturation/stability controls where available | Same protocol as above; preserve failures and exclusions | Specificity and whether apparently inactive means chemically inactive |
| WT and mutant across alternative starting poses | Same experimental label, not extra independent assay observations | A small common pose panel, including opposite local orientations where the signed scissile-amide mapping identifies them; independent seeds | Pose/start sensitivity versus sequence effect |
| WT and mutant across plausible microstates | Same experimental label | Explicit, charge-consistent states evaluated as separate hypotheses; paired state-by-state comparisons | Dependence on the protonation assumption |
| WT and mutant across relaxation choices | Same experimental label | One frozen primary protocol and one bounded sensitivity protocol applied to all controls | Whether preparation, rather than sequence, drives rejection |

Use representative ensembles, not just a best-scoring snapshot. Report
per-replica distributions, individual gate failures, joint occupancy, analyzed
frames/time, and confidence or uncertainty where support permits. Frames,
symmetry-related sites and multiple poses from one variant are not independent
experimental variants. An MD endpoint and a force-converged stationary point
have different acceptance requirements; a completed job is not proof of either.

Preserve chemical/product/connectivity gates: the modeled PA6 L4-to-two-L2
reaction must have explicit atom/product mapping, while remaining a modeling
hypothesis rather than an experimentally established precursor-product pair.
The experimental products are predominantly L3. TI1, when modeled, must retain
the scissile C-N bond and declare the Thr233 linkage, oxyanion and proton
inventory. Do not mix TI1, ES, acyl enzyme,
or precursor maturation chemistry. Preserve tetramer correspondence and
alpha/beta maturation mapping. Distinguish these requirements from soft
preferences such as one H-bond distance, angular bin or rotamer population.

Report how many known positives the current filter retains before changing it.
If positives are rejected systematically, use the failure decomposition to
choose one hypothesis to test: pose/register, assembly, preparation, protonation,
parameter transferability, limited sampling, or an overly specific geometric
surrogate. Parameter quality should be assessed with relevant chemical evidence,
not repeatedly refit simply because occupancy is low. A filter that cannot
separate controls can remain descriptive but should not eliminate designs.
If controls are scarce, report that limitation; do not fit a high-dimensional
classifier. Hold out variants when enough independent labels exist, and freeze
any revised rule before prospectively testing new sequences.

For single-mutant structural controls, a common mature X-ray-derived scaffold
with locally repacked mutant side chains is a useful matched comparison to
independently predicted AF3/Boltz models. It separates mutation effects from
unequal prediction/reconstruction quality. It is not an E-S crystal structure:
the substrate pose still needs an explicitly sourced model, and relaxation/
ensemble checks remain necessary. Use this as one bounded preparation-sensitivity
test rather than regenerating every variant with every available predictor.

## Current PLM roles and access

The old EvolutionaryScale repository now redirects to Biohub. The statements
below reflect official pages inspected on 2026-11-20, rather than the original
2024 release terms. Sizes and availability are source claims, not local runtime
measurements.

| Model | Verified role and input | Appropriate use here | Current local access |
| --- | --- | --- | --- |
| ESMC | Sequence masked language model with logits and embeddings; public 300M, 600M and 6B variants [S1] | Cheap sequence plausibility, mutation scores, representation/outlier analysis; evaluate as a secondary prior | Released weights; smallest variant is the proposed pilot, not a requirement to use 6B |
| ESM3 | Generates across sequence, structure and function tracks; partial prompts can condition generation [S2] | Optional later constrained sequence exploration or structure-conditioned comparison | Open small model is 1.4B; official README lists larger 7B/98B models through the platform, not equivalent local releases |
| ESMFold2 | Actual current all-atom structure predictor with optional MSA; explicitly accepts small molecules, modified residues, DNA/RNA and complexes [S3] | Optional independent structural/pose hypothesis comparison on a few controls | Released ESMFold2 and single-sequence ESMFold2-Fast weights; this is not the 2022 ESMFold model |
| ESM-IF1 | Protein inverse folding from backbone coordinates; official loader extracts N/CA/C [S4] | Protein-backbone compatibility comparison; no default explicit ligand-atom conditioning | Meta repository provides model usage and is archived; qualify dependencies before reuse |
| LigandMPNN | Protein sequence design with ligand atom context; atom-context masking is an explicit ablation [S5] | Retain focused, fixed-length, symmetry-tied ligand-aware generation as the first design baseline | Official code and model parameters are MIT |

ESMC sequence likelihood is neither binding affinity nor catalytic rate.
ESM-IF1's protein backbone conditioning does not acquire ligand chemistry just
because the input PDB contains a ligand. LigandMPNN explicitly sees ligand
context, but its sequence confidence is also not a measured affinity or barrier.
ESM3's function-conditioning track does not constitute a validated Nyl12 turnover
predictor. ESMFold2's ligand capability makes it relevant to pose hypotheses;
fold confidence and benchmark binding-design results do not validate PA6
hydrolysis chemistry or special proton inventories. These are project-specific
inferences from the documented inputs/objectives, requiring calibration here.

### License and resource qualification

The current Biohub repository explicitly states MIT for the released models,
and its license text permits commercial use subject to notice retention [S6].
The current ESM3 open-small model card independently states MIT [S7]. Do not
repeat obsolete blanket claims that current ESMC/ESM3 releases are necessarily
noncommercial. This is not a claim that every historical weight artifact was
relicensed, or that larger hosted ESM3 weights are downloadable.

Nuance: the current ESMC Hugging Face cards show both `mit` and `other` metadata
tags [S1], whereas the official repository's explicit license section is MIT
[S6]. Before a production download, pin the exact revision and retain its
license/model-card files; resolve any artifact-specific inconsistency then.
Platform API access has separate service terms/AUP and is not equivalent to a
local weight license. No API enrollment or payment decision is made here.

The ESMC card recommends GPU execution; no local wall-time or VRAM claim has
been established [S1]. The official file listings show approximately 1.33 GB
for ESMC-300M and 26.7 GB for ESMFold2 [S8]. Download size is not peak VRAM.
Sequence length, precision, batches, chain count and sampling settings matter;
benchmark one control before budgeting a campaign. The smallest sequence model
is a resource-saving inference, not a verified accuracy guarantee.

## Minimal PLM pilot, after dependencies and exact artifacts are approved

1. Freeze a FASTA/control manifest containing WT and all known active/inactive
   mutants, plus an independently curated homolog alignment if available.
   Keep the precursor sequence versus mature alpha/beta-chain representation
   explicit; never insert a synthetic chain linker silently.
2. Start with MSA conservation and ESMC-300M masked mutation log-odds:
   log p(mutant residue | masked WT context) minus log p(WT residue | same
   context). For single mutants this is a direct, inexpensive candidate score;
   define multi-mutant scoring separately because additive scores omit epistasis.
   Logits from an unmasked position must not be mislabeled masked likelihood.
3. Compute only a fixed embedding summary and distance to WT/validated homologs
   as an outlier descriptor. Whole-protein embeddings may barely move for a
   single substitution; no activity cutoff should be inferred from that alone.
4. Compare these scores with assay labels and the frozen geometry-filter results.
   Report rank agreement and known-positive retention; with few variants,
   present the table rather than a misleading trained activity model.
5. Continue only if the inexpensive scores add useful discrimination or flag
   independent sequence risks. Retain Pareto/diverse candidates rather than
   deleting every low-PLM-score variant. A beneficial non-natural catalytic
   mutation may score poorly under an evolutionary prior.

This pilot supplements the existing LigandMPNN baseline; it neither replaces
ligand context nor cures an unvalidated structural filter. The next expensive
calculation should answer a calibrated failure hypothesis, not merely increase
candidate count.

## Primary sources inspected

- [S1] [Biohub ESMC-300M model card](https://huggingface.co/biohub/ESMC-300M)
  and [ESMC-6B model card](https://huggingface.co/biohub/ESMC-6B).
- [S2] [Official ESM3 README](https://github.com/Biohub/esm/blob/main/_assets/ESM3_README.md).
- [S3] [Official ESMFold2 model card](https://huggingface.co/biohub/ESMFold2).
- [S4] [Meta ESM-IF1 documentation](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding)
  and [backbone coordinate loader](https://github.com/facebookresearch/esm/blob/main/esm/inverse_folding/util.py).
- [S5] [LigandMPNN official README](https://github.com/dauparas/LigandMPNN),
  especially ligand atom-context and symmetry options.
- [S6] [Biohub ESM README, licenses](https://github.com/Biohub/esm#licenses)
  and [MIT text](https://github.com/Biohub/esm/blob/main/LICENSE.md).
- [S7] [Biohub ESM3 open-small model card](https://huggingface.co/biohub/esm3-sm-open-v1).
- [S8] [ESMC-300M files](https://huggingface.co/biohub/ESMC-300M/tree/main)
  and [ESMFold2 files](https://huggingface.co/biohub/ESMFold2/tree/main).

Repository context: `docs/nyl12_covalent_resume_20261019.md` (state/product
definitions), `docs/MASTER_PLAN_NYL12_COVALENT_REDESIGN.md` (fixed-length
LigandMPNN baseline), `docs/redesign_backbone_symmetry_and_access_20261022.md`
(ensemble and symmetry controls), and `manuscript/RESULTS_DRAFT.md` (barrier
scale and assay interpretation). This note changes no scientific acceptance
record and certifies no existing topology or trajectory.
