# Focused Nyl12 redesign workflow

Date: 2026-10-16. Canonical host: violet. Initial target: PA6 L4z -> L2.

## Objective and current implementation

Propose a small, diverse sequence panel that preserves Nyl12 maturation and
assembly and improves the probability of a product-compatible catalytic
arrangement. This is a prospective engineering workflow, not an addition to
the commercial-polyamide manuscript's activity claims. PA6 is the first
method-development target because it has the clearest MD-supported reference;
PA6,6 requires its own chemistry, atom mapping, controls, and calibration.

Implemented now: reproducible preparation of three audited PA6 backbones,
physical-chain PDBs, a common focused design mask, explicit LigandMPNN symmetry
constraints, executable argument vectors, hashes, and a diversity diagnostic.
No sequences, new predictions, covalent parameters, quantum calculations, or MD
have been run for this campaign. Later stages below are specifications, not
claims of implemented adapters or completed validation.

Entry point: `scripts/prepare_nyl12_focused_redesign.py`.
Configuration: `configs/nyl12_focused_redesign_pa6_20261016.json`.
Current preparation: `output_nyl12_focused_redesign_pa6_20261016/protected_inputs/`.
The output-root and `pilot_inputs/` snapshots are superseded because their
protection lists were incomplete; their manifests explicitly say DO NOT EXECUTE.

```bash
bash scripts/launch_nyl12_focused_redesign_preflight_20261016.sh
```

This launches scripted tests, protection audits, preparation, and runtime import
checks in tmux session `nyl12_focused_redesign_preflight_20261016`, with logs in
the campaign's `logs/preflight.log`. All later calculation stages must likewise
use scripted tmux workers. Existing sessions/destinations fail rather than
being silently overwritten. The manifest contains six LigandMPNN command
argument vectors; preparing them does not execute them.

## 1. Freeze structural hypotheses before sequence design

Start with the following Amber/OpenMM-relaxed Boltz WT models:

| Model | Role | Tier | d1 / d2 / d3 (A) | Attack angle / projected phi (degrees) |
|---|---|---|---|---|
| 40 | Best geometry under the existing geometry-only inspection ranking | A | 3.05 / 2.88 / 2.90 | 103.0 / 96.5 |
| 23 | Second Tier A arrangement | A | 2.88 / 2.89 / 2.87 | 94.4 / 99.0 |
| 41 | Existing replicated bound-start MD reference | B | 3.29 / 2.85 / 2.82 | 100.3 / -91.3 |

All come from `output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/final/`.
Model 41's three matched trajectories had Tier A/B occupancies of 37%, 60%, and
25%; this supports a useful reference, not perfect persistent catalysis.
After whole-protein C-alpha alignment, focused C-alpha RMSDs within 12 A of
the reference scissile group are 0.954 A (40/23), 0.864 A (40/41), and 0.697 A
(23/41). These quantify local variation, not proof of distinct free-energy
basins. Visual review of ligand conformation, loop states, and clashes is still
required. Do not count conformationally redundant starts as independent evidence.

Full tetramers are retained, with one bound substrate initially. Exported
design PDBs use physical alpha chains A-D (1-232), beta E-H (1-96), and ligand
I. The exporter inserts physical chain breaks; it does not invent a peptide
bond between the processed alpha terminus and catalytic beta Thr1. Original
coordinates, ligand atom names, and heavy-atom CONECT records are retained.
Hydrogens are omitted for LigandMPNN, not chemically reassigned.

Boltz's central scissile atoms are C54/O49 with leaving nitrogen N59. These
correspond to AF3 C12/O4/N3 in the prior mapped PA6 representation. Never copy
atom names between predictors or use the existing redesigner L5 defaults
(C66/O60). Preparation recomputes the three contact distances and verifies them
against the source metrics within 0.03 A, plus carbonyl/amide bond distances.
These checks are not a complete bond-order/chemical-graph validation.

## 2. Focused, symmetry-aware LigandMPNN proposals

Use the union of heavy-atom contacts within 8 A of C54 or O49 across the three
backbones, mapped to equivalent sequence positions. This is a candidate design
mask, not a declaration that every included residue is safe to mutate.
The recovered historical freeze policy is in
`/data/jpc/code/redesigner/archive/20260520_archive_cleanup/root_yaml/redesigned_residues_expanded_protected.yaml`,
lines 3-6. Its catalytic anchors are 113/186; its NylC-homolog maturation core
is 229/232/233/272/274, and it protects the entire 222-232 region. These were
conservative homology-based protections, not a demonstration that every position
has been experimentally established as essential in Nyl12.

| Merged position / identity | Physical-chain equivalent | Protection |
|---|---|---|
| Phe113 | Alpha A-D:113 | Catalytic geometry anchor |
| Asn186 | Alpha A-D:186 | Catalytic geometry anchor |
| Glu229 | Alpha A-D:229 | Historical maturation core; within protected tail |
| Asn232 | Alpha A-D:232 | Pre-cleavage alpha terminus |
| Thr233 | Beta E-H:1 | Catalytic N-terminal nucleophile |
| Asp272 | Beta E-H:40 | Historical maturation core |
| Asp274 | Beta E-H:42 | Historical maturation core |

Restore all those exclusions in all four copies, including alpha 222-232.
The additional beta 2-3 protection is retained conservatively. The initial
mask had mistakenly allowed alpha227 and beta40/42; no design job ran.
`scripts/lib/nyl12_redesign_protection.py` now fails on missing required
protections or incorrect identities in any copy. Each backbone records all
28 core identity checks, and the manifest hashes the recovered source YAML.

The initial common mask contains alpha positions
50, 51, 52, 73, 76, 101, 111, 112, 114, 115, 116, 151, 152, 153, 156, 184,
185, 187, 188 and beta positions 39, 41: 21 equivalent positions,
not 84 independent targets. Some, notably 101, are uncertain in the source
models and require explicit loop/structure review. Relaxed PDB B-factors are
zero; obtain confidence from mapped original predictions, not these PDBs.

Tie A/B/C/D choices at each alpha position and E/F/G/H choices at each beta
position using explicit equal weights. Do not use the automatic homooligomer
flag on this alpha/beta heterooligomer. Retain ligand context and fixed-residue
side-chain context. One ligand plus three empty sites is an occupancy
assumption; symmetry tying does not make their environments equivalent.

Plan three sequence samples at temperature 0.1 and three at 0.2 per backbone:
18 proposals before deduplication, not a large evolution campaign. Hash the
checkpoint, software, source structures, configuration, and emitted commands.
All positions outside the reviewed mask remain fixed. Audit the resulting
sequences for catalytic/protected identities, correct 232/96 lengths, tied
copies, and absence of unintended substitutions. Deduplicate by the complete
alpha/beta sequence pair, retaining every parent/backbone provenance link.
Report changed equivalent positions; mutation count alone is not a quality
score. A conservative subset mask can be frozen as a separate config if the
review finds that sampling all 21 positions is too disruptive.

Use WT, K59D, F101E, and F101H as matched diagnostic controls, recognizing their
lysate phenotypes do not isolate intrinsic catalytic activity. Keep a subset
of existing experimental data for validation rather than tuning every filter
against all eight mutants. New designs require eventual experimental feedback.

## 3. Independent structure prediction and noncovalent assessment

First test WT plus one proposed sequence from each parent backbone. Predict
apo and noncovalent ES states with five explicit independent seeds per sequence:
40 structures for this four-sequence smoke pilot. This verifies atom maps,
sequence handling, assembly, and whether the prediction protocol can recover
the reference geometry before expanding the panel. One diffusion sample per
seed keeps seed accounting simple. Identical seeds, sampling budgets, and
protocols apply to controls and candidates.

After that gate, at most 18 unique proposals plus four controls, two states,
and five seeds would require 220 predictions. Deduplication can reduce this.
Boltz is the initial backend to match established controls. AF3/AF3-LA is a
separate evaluation arm when needed; global confidence is not interchangeable
across backends. Record MSA/template use. A strongly templated fold or imposed
contacts reduce independence from the starting hypothesis.

Run a register-guided hypothesis arm when needed, but report it as guided.
Include a matched unguided arm on WT and a small shortlist to diagnose whether
the geometry depends on guidance. Restraint satisfaction is an implementation
check, not independent evidence for binding or catalysis. Do not discard all
nonproductive predictions before calculating a pass fraction.

Use existing noncovalent Amber/OpenMM relaxation and parameter audits, including
correct charge, atom mapping, and topology. Recompute geometry after removing
external catalytic-contact restraints; document any protein restraints retained.
Evaluate all generated samples. Explicit failures and missing outputs appear
in the completeness ledger; no candidate receives a final rank from an
incomplete or selectively filtered ensemble.

## 4. Gates and ranking: models, sequences, then parent arrangements

Avoid an arbitrary weighted sum presented as a catalytic activity predictor.
Hard gates establish interpretable models; Pareto selection identifies
candidates that trade off the remaining properties. Freeze definitions before
ranking new designs and retain the WT reference distribution.

| Level | Criterion | Decision rule |
|---|---|---|
| Hard gate | Sequence, maturation, assembly, and atom graph | Protected identities and symmetry intact; correct chains, register, termini, charges, bond orders, stereochemistry, and complete topology; no severe unresolved clashes |
| Individual ES model | Full catalytic geometry after relaxation | Existing fixed Tier A/B criteria, including both oxyanion contacts and both angular criteria; inspect lower-distance bounds/clashes separately |
| Sequence, primary | Robustness across seeds | Report A/B fraction with all attempted samples accounted for, plus the full geometry distributions; successful seeds are samples, not activity probabilities |
| Sequence, primary when available | Compatibility across catalytic states | Shared reacting-core placement, supporting interactions, and retained water access; align a fixed local protein core and compare correctly mapped atoms |
| Sequence, supporting | Local confidence and structural consistency | Use ligand/reacting-group, catalytic-region, and interface confidence plus seed dispersion, not global protein scores alone; confidence is not proof of correctness |
| Sequence, supporting | Physical quality | Local clashes, buried unsatisfied polar groups, protein packing, assembly contacts, and loss of key noncatalytic support interactions |
| Sequence, supporting | Relaxation/representation sensitivity | Check dependence on guidance, protonation, terminal groups, ligand occupancy, and source backbone for the shortlist |
| Secondary only | Energetics and effective strain | Same-state noncovalent interface_delta and local strain flag anomalies; do not optimize either as the primary activity objective |
| Parent backbone | Resequence robustness | Fraction of unique proposed sequences that pass frozen per-sequence criteria, with actual counts; six samples per parent are exploratory |
| Experimental panel | Diversity | Select complementary sequences/loop arrangements from the Pareto set, not many copies of one mutation pattern |

Tier A: d1 <=3.2 A, d2/d3 <=3.8 A, attack angle 90-120 degrees,
absolute projected phi 70-110 degrees. Tier B: d1 <=3.5 A, d2/d3 <=4.0 A,
attack angle 85-125 degrees, absolute projected phi 60-125 degrees.
Reuse `src/geometry_tiers.py`; the old external redesigner YAML has different
geometry conventions and L5 atom names and cannot be used unchanged.
These ES criteria do not apply to a covalent linkage or tetrahedral carbon.

Do not transfer arbitrary global pLDDT cutoffs to low-confidence ligands or
assume a confidence cutoff calibrated for protein C-alpha strain validates
substrates. Initially report confidence continuously and use gross uncertainty
to qualify interpretation; calibrate gates on matched controls. Lack of a
biologically calibrated score is an explicit limitation, not a reason to hide
failures. With small ensembles, report raw counts and avoid precise success
probabilities or unstable best-model statistics.

The current preparation implements neither the downstream Pareto ranker nor
multi-state prediction adapters. Their required inputs are seed-level sequence,
parent, state, guidance, chemistry version, geometry, confidence, and validation
records, including failed jobs. They should be implemented against this
explicit ledger rather than guessed from filenames.

## 5. Bounded covalent-state development: acyl-enzyme first

This work makes sense as a reusable capability, but its first deliverable is
one validated WT acyl-enzyme topology and relaxation, not a large covalent
prediction campaign. Compare a deterministic chemically built adduct against
a small AF3-LA prediction pilot. AF3-LA handles leaving atoms; it is not a
force field or a guarantee of physical geometry. Its code/weights are not
installed by this workflow. New production dependencies require confirmation.
Keep the existing contact-restraint fork separate and never modify
`/data/jpc/code/alphafold3`.

For the central PA6 scissile bond, explicitly partition the substrate graph
into retained acyl fragment and departing amine fragment. The stable acyl-enzyme
is a Thr-O-C(=O)-R ester after cleavage of the substrate C-N bond; the departing
fragment must not remain covalently attached. Beta Thr1's free alpha amino
group must remain chemically represented. Acylation removes the side-chain
hydroxyl proton. Select protonation states with explicit proton accounting;
do not transfer the previous deprotonated-Thr ES topology into this state.

The first tetrahedral intermediate retains the leaving C-N bond and has a
tetrahedral carbon/oxyanion with state-specific protonation and chirality.
Deacylation has different water-derived atoms and proton transfers. These
require separate graphs and parameters; never create them by adding a bond
to the intact ES graph. Do not apply serine-His-Asp triad assumptions to this
N-terminal Thr hydrolase. Write down the hypothesized proton relay and water
roles before considering transition-state design.

Parameterization and refinement sequence:

1. Freeze a named-atom graph, bond orders, hydrogen inventory, formal charge,
   state, residue boundaries, and full-system atom map. Build a capped model
   containing free N-terminal Thr and sufficient acyl context; put the peptide
   cap at the downstream backbone boundary, not on the catalytic amino group.
2. Generate several chemically valid starting conformers and preoptimize with
   GFN2-xTB and implicit water, using anchored or restrained boundaries as
   appropriate. A small active-site cluster can add oxyanion donors/waters.
   Existing `src/optimize_hydrogens_xtb.py` freezes ALL heavy atoms; it cannot
   repair a distorted heavy-atom adduct and is not this optimization protocol.
3. For a stable acyl minimum, inspect convergence, stereochemistry, and
   vibrational behavior within the appropriate constraints. A tetrahedral
   intermediate may collapse during unconstrained optimization; that is not
   evidence that stretching restraints should be increased until it stays.
4. Derive and audit a modified N-terminal Thr-acyl residue library and explicit
   junction parameters compatible with ff19SB and the chosen ligand force
   field. Include cross-junction bonds, angles, torsions, impropers, exclusions,
   and 1-4 scaling. Avoid double-counted atoms/charges and retain the peptide
   connection to beta residue 2. Derive charges with an Amber-compatible
   protocol; xTB Mulliken charges are not a substitute for RESP or a validated
   alternative. Audit cap removal and the total integer system charge.
5. Compare MM minima and relevant torsional profiles against the reference
   quantum calculations. Use higher-level QM/ESP/Hessian or torsion data when
   transferable ester terms do not reproduce the linkage adequately. A clean
   `parmchk` report only means terms were assigned, not that they are accurate.
6. Minimize the solvated full complex with staged protein/boundary restraints,
   releasing external catalytic-distance restraints. Inspect the linkage,
   surrounding backbone/side chains, water organization, and strain after
   release. Export the complete graph/CCD/Amber library/frcmod/topology/atom-map
   package with hashes, logs, and a reproducible validation report.

A higher-level QM/MM optimization is NOT required merely to construct the first
plausible stable acyl-enzyme model. It becomes appropriate when environmental
effects/proton transfers cannot be represented by the cluster, or when asking
about reaction barriers and mutation effects on chemistry. Classical MM
minimization and MD cannot establish those barriers. xTB optimization is an
economical preparation step, not quantitative validation of a catalytic barrier.
Check the protocol on a known experimentally characterized covalent enzyme
before treating Nyl12 results as general validation.

Standard Rosetta interface_delta is disabled for covalent states: separating
enzyme and ligand would break a chemical bond, and different reaction states
have different composition/protonation. Do not compare raw MM energies across
ES, acyl, and tetrahedral states as reaction free energies. An intact tether
during MD is not evidence of stronger binding or faster catalysis.

## 6. When MD equilibration is needed

Energy minimization is required for physical refinement. MD equilibration is
required before analyzing a solvated trajectory, and strongly recommended for
shortlisted ES/adduct structures to test local adaptation. It is not required
for every LigandMPNN proposal, cannot repair a wrong topology, and cannot
validate a reaction mechanism or measure a chemical turnover rate.

After topology/parameter validation, begin with WT plus two diverse finalists:
three independent velocity replicas per system. A provisional screening
schedule is restrained minimization, 0.1 ns NVT heating, 0.5 ns restrained NPT,
2 ns gradual restraint release, and 20 ns unrestrained production per replica.
Require stable density/temperature, no numerical failures, sensible linkage
geometry, and inspection of structural adaptation; extend equilibration if
needed rather than calling a fixed elapsed time equilibrium. All durations
are initial engineering budgets, not convergence guarantees.

Use the same full assembly, one-ligand occupancy, explicit solvent, salt,
protonation treatment, timestep, and restraints for WT and designs. The
historical matched protocol used ff19SB/GAFF/TIP3P, 150 mM salt, and 300 K.
Use 300 K for an initial method-comparison pilot; a separately labeled 338.15 K
(65 C) experiment is relevant to the purified-polymer assay. Do not change
temperature only for selected designs or pool the two ensembles. Reassess ES
Thr1/protonation rather than automatically inheriting the historical activated
OG1 microstate. New covalent parameters remain separate from ES parameters.

Nine ES trajectories of 20 ns give 180 ns of initial production. Acyl MD is
another 180 ns only after its independent parameter gate. Extend a shortlist
to 3 x 100 ns if it answers a specified question and the shorter pilot is
informative. Report replica-level productive geometry/contact distributions,
water access, loop states, drift, and assembly integrity. Do not count frames
or the four subunits as independent replicates, or interpret bound-start
occupancy as the equilibrium probability of spontaneous productive binding.

## 7. Practical execution, acceptance, and model choice

Preparation/provenance stays on violet. Future AF3/AF3-LA tests should begin on
violet's known-working Ada GPUs. Use azalea for Amber GPU MD; use its Blackwell
GPUs for prediction only after a small compatibility check. Run long stages
in tmux with one worker per available GPU, explicit locks/checkpoints, and
failure propagation. Freeze run manifests before execution. Never create
competing workers against the same output directory.

Near-term next action: review the three prepared structures and the 21-position
mask, validate the existing LigandMPNN runtime, then execute the six frozen
commands and audit the sequences. In parallel conceptually, develop a single
WT acyl graph/parameter test; no broad covalent or MD jobs are queued now.
New parameterization/QM software choices and production dependencies need
approval before installation. The previous note to extract a new reusable
AF3/OpenMM/Rosetta repository remains in force; do not silently replace or
rewrite the separate existing redesigner repository.

GPT-6.1 Sol Medium is a reasonable default for orchestration, adapters, tests,
provenance, and routine analysis. For the first chemical-state specification,
charge/boundary decisions, parameter audit, and QM/MM interpretation, use a
higher reasoning effort (Sol High or Astra) and ideally independent chemistry
review. This task allocation is scientific judgment, not a published benchmark
of force-field expertise. A stronger language model does not validate a force
field. Official OpenAI guidance describes Sol as the cost/quality balance for
complex coding and Astra as the highest-capability option; test representative
tasks rather than assuming Medium cannot handle them.

## References and supporting implementations

Verification on 2026-10-16: Ruff and shell syntax checks passed; 31 focused
tests passed with 93.25% aggregate coverage (100% protection-module coverage).
The scripted tmux preflight completed with exit status 0, prepared the real
three-backbone inputs with 21 designable equivalent positions, and confirmed
that the installed LigandMPNN CLI imports successfully. This is not an
inference/GPU smoke test. No sequence-generation, quantum, or MD job was run.
Failed-attempt logs were retained; `--restart-failed` explicitly retries only
a dead nonzero-status preflight pane and refuses live/successful panes.

- Z. Wu et al. Designing enzymes for new-to-nature chemistry and non-natural
  substrates with AlphaProtein Novo. bioRxiv preprint (2026), posted October 5.
  DOI: https://doi.org/10.64898/2026.10.01.756017. Mechanistic multi-state,
  multi-seed, and multi-resequence evaluation is the main transferable idea;
  its esterase results do not validate Nyl12 amide or polymer catalysis.
- J. Dauparas et al. Atomic context-conditioned protein sequence design using
  LigandMPNN. Nature Methods (2025) 22: 717-723.
  DOI: https://doi.org/10.1038/s41592-025-02626-1.
- Official AP Novo release, including AF3-LA leaving-atom handling:
  https://github.com/google-deepmind/alphaprotein-novo.
- Official LigandMPNN options and explicit symmetry examples:
  https://github.com/dauparas/LigandMPNN.
- xTB geometry-optimization documentation:
  https://xtb-docs.readthedocs.io/en/latest/optimization.html.
- Amber modified-residue tutorial:
  https://ambermd.org/tutorials/advanced/tutorial1/section1.php.
- OpenAI model selection and the exact requested model:
  https://developers.openai.com/api/docs/guides/model-selection and
  https://developers.openai.com/api/docs/models/gpt-6.1-sol (checked 2026-10-16).
