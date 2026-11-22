# Nyl12 sequence-engineering pilot

Status: PI approved a single-parent, 24-sample Model23 LigandMPNN pilot on
November 21, 2026. WT40/41 fail the native-Thr seed checks; among seven Tier
A/B WT models, only Model23 passes occupied-site native-Thr, observed-reference
identity and overlap checks. The CPU tmux pilot completed successfully: 24 samples, eight unique sequences,
16 duplicates and four to six substitutions per precursor. All symmetry,
protected/fixed-position and tensor/FASTA checks passed. Results and mutation
sets are in `design/nyl12_ligandmpnn_pilot/README.md` and
`accepted_parent/run/results/`. The omitted bundled OpenFold startup failure
is preserved; the complete-runtime run used existing installed support without
adding dependencies. Candidate structural/functional validation remains open.
Existing QM/MM and US continue unchanged.

Use **violet and azalea only**. Geranium, tulip and aster were inspected and
are occupied by other users; no transfer or job was launched there. Preserve
current QM/MM, US and WT MD first, followed by conditional pilot execution and
the already queued post-MD analysis. No new mechanism branches/LLM benchmarks.
The existing violet CPU runtime is verified and can run the pilot on one
low-priority CPU with GPUs hidden, with the accepted parent specification.

Preflight evidence: `design/nyl12_ligandmpnn_pilot/native_preflight/`.
WT40 occupied-site OG1/N shifts are 1.261/2.149 Angstrom, psi delta -56.03
degrees; WT41 has 1.480/2.530 Angstrom and -38.17 degrees. Model23 has
0.188/0.705 Angstrom, psi -30.09 and chi1 +4.52 degrees. Cutoffs remain
0.75/1.25 Angstrom and 35/20 degrees. Other copies are measured diagnostics,
not required to share identical instantaneous geometry with the occupied site.
All four copies must retain identity and structural integrity. The CIF has 313
observed author-numbered residues, with unresolved positions 1–7, 228–232 and
326–328; it is not a different 313-residue construct to align end-to-end.
The earlier Model23 QM/MM deterioration remains a negative result; passing the
static seed checks neither erases it nor certifies catalysis.

## What the first calibration establishes

`analysis/nyl12_filter_calibration/` retains a scripted, source-frozen static
audit of 900 relaxed Boltz-2 models: WT and eight mutants, 50 models per variant
in each of two separate PA6 and PA66 modeled-register campaigns. All historical
ensemble counts reproduce, with zero Tier A/B recomputation mismatches.

Every variant retains at least one A/B model in each campaign. This panel does
not reproduce complete variant rejection by an any-passing-model rule. It does
show why model pass fraction is not established as an activity-ranking score:
PA6 WT retains 7/50, K59D 18/50 and V111C 1/50. The corresponding legacy L2
signals were 1, 3.79 and 4.77 relative to WT; these are noisy product signals,
not verified rate enhancements. All PA6 mutant median interface scores are
less favorable than WT in this panel. Neither score should eliminate designs
as if it were a validated predictor of polymer hydrolysis.

The Burgi-Dunitz angle gate rejects 318/450 PA6 models; in 186 models it is the
only failed inclusive-B criterion. This identifies a preparation/pose/geometry
diagnostic, not evidence that its cutoff should simply be widened.

These results concern static historical relaxed models, not the newer matched
RESP/GAFF2 MD ensembles. Confidence/strain and independent mutant MD comparisons
are not validated by this geometric-count audit.

## Experimental evidence must stay provisional

The PI confirms crude-lysate, MS-based assays and considerably discordant
duplicate **K59E** experiments. The workbook is
`/data/jpc/enzyme-substrate-modeling/config/cb150_nyl12_ssl_sites_1_4_pa6_pa66_pmaxa.xlsx`.
Keep both source entries, matched WT/blank controls and product identities.
Do not choose a K59E run, average silently, fit thresholds or treat a derived
fold difference as a clean catalytic phenotype. Determine whether repeats are
independent experiments, biological replicates or repeat MS measurements.

The curated static/presentation analyses hardcoded one K59E value per product.
An older workbook-based MD analysis retained both entries and their arithmetic
mean; neither operation establishes a reliable underlying activity. The new
workbook audit will retain individual values without choosing or averaging.
Raw MS intensity/background and active-enzyme normalization may remain absent.

The two K59E PA6 L2 workbook values are 0 (B5) and 1.1465832774 (B6): their
mean is the legacy 0.5732916387. Thus the hard-coded legacy value represents
an average, not evidence that one experiment was selected. The disagreement
must remain visible. The workbook is a relative-product table, not raw MS
spectra or a complete replicate-level quantification archive.

The uploaded commercial-polyamides manuscript, Figure 7 and Methods, supplies
useful experimental context: independently re-screened hits (n=3 stated),
25 microliters lysate, 5 mg polymer, 200 mM Tris pH 8, 75 C for 6 h. It reports
background subtraction, propranolol normalization, product-specific 12-point
calibration curves and triplicate MS measurements per sample. These controls
do not normalize mature active-enzyme concentration or establish which source
rows correspond to the different replication levels. The newer PI-provided
Figure 7 image and the older cb150 workbook remain distinct source versions
until their numerical relationship is confirmed. The PI believes they use the
same data but cannot confirm with the experimentalist for a few days; do not
assume new measurements or delay structural analysis on that account.

The figure supports evaluating product profiles: F101E/F101H/P105K show large
PA66 monomer responses without comparable dimer responses. K59D/V111C show
increases across several PA6 products. These are reported substrate/product
effects, not proof of a shared change in intrinsic catalytic rate. Do not add
WT-relative folds across products; total measured-product release requires
the corresponding absolute concentrations and consistent units.

F101E/F101H have zero legacy L2 signal but positive alternative PA6 L3 signals;
they are not confirmed inactive controls. Polymer activity is not proof of
isolated L4 turnover. No specificity, ROC/AUC or prospective accuracy claim is
supported by this small, noisy positive-control panel.

## Bounded first design batch

| Item | Proposed specification |
| --- | --- |
| Generator | Existing LigandMPNN implementation/checkpoint, fixed backbone and WT chain lengths |
| Context | Complete mature tetramer, one declared PA6 ligand, all equivalent sites analyzed after prediction/relaxation |
| Parent | WT Model23 only; WT40/41 failed the occupied-site native-Thr gate, and no second parent is included |
| Budget | Two temperatures (0.1/0.2) x two independent RNG seeds x six samples = 24 raw samples; report duplicates and actual unique count |
| Editable positions | Alpha precursor positions 50, 51, 52, 73, 76, 101, 111, 112, 114, 115, 116, 151, 152, 153, 185, 187, 188; beta positions 39, 41 |
| Sequence ties | Explicit A/B/C/D and E/F/G/H ties, each with weight 0.25; heterodimer families and AE/BF/CG/DH sites |
| Mutation reporting | Count unique precursor positions, not four symmetry copies; aim for modest changes without assuming fewer changes means higher activity |
| Controls | WT and all eight characterized substitutions remain in retrospective comparison; retain K59E's inconsistent experiments explicitly |

The accepted manifest is `design/nyl12_ligandmpnn_pilot/accepted_parent/run/preparation_manifest.json`. The native assembly is four alpha/beta repeats, not eight identical chains. Use explicit `--symmetry_residues`/`--symmetry_weights` and `--homo_oligomer 0`. Preserve alpha length 232 and beta length 96. Protect alpha 113, 156, 186, 222–232 and beta 1–3, 40, 42. No indels are allowed. This conservatively
protects oxyanion-hole/maturation context, Thr233, Lys156 and Asp272/Asp274;
it does not assert that every protected residue has a proven catalytic role.

This is a single-parent E-S exploratory sequence pilot, not a validated-TS design campaign. Candidate relaxation, structure prediction, matched MD and experimental activity remain unvalidated downstream gates. Working preparation code already exists in the older repository:
`scripts/prepare_nyl12_focused_redesign.py` and
`scripts/lib/nyl12_focused_redesign.py`, with tests. Reuse it after source/mapping
review. Model 23 remains an informative negative preparation/QM-MM control;
its static Tier A status is not a reason to label it a validated design seed.

The manuscript's experimentally screened positions (24 positions within a
7-Angstrom modeled pocket) provide a position-level prior. Its alpha 95–106
loop, containing F101/P105, is a particularly relevant later remodeling region.
Retain the mapping caveat: precursor 289–300 maps to beta 57–68, while manuscript
chain letters refer to a different representation than the split A–H assembly.
These regions are candidates for review, not newly authorized backbone masks.

Hard checks remain identity, fold/assembly integrity, correct substrate graph,
valid preparation and catastrophic clashes. Treat confidence, strain, Tier
occupancy, interface scores and sequence likelihood as descriptive/Pareto
information until controls support stronger use. Compare model budgets and
report every rejected control and missing field. Final selection must include
sequence diversity and retain substrate-access/product-release compatibility.

## Later scaffold-constrained theozyme comparison

An optional later route is local active-site/backbone remodeling inside the
native fold, using fixed catalytic/maturation anchors and retained oligomer
context. RFD4 is the candidate generator, subject to its actual supported
conditioning and this assembly's mapping checks. Begin with exact WT lengths;
do not enable indels or redesign maturation in the first test.

Formal tetrahedral/acylenzyme geometries can be constructed as explicit design
hypotheses without claiming a stable TI1. They cannot supply missing evidence
for a particular proton relay or a true saddle point. Generic cluster motifs
and alternative proton-delivery geometries require embedded-pocket validation.
Do not require an artificial TI1 well or optimize binding of one intermediate
at the expense of adjacent TSs and product release. The scaffold does not become
mechanistically valid merely because a generator confidently rebuilds it.

This remains secondary to the current calibration-first sequence route. Slow
or incomplete polymer conversion alone does not localize the bottleneck to
chemical catalysis; substrate access and mature active-enzyme abundance matter.

## Future mask expansion proposal — not authorized for launch

The PI proposed a near-global editable mask while protecting catalytic
residues. A broader search is reasonable after the current 45-model validation
continues, but fixed-backbone LigandMPNN does not itself generate backbone
plasticity. Stage expansion from the current mask to second-shell positions,
then consider more distant positions only after matched WT calibration.
Evaluate an ensemble of open and closed backbone states. Protect the
maturation/catalytic core and essential assembly interfaces; tie equivalent
alpha copies to each other and beta copies to each other. Symmetry ties are
necessary for the intended assembly design but do not guarantee assembly
stability. No expanded campaign is authorized or queued; the existing 45 models
continue under their current scope.

## Separate future ideas

- **NylB/poly(MAXA):** investigate a NylB-derived catalyst that releases MAXA
  selectively. Define substrate graph, repeat boundaries, desired scissile
  bonds and competing products before modeling. Novel product intent is not
  evidence that Samsara has not investigated it or that existing claims do not
  apply. Keep it as a separate hypothesis requiring institutional IP review.
  Define exclusivity using product standards, detection limits and mass balance;
  a targeted MAXA signal alone cannot exclude other products.
- **Mechanism-agent:** retain reproducible tool adapters, evidence indexing,
  source/charge/mapping checks and human decisions. The tested local models
  are not qualified for autonomous mechanistic planning. Bounded extraction,
  reporting and proposals remain useful; stronger model claims require held-out
  qualification. No new model installation or autonomous dispatch is authorized.


## Approved structural validation, 2026-11-22

WT and all eight unique variants now have a matched 45-model Boltz/Amber
campaign: four zwitterionic PA6 L4 copies per tetramer and five shared prediction
seeds per sequence. See `design/nyl12_ligandmpnn_pilot/boltz_validation/README.md`
and its live status. Violet/tulip predict; geranium relaxes with released
RESP/GAFF2, ff19SB, explicit TIP3P/150 mM salt and neutral Thr / two anionic Asp
baseline. Generation contacts are released; all atoms move in minimization.
Graph-based ligand renaming precedes parameter assignment. Four site metrics
are nested within each seed/model. This implements the previously pending
prediction/relaxation gate, without authorizing new candidate MD or activity
claims. Existing mechanism/US jobs remain unchanged.
