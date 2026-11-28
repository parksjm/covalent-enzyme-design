# Nyl12 TS-ensemble redesign protocol

Updated 27 November 2026. This is the authoritative prospective protocol for
the bounded pilot. JSON acceptance receipts establish what has actually run;
the mechanism roadmap retains results and superseded approaches. The objective
is improved PA6 cleavage, without prescribing a particular oligomer product.

## 1. Select chemically coherent ensemble representatives

Use the relay PMF together with its physical coordinates, proton inventory,
overlap, statistical support and temporal/initialization checks. Cluster
decorrelated structures near the dominant free-energy maximum and retain
representatives of any earlier competing bottleneck. Preserve each source
window, block, frame, bias, time, weight and reaction-atom map. Do not select a
frame merely because its instantaneous potential energy is high.

Until sampling and reaction-coordinate qualification are sufficient, label
these structures **TS-like design templates**. A PMF maximum alone neither
establishes a first-order PES saddle nor a committor probability of one-half.
Conversely, a discrete TI1 minimum is not required for hydrolysis. A verified
PES TS can complement the ensemble; it is not interchangeable with a sampled
free-energy bottleneck. Preserve both categories in metadata.

The current baseline uses the two original chemically coherent MD-derived relay
representatives and their existing 12 LigandMPNN designs. Historical RFD3
remodels are preserved as diagnostics, outside the current comparison. This
does not claim a converged transition-state ensemble. Provenance is in
`design/nyl12_relay_ts_ensemble_20261126/templates/template_manifest.json`.
The 101-protection pilot has completed generation and scoring. Its
`original_backbone_triage_20261127/` reuses the two original TS-parent scores
to rank the 12 baseline candidates; no new sequences or predictions were
generated. The dated `md_baseline_decision_20261127.json` supersedes the
earlier universal symmetry-hold receipt without rewriting it.

## 2. Export, canonicalize and qualify before design

Convert matching Amber topology/restart/frame sources to canonical precursor
identities, retaining mature alpha 232/beta 96 chains and four equivalent copies.
Retain explicit canonical ↔ mature-local ↔ Amber atom/residue maps, complete
protein backbone, PA6 atom names/connectivity and source coordinates. Distinguish
the reactive substrate fragment from artificial QM link atoms. Verify the
native PA6 graph, charge/protonation inventory and fixed catalytic geometry.

Never infer identity from a historical Amber residue 929 offset. Before any
RFD3 export after Amber, restore canonical identity using the actual source
map. Existing parent exporters, graph audits and native LigandMPNN parser
checks are reused. Loss of a fixed atom or alteration of substrate chemistry
fails qualification; it is not repaired by silently changing the method.

## 3. Compare direct sequence design with optional partial diffusion

RFD3 is outside the current critical path. Review the existing 12 original
MD/TS-backbone LigandMPNN designs from the two chemically coherent MD-derived
parents, retaining the 101-position protection mask and mandatory four-copy
sequence ties. No RFD3 GPU qualification or 12-design replacement is required.
Identical sequences need not have identical thermal conformations; do not force
exact D2 coordinates onto an MD snapshot. Multiple chain conformations alone
also do not establish TS compatibility: retain reaction mapping, proton inventory,
ligand graph and catalytic geometry for each parent.

Fixed-backbone LigandMPNN design followed by independent template-free,
paired-MSA Boltz prediction and matched Amber relaxation can assess
sequence-driven repacking and assembly; it does not promise improved activity.
A future optional RFD3 arm may evaluate coupled perturbations of A/E, B/F,
C/G and D/H, with independent physical-movement, clash, assembly, fixed-context
and chemical checks. It is not implemented or qualified for production, and
exact coordinate symmetry is not a universal acceptance requirement.

No loop insertion/deletion or new capped-substrate parameter campaign belongs
to this pilot. See [the retained qualification audit](RFD3_SYMMETRY_QUALIFICATION.md).

## 4. Run a bounded, symmetry-tied LigandMPNN experiment

The immutable pilot sampled six sequences on each of four historical parents.
Review the 12 samples on the two original MD-derived parents; the other 12
remain diagnostic. Retain temperature 0.1, no native-identity bias and no forced
mutation count; no replacement generation is required.
Tie equivalent alpha positions across A–D and beta positions across E–H.
Preserve the current 101-position mask in all copies, including K156 and L235.
Keep the entire 221–232 alpha tail and the established catalytic identities.

The mask comprises 24 core policy positions plus 77 interface-only positions,
with no detected disulfide constraints. Interface protection is a conservative
assembly heuristic, not evidence that every protected site is catalytic or
individually essential. The exact rule, all 101 identities, per-position reasons
and source structures are in the [SI inventory and figure](figures/nyl12_sequence_protection_20261127/SI_caption_and_inventory.md).
The restriction can limit active-site loop redesign; future interface
co-design should be a controlled comparison, not a silent relaxation.

In particular, the strict contact set includes T53, F55, K59 and F101, where
previous lysate assays reported active variants. Those data are noisy (K59E
has discordant duplicates), but they preclude describing these sites as
universally immutable. The 101-position mask is a conservative baseline, not
a claim of the maximum biologically permissible design space. A later core-
only or interface-co-design arm should retain the catalytic/maturation
constraints and explicitly test assembly; symmetry ties alone do not
guarantee favorable protein–protein interactions.


The [primary-source Negoro audit](nylc_critical_residue_audit_20261127.md)
checks catalytic, maturation, assembly and engineering sites with an explicit
gap-aware NylC–Nyl12 map. NylC Lys189 mutagenesis supports protecting Nyl12
Lys156; exact Lys chemistry is still proposed. Mapped NylC Phe134 and Leu137
correspond to already-protected F101 and M104, respectively. Stability and
processability variants are not automatically reclassified as essential
catalytic residues. The audit found no additional implicated counterpart
outside the current 101-position mask requiring protection. NylC results do
not establish direct Nyl12 essentiality at every homologous site.


Deduplicate complete sequences, preserve every parent/sample lineage, and
assign a new campaign-specific identity. Audit the entire sequence against the
current mask, not only the requested mutation list. An earlier 100-position
diagnostic pilot allowed K156 changes; all 24 of its candidates are now
ineligible. Preserve them as diagnostics and never silently revert them.

Cross-score WT and candidates on all four TS-like parents and the matched
E–S reference with identical masks, model weights and decoding orders. One
coherent protection implementation must serve generation and scoring. Treat
likelihoods as sequence/structure diagnostics, not activation energies or
predicted catalytic rates. No structure-prediction campaign follows
automatically from a favorable likelihood alone.

## 5. Predict a small matched E–S control/design cohort

Select at most two initial design representatives, balancing original/remodeled
parentage, sequence diversity and cross-template compatibility. Predict WT and
both designs with the same two declared seeds (six models initially), using
the qualified precursor-derived paired protein MSAs and template-free Boltz-2
condition. Four PA6 L4 molecules remain unguided: no substrate contact restraints,
prescribed amide register or TS-derived protein template is supplied.

The current calibration has one successful WT seed with four relaxed positive-phi
Tier A sites; the second seed failed assembly. This establishes a usable
single-seed baseline, not reliable prediction across seeds. Retain failures and
use matched WT seeds when assessing designs. AF3 is an optional separately
qualified comparison, not a compulsory step for every candidate.

## 6. Postprocess, reject pathological assemblies, then relax

Canonicalize raw predictions; verify sequence symmetry, chain/residue identity,
all four ligands and PA6 connectivity. Reject severe nonbonded interchain
overlaps before spending time on minimization. Evaluate all 48 combinations
(four sites × four ligands × three internal amides); retain every measurement
and a one-to-one eligible assignment. Do not assume a particular amide or
residue 929 identifies the occupied or best site.

The currently matched relaxation backend is **Amber**, using the released
PA6 RESP/GAFF2 parameters with ff19SB, TIP3P, 150 mM salt plus neutralization,
and explicitly declared protonation states. The recent controls used fully
unrestrained minimization, at most 10,000 cycles and RMS target 0.05
kcal mol⁻¹ Å⁻¹. A zero exit code is not force convergence. Recanonicalize the
matching topology/restart export and verify its roundtrip before remeasurement.

An existing OpenMM adapter can consume validated Amber topologies. It is an
alternative backend, not a necessary additional relaxation. Do not switch
only designs to OpenMM while retaining Amber WT controls. Qualify the same
topology/force settings and rerun matched WT/design controls if adopting it.

## 7. Filter geometry, then validate only a few sequences dynamically

Apply explicit positive signed Herschlag phi and Tier A/B criteria to all
candidate site/ligand/amide rows before selecting assignments. Report all four
sites before and after relaxation. Native-Thr similarity remains diagnostic,
not a hard exclusion. Report scissile unit indices separately from finite-model
product lengths: a terminal internal amide of L4 gives L3+L1, not L2+L2.
An alternative E–S register cannot inherit the existing relay QM atom map.

Advance at most two or three diverse finalists to matched independent MD
replicas with WT controls. Reuse the established preparation, topology-charge,
counterion, two-replica and trajectory-analysis procedures; retain 300 K for
matched comparisons. Analyze catalytic occupancies, distributions, substrate
retention, hydrogen-bond networks, assembly and native-Thr diagnostics. Four
sites within a tetramer are correlated observations, not four independent
replicas. Geometry occupancy alone is not a catalytic-efficiency prediction.

## 8. Nominate experiments and reserve detailed chemistry for a shortlist

The PI may nominate candidates for experimental testing on sequence, assembly,
structural and dynamic evidence, with expression, maturation and activity
checks against WT. Neither a QM/MM PMF nor a verified PES TS is a prerequisite.
Detailed chemistry may precede experiments if resources permit, or follow
experimental improvement to investigate its mechanism. Experimental improvement
alone does not establish a computed barrier improvement.

For a few persistent MD candidates, optional chemistry uses MM minimization,
matched QM/MM optimization and inexpensive reaction-path tests. Regenerate
candidate-specific
atom maps, QM charges and input hashes; preserve the baseline Hamiltonian,
protonation, movable region and mechanism for comparisons. Reserve expensive
US/MBAR for the best supported candidates, using the same equilibration,
analysis, uncertainty and convergence criteria as WT. Report a mutation effect
on a barrier only with a matched WT reference and statistical support.
Post hoc WT/design protocols must remain matched
and must not be tuned to reproduce assay outcomes.

Any later SolubleMPNN layer acts on a selected catalytic parent while retaining
catalytic, maturation, binding and assembly constraints. Revalidate the combined
sequence; do not assume separately favorable layers are additive. Experimental
triage includes expression, maturation and purified WT/few-design measurements
where feasible. Interpret noisy lysate/MS assays and their limited product
standards accordingly.

## Existing implementation and evidence

| Step | Reused entry point or authoritative record |
| --- | --- |
| Persistent execution | `scripts/run_nyl12_job.py` |
| RFD3 parent qualification | `scripts/run_rfd3_parent_qualification.py`, `scripts/export_rfd3_qualified_parent.py` |
| Optional exact-symmetry experiment admission | `src/nyl12_covalent/rfd3_symmetry_acceptance.py`; [qualification contract](RFD3_SYMMETRY_QUALIFICATION.md) |
| Broad mask and sequence budget | `scripts/prepare_relay_broad_pilot.py`, `scripts/run_relay_broad_pilot.py` |
| Shared protection and audits | `src/nyl12_covalent/nyl12_sequence_protection.py`, `relay_broad_pilot.py`, `relay_ensemble_sequences.py` |
| Native sequence injection/scoring | `scripts/score_relay_ensemble_sequence.py` |
| Prediction and native model acceptance | `scripts/run_relay_broad_predictions.py`, paired-WT campaign receipts |
| Canonical export and all-amide inspection | `scripts/export_wt_control_inspection.py`, `src/nyl12_covalent/design_case_postprocessing.py` |
| Matched Amber relaxation | `src/nyl12_covalent/design_model_relaxation.py` |
| Optional OpenMM control | `src/nyl12_covalent/design_openmm_control.py` |
| MD analysis and controls | `scripts/analyze_design_md_geometry.py`, `scripts/analyze_matched_wt_occupancy.py` |
| Ensemble source/protocol | `design/nyl12_relay_ts_ensemble_20261126/` manifests and acceptance receipts |

Every stage retains runnable commands, explicit environments, input/source
hashes, immutable source snapshots, logs and process status. Scientific
acceptance is separate. Inspect failure receipts and resume the relevant stage;
never rerun successful generation to repair a scoring or publication problem.
