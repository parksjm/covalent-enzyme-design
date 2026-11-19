# Manuscript source and claim ledger — 2026-10-20

## Completed constrained TI1 references — 2026-10-23

The four original optimization/frequency/conditional-curvature records under
`../output_nyl12_ti1_parameters_20261022/pa6_local_bonded/` all pass their
declared constrained-reference checks. R cap length sensitivity fails0.03A;
S length sensitivity passes. This is a failed native-transfer prerequisite,
not a failed Gaussian calculation. Data and model-scope interpretation:
`../docs/ti1_qm_results_20261023.md`.

The user's AF3-water concern and our recommendation are preserved in
`../docs/af3_water_placement_decision_20261023.md`. Earlier failed water
diagnostics do not establish a requirement for improved AF3 guidance; the
validated-TI1-MM-first route remains to be tested. No fork modification occurred.
The final strict standalone HF ESP checker has27 passing tests and100%line
coverage; full regression has447 passes/13 optional workstation-fixture skips. Prospective
surrogate-only ESP and bond-constant diagnostics do not certify a native fit.

## Ranking interpretation and symmetry verification — 2026-10-22

The user reports lower confidence and increased effective strain for many
single mutants and expects coordinated LigandMPNN designs to fold and bind
well. This is retained as a project observation/design hypothesis, not a
universal biological rule. Folding confidence, folding stability, binding
affinity and activity remain distinct. Matched no-MSA prediction is proposed,
not executed. Ranking and iteration definitions are in
`../docs/design_ranking_and_iteration_20261022.md`.

Official LigandMPNN GitHub README/run.py confirm all-chain automatic ties;
local commit `26ec57ac976ade5379920dbd43c7f97a91cf82de` and the recent
October16 protected source support explicit separate alpha/beta groups.
The existing focused-redesign workflow's section4 supplies gate/Pareto ranking
precedent. The old redesigner interface_delta helper was inspected directly:
frozen-coordinate complex-minus-protein-minus-one-ligand scoring, with an
unsafe missing-LIG fallback that must not be reused. Neither helper nor an
implemented new ranker was run. Direct LigandMPNN publisher/author-PDF retrieval
failed; primary publisher-indexed method/result excerpts were available.

## Interim reference geometry checks — 2026-10-22, 19:45 UTC

`../output_nyl12_parameter_validation_20261022/progress_audit_01/` retains
frozen Gaussian log prefixes, their hashes, the executed analysis source and
7/7/7/8-frame exports. No optimizer completed at this check. Latest frames
pass interim geometry screens; all states retain 37 atoms and 20 explicit H.
Both S cases show capped Thr N–H/oxyanion contacts of 1.77–1.93 A. S/-120
frame 4 has a failed C*–O screen (1.225 A) and is retained. No proton transfer
is detected in the latest frames. These checks do not validate the force field.

Actual PyMOL all-state coordinate/hydrogen checks and repeatable reloads are
under `../output_nyl12_parameter_validation_20261022/progress_viewer_check_03/`.
The earlier checker failure was an API mismatch for the textual `off` setting,
not missing hydrogens or a coordinate change. Its source/logs remain available.
The new optional laptop inventory has 34 focused passing tests and 97% module
coverage; the complete regression suite passed 371 tests. Versioned small
artifacts and workflow snapshots are in
`../provenance/ti1_progress_validation_20261022/`; raw logs remain on violet.

## Capped TI1 reference launch — 2026-10-22

The user requested an ACE/Thr/NME anionic bonded reference and selected
`/data/jpc/covalent-enzyme-design/cif/n12c.cif` as its torsion source. Exact
source SHA-256: `e318fddc7142be1efabd6ebbd80ad6d45bf98f6d83f1ea7c1e60add9646d7b8c`.
Author chain A Thr233 has psi 134.110276428 and chi1 -58.964417878 degrees;
native phi is undefined. Both declared cap phi values (-120/-60) and both
attacked-carbon configurations (R/S) are retained. This 37-atom, -1 singlet
surrogate supports local bonded fitting, not native NH3+/oxyanion charges.
The CIF entry label is `NYL12_CRYO-READY-FOR-MD`; no `_exptl` method record
is present. Preserve the supplied coordinates without inferring an experimental
structure determination method from the filename or format.

All four starting geometries passed the predeclared contact/connectivity
screen. Four actual 8-core/16-GB Gaussian optimizations launched in tmux;
launch is not completion or scientific acceptance. Source/CIF snapshots,
inputs, screens and process records are under
`../output_nyl12_ti1_parameters_20261022/pa6_local_bonded/`. Protocol:
`../docs/ti1_parameter_protocol_20261022.md`. Native charge fitting, conditional
curvature, held-out validation and the protein-refinement gate remain open.

An existing neutral final-HF checkpoint was converted with formchk and used
to check actual text-record parsing. The first reader failed on fixed-width
character arrays; its correction now reaches the deliberate missing-Hessian
rejection because this HF checkpoint contains no Cartesian Force Constants.
Neither probe is a validated TI1 Hessian. Retained attempts:
`../output_nyl12_parameter_validation_20261022/real_fchk_probe_01/`.

## Local water repair evidence — 2026-10-21 evening

User-directed prerequisite after the diagnostic: validate the coupled TI1
force field, then classically relax the AF3 model, then perform DFTB3/MM.
See `docs/ti1_forcefield_validation_gate_20261021.md`. Current reactive terms
remain provisional. Reactant-like movement is a hypothesis, not a validated
ES result. A trajectory-only replay reproduced the failed endpoint exactly;
all100 full-system frames are retained with source hashes.

Terminal outcome: attempt 05 completed its first QM/MM stage and failed TI1
acceptance (OG1–C54 2.195122 A; directional water bridge lost; RMS 0.63057).
Water steric screening passed with maximum shortfall 0.029899 A. Later stages
were withheld. No accepted minimum or active calculation remains.

`docs/local_repair_checkpoint_20261021.md` records the current protocol and
limitations. Actual immutable inputs, source/engine-data snapshots, commands,
atom permutations, process logs and acceptance records are under
`output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/`. Attempts 01–04 preserve
the superseded displacement, SCC failure, GB order rejection and ParmEd copy
failure. Attempt 05 uses original AF3 oxygen plus approved TIP3P water-only
preconditioning. Nonwater displacement was 0.0 A; O moved 1.408215 A through
minimization, which exhausted 200 cycles above its gradient target. DFTB3/MM
subsequently produced converged SCC evaluations. This is work in progress,
with no accepted minimum or manuscript-level mechanistic result claimed.
Final verification is `verification_08/` under the same output root: 310 tests
passed; each new scientific module has at least 84% measured line coverage.

Paths below are source records, not interchangeable datasets. File hashes and
external repository revisions are in `../provenance/manuscript_sources_20261020.json`.
The original manuscript files were read and remain unchanged.

| Key | Source | Reuse and limits |
|---|---|---|
| S1 | `/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/manuscript_modeling_section_20261011/METHODS_SECTION_DRAFT.md` | Mature assembly, product/register framing and historical modeling methods. Its Boltz, GAFF/AM1-BCC, MD and Rosetta protocols are historical and must not be silently substituted for the new covalent workflow. |
| S2 | `/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/NYL12_PA6_PA66_FINAL_INTERPRETATION_AND_NEXT_PLAN_20260928.md` | Bound-start MD and score/activity limitations; invalidated PA66 histories. No quantitative activity conclusions imported into the new Introduction. |
| S3 | `/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/NYL12_WT_MUTANT_PROTOCOL_NARRATIVE_AND_AF3_RESTRAINT_TEST_20261009.md` | Historical AF3 restraint work and interpretation. Its evolving dated sections must be read chronologically. |
| S4 | `/data/jpc/enzyme-substrate-modeling/output_nyl12_af3_contact_pa6_l4z_l2_four_site_static_20261009/` | Production RUN_NOTES, input/canonical manifests, and pre/post-relax geometry summaries corroborate 90 models, 360 sites and the implemented AF3 capability. Historical single-site atom naming differs from current TI1. |
| S5 | `../docs/standalone_handoff_20261020.md`, `../docs/ti1_attempt_protocol_20261020.md`, `../provenance/resume_20261020.json` | Current audited topology/preflight status; original water clashes still unresolved. Latest dated checkpoint supersedes prior job-state statements. |
| S6 | `../output_nyl12_resume_20261020/confidence_audit_v2/report.json` | Literal MAXA|MAXA graph/products and ligand-local confidence with original AF3 hashes. Initial failed precision-check attempt retained separately. |
| S7 | `../output_nyl12_resume_20261020/active_site_baseline/pa6_l2.json` and `poly_xama_maxa.json` | Descriptive unrefined geometry; includes hashes of current metric implementation and original projected-phi math. |
| S8 | `../output_nyl12_resume_20261020/neutral_thr_nme_wb97xd/` | Prepared Gaussian inputs, source snapshots, job status and live calculation logs. Launch is verified; acceptance remains pending. |
| S9 | `/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/sander/dlfind_module.F90`, `sander.F90`, and adjacent `libdlfind/src/dlf_formstep.F90` | Native DL-FIND capability, numerical Hessian and mode-selection evidence. Source inspection and library resolution are verified; runtime saddle-search validation is pending. |

## Newly supplied information and verification

- **2026-10-22 viewer evidence:** corrected viewers were checked twice in
  actual PyMOL 3.1.0. All 100 movie frames retain water O/H1/H2 and two O–H
  bonds; both H are visible. A comment semicolon caused the reported parser
  error. Coordinate files are unchanged. Frozen reports:
  `../provenance/pymol_movie_water_check_20261022.json` and
  `../provenance/pymol_final_water_check_20261022.json`.
- **2026-10-22 parameter inventory:**
  `../provenance/ti1_parameter_inventory_20261022.json` records the provisional
  topology digest and reactive terms. The tetrahedral center retains GAFF
  carbonyl type `c`, a 1.2183 A C54–O49 equilibrium target and several angles
  near 123 degrees; it is not a validated TI1 classical force field. These
  observations motivate reparameterization and do not establish why the
  DFTB3/MM diagnostic dissociated.
- **User clarification, 2026-10-22:** preserve a real charged PA6 terminus for
  L2-directed design. For local copolyamide chemistry, consider the shorter
  capped XA/MA model used in the Dishner study. The user supplied its explicit
  SMILES; see `../docs/substrate_termini_decision_20261022.md`. A scripted
  RDKit/tmux audit confirms neutral acetyl-X-A|M-A-N-methylamide, C29H47N5O5,
  with slash C17/O36/N18 in zero-based original input order. Protein peptide caps and substrate
  caps are distinct choices. This is a planned reference specification,
  not an accepted parameter set or demonstrated product selectivity.
- **Further user assay clarification, 2026-10-22:** PA6 L4 and L2 products were
  both detected. L4 turnover to two L2 molecules is a substrate hypothesis;
  retain the doubly charged L4 reference for that purpose. One-end recognition
  on a long polymer requires a different substrate representation. The assay
  detects only certain soluble, low-molecular-weight oligomers; distinguish
  those measured products from total cleavage and inferred endo/exo mechanism.
- **User update, 2026-10-22:** the AF3 contact-restraint fork is being improved
  independently. No new AF3 run or fork modification is part of this checkpoint.

- **User instruction, 2026-10-20:** the contact-restraint fork exists at
  `/data/jpc/code/alphafold3-contact-restraints`, has tests and successful
  production use, and remains available even though independent ES prediction
  should normally omit restraints. Verified fork HEAD:
  `ef5dc91d65b4845f63e95ea03316f1613c6dfb5e`, clean working tree.
- **Production evidence:** S4 pre-relax summary reports 90 models, 360 sites,
  250 Tier A/B sites and 35 all-four-pass models. Post-relax summary reports
  193 Tier A/B sites and 11 all-four-pass models. Both carry original AF3
  all-restraints-satisfied diagnostics for 90 models; the latter is not a
  fresh post-relax contact-restraint measurement. These are historical campaign
  outputs, not new validation of enzyme turnover.
- **User instruction, 2026-10-20:** an unrestrained substrate may choose any
  of four equivalent sites. Verified code maps AE, BF, CG and DH. Future
  independent ES scoring must be symmetry-aware; the current AE-only TI1
  metric CLI is insufficient for that purpose. Record all sites and avoid
  treating within-model sites as independent replicates.
- **Atom-name caveat:** historical restraint analysis uses C12/O4/N3/C11;
  current PA6 TI1 uses C54/O49/N59/C62. The earlier WT RUN_NOTES describes the
  graph-checked carbonyl crosswalk. Preserve an explicit mapping for every
  state; do not copy atom selectors by assumption.
- **User preference, 2026-10-20:** wB97XD/def2TZVP geometry, no extra dispersion;
  HF/6-31G(d) ESP for planned RESP. DFTB for initial paths; calibrated chemistry
  and Hessian/endpoint validation are required for TS claims.
- **User instruction, 2026-10-20:** maintain manuscript-grade records and start
  Introduction/Methods using existing Nyl12 text; compare newer saddle
  optimizers, and commit/push at the evening stopping point. This draft and
  ledger implement the writing/record-keeping portion; optimizer integration
  and TS searches remain future work.

## Publication readiness gaps

No accepted TI1/TS1/ES geometry, full reaction profile, completed RESP fit,
validated new FF, RFD3 redesign campaign or new sequence prediction campaign
can yet be described as completed. Before including historical quantitative
activity comparisons, recover their original assay definitions and raw tables.
Before release, complete the bibliography and archive accepted input/output
artifacts with a repository revision and retrievable storage location.
# 2026-10-21 execution update

The two original Gaussian runs under
`output_nyl12_resume_20261020/neutral_thr_nme_wb97xd/calculations/` completed
normally and passed geometry/frequency/ESP prerequisites. Their original
acceptance and process files remain unchanged. Independent rechecking and
optimized-coordinate exports are in
`output_nyl12_resume_20261021/qm_conformer_audit_v2/`, with source snapshots,
input/output hashes and a separate tmux process record. The initial analysis
attempt is retained as failed; no scientific result is taken from it.

The portable log converter produced 64/24 selected Input-orientation snapshots
under `output_nyl12_resume_20261021/qm_trajectories/`, with per-file provenance.
These snapshots are visualization artifacts, not a new dynamical calculation.
Relay oxygen feasibility reports are in
`output_nyl12_resume_20261021/relay_feasibility/{pa6_l2,poly_xama_maxa}/report.json`;
the optimistic atom exclusions were diagnostic only and did not alter any
protein coordinates. See the current handoff and
`provenance/resume_20261021_checkpoint.json` for consolidated evidence.


## Parallel preparation — 2026-10-22, approximately 20:00–21:00 UTC

- User inspected the interim TI1 optimization trajectories and judged them good.
  Recorded as visual review only; convergence and cap-sensitivity gates remain.
- Implemented nonlinear constrained Lagrangian curvature, optional actual fchk
  gradient parsing, analytic derivative tests and surrogate/native crosswalks.
  All 433 tests, Ruff and changed-file formatting passed in tmux; changed
  module coverage 84–96%.
  Process/coverage: output_nyl12_parameter_validation_20261022/
  parallel_final_verification_01/. No real TI1 Hessian acceptance is implied.
- Source-matched NTHR/NME/TP3 boundary proposal and graph correspondence:
  parallel_prerequisites_01/result/ plus independent native_library_extract_01/
  and native_library_verification_01/ under that same output root. Preserve
  first whole-library parser warnings and subsequent clean unit verification.
  No native coordinates, ESP fit or parameters were generated by these checks.
- Mathematical documentation consulted: SciPy NonlinearConstraint official
  API documentation, https://docs.scipy.org/doc/scipy/reference/generated/
  scipy.optimize.NonlinearConstraint.html, which documents the Hessian of the
  multiplier-weighted constraint function. This is an implementation reference,
  not a new software dependency or a molecular parameterization benchmark.


## First completed constrained TI1 references — 2026-10-22 evening

S/-60 and R/-60 completed normally and passed the declared geometry gate.
S/-60 independently re-audited and exported 37 atoms/20 H; its 1.581748 A
capped-backbone N-H/oxyanion contact remains a cap-sensitivity concern.
No real Hessian or native transfer/parameter validation is accepted yet.
Small evidence and explicit failed/repaired audit attempts are in
provenance/ti1_completed_reference_20261022/. Installed MCPB bond extraction
passed eight analytic convention checks; modified-angle and real-QM fitting
validation remain outstanding. These checks do not supply TI1 parameters.


## Standard Amber workflow clarification — 2026-10-22 evening

Installed topology.py confirms original GAFF, AM1-BCC/Gasteiger provisional
scaffold charges and a hard-coded neutral intact substrate; it is not the
new coupled-state Gaussian/RESP fit. The planned candidate is ff19SB/GAFF2
with explicit native-state RESP controls and targeted refits only if needed.
Official antechamber documentation and tutorial4b were consulted at
https://ambermd.org/antechamber/ac.html and
https://ambermd.org/tutorials/basic/tutorial4b/index.php; installed CLI was
checked for current GAFF2 support. MCPB manufactured angle checks passed
96 cases and three shared-angle examples; neither they nor the prospective
22-point local validation plan supplies or certifies TI1 parameters.


## Force-field Methods drafting and export — 2026-10-22

The user requested a force-field parameterization Methods section with references
and inclusion in the rsync pull. The integrated Methods section distinguishes
executed graph/QM preparation from pending ff19SB/GAFF2 typing, state-specific
HF ESP/RESP, targeted refitting and independent validation. It preserves the
native/surrogate charge distinction and labels water/boundary RESP constraints
as proposals. No new fitting, topology migration or model refinement occurred.

References R15–R22 cover original GAFF, Antechamber, AmberTools, the exact GAFF2
distribution, TIP3P, Gaussian C.02 and ESP-fitting lineage. The Antechamber DOI
was corrected to 10.1016/j.jmgm.2005.12.005 during primary-source checking.
The GAFF2 version comes from the installed file header, not the 2004 GAFF paper.
Gaussian's 2019 software citation was recovered from the actual C.02 output.
Publisher-access failures and unresolved full-text review are explicit in the
bibliography; no missing software citation year was guessed.

The inspection exporter selects every manuscript/*.md and docs/*.md. A fresh
package will preserve existing structures, include updated Methods/References/
Source ledger and verify every packaged digest before promotion. Exact package
job, source inventory and verification records accompany this checkpoint.


## Backbone, loop-access and tetramer design decisions — 2026-10-22

The user requested careful comparison of optional RFD3 versus LigandMPNN/AF3/
relaxation, explicit length and indel handling, gate opening and tetramer symmetry.
Working LigandMPNN symmetry examples exist and the user will point to them when
needed. The recommendation is a WT-length LigandMPNN baseline; optional later
exact-length local RFD3; indels only after explicit mapping/postprocessing tests.
Detailed contract and open/catalytic-state screening limits are in
`docs/redesign_backbone_symmetry_and_access_20261022.md`.

A scripted rigid-transform audit of the retained AF3-derived WT model found
three nearly180degree rotations about orthogonal axes (CA RMSD0.051–0.116A),
consistent with D2-like geometry. This is not a crystal biological-assembly
assignment. Preserve alpha and beta symmetry ties separately and verify the
experimental operators before configuring RFD3. Existing WT-specific exporters
contain fixed lengths/Amber indices and are not certified for indel variants.
No new backbone/sequence design, AF3 inference or physical gate sampling was
run. Source/audit records are retained at
output_nyl12_parameter_validation_20261022/repeat_symmetry_audit_01/.


### Recovered recent symmetry source and GAFF2 inventory

The user corrected the preferred source to /data/jpc/enzyme-substrate-modeling.
Recovered the October16 protected LigandMPNN preparation, explicit alpha/beta
ties (weights0.25, homo_oligomer0), protection rules and tests. Preserve its
preparation-only status and superseded-input guards. No literal earlier
point-group statement was located; the separate D2-like geometric audit is
not misattributed to that earlier prose.

Charge-deleted GAFF2 typing and parmchk2 inventory preserved the capped graph.
No ATTN placeholders, but empirical/analog estimates remain. Generic c3-os
is shared by the reactive junction and ordinary CB-OG1; targeted corrections
must avoid global collateral changes. This is not an accepted hybrid force
field. Exact successful/failed attempts and source snapshots are retained in
provenance/redesign_symmetry_and_gaff2_20261022/.

## Relay references and latest design-method review — 2026-10-24

Primary relay source summaries: ../docs/water_relay_mechanism_evidence_20261024.md.
Finite packing evidence and exact screening limits:
../docs/front_water_feasibility_20261024.md. Recent-paper comparisons and
project-specific proposals: ../docs/recent_enzyme_design_lessons_20261024.md.
Stable reference keys R23–R33 and DOI-only BibTeX records are now retained;
complete author lists remain pending editorial formatting. No new QM or design
campaign was launched by this review. Do not turn proposals into executed methods.

Retrieval distinctions: Science2025 final PMC12288761 vs preprintPMC11384011;
Nature metallohydrolase onlineDec2025 vs issue2026; Oct2026 cysteine/protease
final abstract/extended data and official code vs metalloprotease detailed
motif provenance verified only in primary preprint. Final RFD3 bibliographic
record/abstract and official tutorial verified; final full Methods/SI could
not be retrieved. RFD3 hydrogen stripping means heavy-atom donor/acceptor
conditioning is not a calculation of proton locations or reaction barriers.
Primary Nature source for Listov2025 was indexed/searchable; direct full-page
retrieval redirected to authentication. PET2025 full Methods were available.

Approval-review timeouts occurred before two root/agent read commands; no
scientific input changed. A competing-literature subagent terminated with a
model content-access error; its partial suggestions were not treated as verified
method facts. Root independently checked benign native-fold/PET primary sources.

## Proposed PA6 mechanistic project and prior MD — 2026-10-24

The user's three 100 ns WT PA6 L4 persistence report and AM1-BCC provenance are
retained in ../docs/nyl12_mechanistic_project_scope_20261024.md with the
corroborating October 9 original-project synthesis. That source distinguishes
local PA66 geometry loss from global dissociation and preserves corrected versus
invalidated histories. Exact MD topology/protonation hashes await a reuse audit.
DFTB sampling with sparse hybrid DFT PES checks is a proposal only; no new MD
or QM/MM simulation was started. No NylB pivot is authorized. Primary references
R34/R35 document existing NylB QM/MM work and approximate-QM proton-transfer
benchmarking. Source access limits are recorded in REFERENCES.md.

## Recovered WT PA6 ES data and DFTB3 installation QA — 2026-10-24

The source system is the August 18 PA6 pilot, WT Boltz model41, replicas01–03.
`output_nyl12_pa6_es_audit_20261024/reuse_audit_attempt1/reuse_audit.json`
records original file SHA256, all three completed production logs, actual
protonation/connectivity, native CIF identity checks and nine extracted seeds.
Original ligand lib/frcmod match the preparation-cache copies exactly. GAFF1,
ff19SB and TIP3P are directly supported by LEaP input. AM1BCC intermediates
support the charge workflow; exact historical invocation/version is a gap.
Original small charge residuals remain unchanged. Nine selected near-reactive
seeds are not independent population samples or optimized TI1/TS structures.

`output_nyl12_dftb_engine_checks_20261024/` retains copied shipped inputs,
references, executable/CHNO parameter hashes, process logs and numerical energy
comparisons. Water minimization and NMA/water regression series match at printed
precision. The supplied NMA regression solvent differs from the enzyme dataset;
it is an installation test, not an adopted water model or chemical benchmark.

Additional primary implementation reading: Rutgers/Amber QM/MM PMF tutorial,
https://ambertutorials-rutgerslbsr-c744272d5a9c1169e0dc9e19b8d800019105.gitlab.io/workshop/04_qmmmMTR1/QMMM_PMF.html .
It illustrates DFTB3, multidimensional sampling and unrestrained QM hydrogen
motion; its ribozyme settings are not automatically Nyl12 settings. FEBS R25
mature-substrate section was reread separately from precursor autocleavage.
Native Asp272 proximity to Thr N is measured; NylC literature correspondence
remains unverified, and no catalytic mechanism is assigned by that distance.

Full saved-frame hydration analysis completed under
output_nyl12_pa6_hydration_20261024/: three 2,000-frame records, no water O
within 4 Å of occupied-site Thr N/OG1 and no simultaneous 3.5 Å contacts.
Independent geometry/water checks of all nine seeds agree. The assay is descriptive and
conditional on the old ES/protonation model, 50 ps sampling and bound starts;
it cannot rule out a water relay or determine deacylation hydration.
Source hashes are inherited from the preceding full-file audit, with size/mtime
checks during this read-only analysis. Passed full project checks: 559 tests,
13 skipped, 95.12% coverage of all new code. Ruff lint/format passed before
coverage-import cleanup; the final formatting check is retained separately.

Reference-origin clarification (user communication, 24 October 2026):
`cif/n12c.cif` was supplied from a Nyl12 manuscript currently in revision,
with the user as a coauthor. The formal citation will be updated after
publication. No public accession, experimental method or editing history was
provided by this clarification. Describe it as the supplied structural
reference; do not infer those missing details from its filename/header.

## Catalytic-reference alignment and hASNaseIII precedent — 2026-10-26

The current solvated E-S endpoint and its preserved raw MD seed were aligned to
the supplied `cif/n12c.cif` using the N/CA/C/O atoms of Phe113, Asn186, Asp272
and Asp274; Thr233 was excluded from the fit. The raw MD seed already contained
the non-native-like Thr233 backbone/N-terminal displacement later seen after
QM/MM minimization. Reproducible reports and overlays are under
`qmmm/pa6_asp272/00_es_unrestrained/catalytic_alignment_*_to_n12c/`.

The historical WT PA6 L4-to-L2 static campaign contains exactly two
OpenMM-relaxed Tier A models, 23 and 40. Both were evaluated with the same
held-out Thr233 alignment. Model 23 retains the more reference-like Thr233
geometry and is the leading replacement seed; model 40's otherwise tighter
scaffold fit does not compensate for its 56-degree Thr233 psi difference.
Results are under `qmmm/pa6_static_es_seed/`. No QM/MM calculation has been
launched from either static model. Before use, the selected PDB still requires
reproducible topology/protonation preparation and visual/register inspection.

Reference R41 was read as a mature-substrate hydrolysis study. Its direct
N-terminal-Thr self-activation path and DFTB3-screen/DFT-refinement hierarchy
motivate a new Nyl12 hypothesis; they do not establish an Asp-independent Nyl12
mechanism, transfer numerical barriers, or describe Nyl12 autocatalytic
cleavage. Gaussian 16 is the available parallel DFT/MM backend in the present
Amber build; installed QUICK is serial. No DFT/MM job was launched during this
audit because the prior E-S seed failed the structural gate.

## Direct proton-transfer embedded-method comparison — 2026-11-10

`qmmm/benchmarks/03_direct_pt_embedded_method_comparison/` retains five exact
coordinate-drive restarts, all DFTB3/3OB, GFN2-xTB and
omegaB97X-D/def2-SVP outputs, tmux job records, source hashes, a separately
hashed parser-only retry, and the final correction table/plot. All 15 engine
calculations terminated normally. The first automatic analysis failed only
because an unprinted Sander field was requested; the retained retry used
`EPtot` without rerunning calculations.

Window-specific 20-A point-charge counts are 8,104, 8,105, 8,107, 8,095 and
8,088, which invalidates the raw cutoff-based reaction profile. Accepted claims
use only method double differences and the original periodic DFTB3/MM profile.
DFTB3 versus DFT RMSE is 1.49 kcal/mol; GFN2-xTB versus DFT RMSE is 3.91
kcal/mol. The diagnostic DFT-corrected periodic hybrid endpoint is +15.85
kcal/mol and the five selected points remain monotonic. No DFT/MM stationary
point, barrier, or free energy is claimed.

## Corrected coupled-path and fast-QM decision — 2026-11-10

The old asynchronous continuation after window 24 is scientifically excluded:
it satisfied the scalar coordinate by transferring terminal-N H1 back to OG1
rather than forming the OG1--substrate bond. The exact window-24 restart remains
valid and seeds a corrected attack-only branch under
`reaction_paths/02_asynchronous_es_to_ti1/03_corrected_attack_from_window24/`.
All three N--H identities are now explicit flat-bottom gates; the local 12 A
active region remains movable and there are no global positional restraints.

Attempt 1 is retained as a process-only failure: the minimal source snapshot
omitted `path_distances.py`, imported by the VMD exporter, and therefore exited
before Sander. Attempt 2 added the complete transitive import set, passed local
and azalea SHA-256/import checks, and launched in tmux session
`nyl12_corrected_attack_window24_attempt2_20261110`. The biased endpoint is not
called TI1. Accepted geometries will support a coupled DFTB3/PM6/PM6-DH+/
GFN2-xTB versus omegaB97X-D/def2-SVP method comparison before any new production
umbrella campaign.

Attempt 2 accepted three windows and then stopped at `attack_04`, despite native
convergence and correct HG1/H1/H2 identities, because OG1--C lagged its nominal
target by 0.115 A versus an arbitrary 0.10 A hard gate. The original acceptance
record is preserved. A separate `acceptance_revision.json` declares target lag
diagnostic, promotes the converged endpoint, and attempt 3 resumes from it in
the same calculation directory. Proton identity and convergence remain hard
gates.


## Completed water-relay diagnostic and next decision — 2026-11-17

The 25-window result is authoritative in `qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/protocol.json`, `energy_ledger.csv`, and `window_025/acceptance.json`. Final physical/bias-inclusive relative energies: +54.0235/+58.1547 kcal/mol; OG1–C 1.852 Å; both relay protons transferred, both Asp residues protonated, Thr N neutral. Final status `endpoint_negative_or_incomplete`; no unbiased endpoint release. Executed Methods and bounded Results interpretation were added to the manuscript. All 25 force-converged windows are path observations, not evidence of an index-one saddle or a stable TI1.

The controlling roadmap now proposes one 5×5 activation/attack map from the accepted E–S seed, sparse matched-embedding method checks, and candidate releases. This is planned work, not executed data. It supersedes earlier pending instructions to extend this water-relay scan. Existing method-error context: `docs/coupled_method_panel_results_20261110.md`. Mechanistic precedent remains Andjelkovic et al., R41, DOI 10.1021/acs.jcim.3c00900; the precedent does not establish mechanistic equivalence for Nyl12.

## Saved gradients and occupancy callback — 2026-11-18

Primary computational evidence is under
`qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid/`:
`reference_checks/gradient_analysis/acceptance.json` verifies 110-atom mapping
and common 104,824-charge embedding for all three saved reference states.
`results.json`, `projections.csv`, `arrays.npz` and `directional_derivatives.png`
retain the projected analytic gradients. Scripts, input/output hashes and
source snapshots are preserved there. No new QM evaluations were performed
for this analysis. Gaussian/xTB embedded-QM and native total-QM/MM derivatives
remain explicitly distinct; this limited comparison does not certify a method.

`occupancy_callback_pilot/callback_proof/scientific_acceptance.json` records
all seven actual-optimizer coupling gates passed, including expected/applied
0.20 kcal/mol bias and analytic force addition. Four native calls and two
bounded iterations were used. `job_callback_proof/status.json` records exit 0.
The separately authorized four-point diagnostic was launched only after these
gates passed; its immutable protocol and source copies are in
`occupancy_callback_pilot/pilot/`, with launcher records in `job_pilot/`.
The native NFE/PMD failure is retained and the Amber installation is unchanged.

### Numerical line-search repair — 2026-11-18

Retained-trial evidence and reproducible scripts are in the grid's
`occupancy_callback_pilot/line_search_diagnosis/`: `audit.json`, `steps.csv`,
`cutoff_pairs.json`, and `pair_energy.json`. Ordinary MM pair 11272/11330
crosses the 10 Å cutoff; pair-removal energy and force predictions explain the
observed jump within native output precision. No new QM evaluation was used
for this diagnosis. The installed `dlf_lbfgs.F90` implementation is linked
through `src/nyl12_covalent/dlfind_lbfgs.py`; library/source hashes are retained
in `pilot/continuation_001/protocol.json` and the new point process records.
Original stop records are preserved in named history folders. A passed the
unchanged force gate after three further optimizer steps; B then started.
The wrapper's real-library and input/resource tests passed (12 tests, 95%
statement coverage); together with analytic-coordinate/bias checks, 43 focused
tests passed. Continuation-script line coverage was not measured. No native
installation, dependency, chemical state, Hamiltonian or force criterion changed.

### Complete occupancy pilot and endpoint-D tests — 2026-11-18

`occupancy_callback_pilot/endpoint_audit/` retains exact endpoint/input
coordinate and hash checks, native-energy/bias checks, all 54 proton owners,
heavy-bond retention and carbonyl/amide geometry. The common E–S energy
(−450493.0781 kcal/mol) was checked against its actual native mdout.
`audit.json`, `summary.csv`, `supplement.json` and `sourcehashes.json` are the
primary evidence for the four-point table and oxyanion-contact observations.
No additional QM evaluation was used for that audit.

The one-candidate zero-bias release is in `pilot/point_D/unbiased_release/`,
with entry point `release_endpoint.py` and `job_D_release/` launch records.
The separate exact-D electronic check is in `d_reference/`, with immutable
`preflight.json`, source copies and `job_D_reference/` launcher records.
Its nonperiodic embedded-QM scope is explicitly distinct from native periodic
total QM/MM. No TI1, new barrier, or method acceptance is asserted at launch.

### Literature-grounded mechanistic reframing — 2026-11-18

The terminal D release/reference files were reviewed: the release force-
converged to reactant geometry, and the fixed-embedding DFT/xTB gradients both
resisted additional attack in the tested direction. Local primary-paper
extraction and source hashes are in
`analysis/nyl12_mechanism_evidence_20261118/local_evidence.md` and
`local_provenance.json`; reproducible acid-reach coordinates and hashes are in
`acid_reach.json`/`extract_acid_reach.py`. No new electronic calculation was
performed for this review.

New references R44–R48 were checked against primary publisher/author/PubMed
records. R48's indexed author-manuscript passages on pp. 793 and 795 explicitly
state reversal of early proton transfer after release and lack of metastable
intermediates within the calculated FAAH acylation/deacylation path. Its direct
PDF requests returned 429/redirect errors, so indexed primary text and the
PubMed abstract were used. R46 is a caution: its concerted deprotonation/attack
still forms a tetrahedral intermediate. R44/R45 support scissile-NH/nitrogen-
inversion hypotheses and experimental amidase-specificity engineering, not
a universal fixed reaction sequence. R47 is the published 2025 Chem Catalysis
version of the earlier homolog-screening preprint. Primary mature NylC Fig. 10
(R25) is explicitly BapA-guided and does not demonstrate a TI1 lifetime.

The roadmap now prioritizes complete proton delivery and C–N cleavage rather
than a required stable TI1. This is an evidence-based prospective revision;
no replacement pathway campaign or enzyme mutation was launched. The actual
proton donor, resting microstate and existence of a TI1/acylenzyme basin in
Nyl12 remain unresolved. Low polymer conversion is not a measurement of kcat
or of a unique chemical free-energy barrier.

### Existing-MD leaving-group relay screen — 2026-11-18

Primary evidence is
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/00_screen/`:
`bounded_screen_acceptance.json`, `bounded_screen_counts.csv`,
`primary_attempt1/screen.json`, `controls_attempt1/screen.json`, and
`diagnostics_attempt1/{diagnostics.json,primary_nearpasses.csv}`. Each screen
attempt retains input/source/trajectory SHA-256 records and source snapshots;
`primary_job/`, `controls_job/` and `diagnostics_job/` retain exact commands,
environments, logs and successful process statuses. The primary null
`selected.json` was preserved after the controls stage; no microstate was
promoted or restart selected.

All 60,000 source frames were represented in the cached geometry tables;
18,958 passed unchanged Tier A/B plus positive φ, 51 passed all structural
gates, and none had the Nalpha–water-O–Nleave heavy bridge needed for either
hydrogen-bond class. Native OG1/Nalpha displacement was the main rejection
source. Exact thresholds and six-case counts are recorded in the executed
Methods section and campaign README. Nine focused tests passed, with measured
coverage 94% for the screen, 99% for reusable geometry and 96% for diagnostics
(96% combined); Ruff passed. The negative is bounded to sampled mature-
substrate E–S conformations; it neither excludes transient/multiple-water
reorganization nor assigns an acid state or catalytic free-energy barrier.

## pKa consolidation — 2026-11-19

Corrected mature NylC PROPKA uses
`analysis/negoro_2023_protonation/crystal_propka_20261116/attempt3/results/asp306_asp308_pka.csv`;
precursor rows only use the attempt2 counterpart. Mature N-terminal values
are audited in `terminal_pka_20261116/before_after.json` under the same
analysis root. Nyl12 uses
`qmmm/pa6_flipped_model30/09_propka_static_tier_a/screen_attempt1/catalytic_pka.csv`.
The student `evidence/pka/source_mapping.json` records exact source hashes
and filtered precursor rows. Manuscript pKa Methods/Results retain corrected
values, unmodified precursor-mutation caveats and explicit MCCE rejection
status. No new pKa calculation was performed for this consolidation.
