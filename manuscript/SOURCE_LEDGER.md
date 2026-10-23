# Manuscript source and claim ledger — 2026-10-20

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
