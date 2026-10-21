# Standalone pipeline handoff — 2026-10-20

## Latest checkpoint — 2026-10-21, before 18:00 Eastern

- User authorized work until 18:00 Eastern and a minimal AF3 water test that may
  run after disconnect: one alpha A, beta E, covalent PA6 L4 TI1 substrate I,
  and explicit HOH ligand J. Full-sequence Thr233 is E1.
- Launched `nyl12_af3_water_pa6_l4_20261021` in tmux on GPU0. Artifacts:
  `output_nyl12_af3_water_test_20261021/pa6_l4_alpha_beta/attempt_01/`.
  Check `job/status.json`, `job/job.log`, `feature_acceptance.json`, and
  `prediction/` on resume; do not infer success from session absence.
- AF3 inference also completed and wrote one predicted structure; the output
  water/covalent geometry has NOT yet been independently audited.
- Actual feature preflight PASSED: one water O and both contact restraints
  retained. Explicit HOH ligand has a non-polymer chain type. Restraints connect
  E1 N and OG1 to J1 O, range 2.6–3.2 Å, water endpoint mobile. One seed
  20261021, one sample, three recycles; guidance 0.5/interval1/20steps.
  Existing MSAs/templates/CCD reused. Source/input snapshots and hashes retained.
  Heavy-atom feature success does NOT establish NH3 protonation, water H-bond
  directions, TI1 geometry, or catalytic accuracy. No QM/MM refinement launched.
- Scripted local relay repair remains experimental and has NOT run on real
  models. Finish search-path tests/review before launch. Geometry/helper and
  CLI/export tests exist; 249 project tests passed at this checkpoint before
  final newly added tests. Full new-code coverage is not yet established.
- Canonical inspection labels are `pa6_l4` and `maxa_maxa`; original identifiers
  remain provenance aliases. Viewer source fixes `resi 113+186` and adds
  `set valence, 0`. Existing protein coordinates remain unrefined/clashy until
  a separately accepted replacement exists. Do not describe renaming as repair.
- Active journal: `output_agent_actions/20261021.jsonl`. Repeated automatic
  approval-review timeouts delayed tools; no safety rejection was inferred.
  See live journal for final export/commit/push outcomes.

## Resume checkpoint — 2026-10-21

This section supersedes historical job-state statements below.

- Both overnight neutral free-NH2 L-Thr-NME Gaussian calculations finished
  normally and passed independent rechecking: 21 atoms, CA S/CB R, 57 positive
  frequencies, minimum frequencies 45.9612 and 42.2567 cm^-1. Opt/Freq used
  wB97XD/def2TZVP; separate HF/6-31G(d) ESP grids contain 5636 and 5767 points.
  They are QM prerequisites, **not fitted force-field parameters**.
- The read-only paired audit and optimized SDF/XYZ exports are under
  `output_nyl12_resume_20261021/qm_conformer_audit_v2/`. Conformer 02 is
  4.6017 kcal/mol higher in electronic DFT energy; same-order heavy-atom RMSD
  is 1.3993 A. No thermal populations are inferred. The first audit attempt is
  retained: it failed on repeated optimization completion during Opt/Freq;
  parser correction and regression tests preceded the separate v2 attempt.
- `scripts/gaussian_log_to_xyz.py` converts logs without third-party Python
  packages. Both current logs were converted in tmux under
  `output_nyl12_resume_20261021/qm_trajectories/`: 64 and 24 Input-orientation
  snapshots, respectively, each with 21 atoms. Repeated frequency/ESP frames
  remain; these are not MD trajectories or one-frame-per-iteration promises.
- The finite-grid relay diagnostic ran in tmux under
  `output_nyl12_resume_20261021/relay_feasibility/`. Neither original model has
  an eligible oxygen site on the 5832-point grid (best maximum shortfall
  0.2510 A PA6 / 0.2679 A poly; gate 0.10 A). Optimistically excluding only
  Asn186 side-chain atoms and Thr929 N hydrogens exposes 9 / 5 eligible sites,
  respectively. No atoms were actually moved. This supports trying the
  documented rotamer/NH3 orientation repair; it does not prove a valid water
  orientation, a clash-free rotamer or an energy minimum.
- RESP still requires an explicit shared multi-conformer fit, cap constraints
  and transfer validation. Atom names repeat across THR/NME: use residue plus
  name or mapped indices. The current nonreacting beta-N charge -0.6254 is a
  charge-balance placeholder, not fitted evidence. Do not change only that
  charge or transfer neutral-template charges onto TI1 NH3+.
- User requested descriptive cache directories, a small laptop rsync package,
  a Gaussian trajectory converter, complete action logging, and commit/push.
  See `docs/laptop_inspection_20261021.md` and
  `docs/action_logging_and_rollback.md`. Logs and original scientific artifacts
  are preserved. Cache identity hashes remain metadata, not normal directory
  names. The inspection package distinguishes optimized small-molecule QM
  models from the still-unrefined protein TI1 starting models.

Next scientific work: explicit Asn186 rotamer/NH3 orientation sampling with
full-environment clash checks, followed by unchanged TI1 preflight before any
Sander launch; independently finalize RESP cap/transfer constraints. No repaired
TI1 minimum, TS1, ES model, reaction path or fitted FF is claimed by this resume.

Verification at this checkpoint: **212 tests pass**; Ruff lint, changed-file
formatting, shell syntax and diff checks pass. New measured modules have
80–97% coverage (89% combined); see
`provenance/resume_20261021_coverage.json`. Both cache directories were migrated
with all file and manifest digests unchanged; exact old/new paths are in
`provenance/cache_directory_migration_20261021.json`. Final code and documentation
are to be committed/pushed to the existing `codex/initial-import` branch at the
user's explicit request; the live action journal records the actual outcome.

## Evening stopping point — 2026-10-20, America/New_York

User requested a wrap-up for the night and authorized commit/push. Gaussian
jobs from the next section remain running independently in tmux on violet;
the laptop can disconnect. No new scientific jobs were launched in this
documentation/TS-capability assessment. Inspect current job records on resume.

- Started `manuscript/INTRODUCTION_DRAFT.md`, `METHODS_DRAFT.md`, `REFERENCES.md`
  and `SOURCE_LEDGER.md`. Existing Nyl12 manuscript sources were reviewed and
  left untouched. Completed preparation and proposed methods are separated;
  no TS, fitted parameter set or improved enzyme is claimed. Source hashes
  and revisions are in `provenance/manuscript_sources_20261020.json`.
- Read `docs/ts_search_and_symmetry_decisions_20261020.md` before TS work.
  Installed Sander supports native DL-FIND via `imin=1,ntmin=5`, with P-RFO,
  two-point numerical Hessians, Bofill updates and Hessian export. Source and
  library linkage were inspected; no runtime saddle-search test has been run.
  The interface follows the lowest mode and disables microiterations. Require
  a fresh final Hessian and endpoint tests in the declared mobile subspace.
- DL-FIND is freely available under LGPL; the user accepts it on that basis
  and requested newer alternatives. Sella is the preferred newer benchmark
  candidate (iterative Hessian diagonalization/RS-P-RFO); pysisyphus and
  geomeTRIC were also assessed. No new dependency was added or optimizer
  silently substituted. A benchmark proposal is documented, not executed.
- AF3 contact guidance already exists in
  `/data/jpc/code/alphafold3-contact-restraints` at
  `ef5dc91d65b4845f63e95ea03316f1613c6dfb5e`. Tests/source and actual historical
  production manifests/geometry summaries were located. Reuse this capability
  if needed; guided and independent prediction arms must remain distinguishable.
- Unrestrained ES ligands may occupy any equivalent site AE/BF/CG/DH.
  Score all four using a fixed register, permit symmetry-equivalent matches,
  and treat sites as nested within each model. The historical four-ligand
  analyzer's fixed I/J/K/L assignments and old atom names must not be applied
  blindly to single-ligand ES or current TI1 models. Symmetry-aware ES scoring
  remains to be implemented.

Immediate restart: check both running Gaussian jobs and acceptance records,
then continue parameter validation and TI1 geometry repair. Do not restart
duplicate calculations or claim the old water clashes have been resolved.

## Active Gaussian jobs and design clarification — 2026-10-20, evening

This section supersedes the job-state statements in earlier checkpoints.
See [method, substrate audit and multistate design plan](qm_method_and_substrate_checkpoint_20261020.md).

- User approved **wB97XD/def2TZVP** geometry/frequency calculations, with no
  additional dispersion correction, and separate HF/6-31G(d) ESP for RESP.
  Two neutral free-NH2 L-Thr-NME conformers (CA S, CB R, charge 0/singlet)
  are running on violet in tmux: `nyl12_thr_wb97xd_01_20261020` and
  `nyl12_thr_wb97xd_02_20261020`, each four cores/4 GB. Gaussian has produced
  SCF output; optimization/frequency/ESP acceptance is not yet established.
  Output root: `output_nyl12_resume_20261020/neutral_thr_nme_wb97xd/`.
  Check `jobs/conformer_*/status.json`, `calculations/conformer_*/gaussian.log`
  and eventual `process.json`/`acceptance.json`. Do not rerun blindly.
- Preparation, source snapshots, input/tool hashes, software versions, logs,
  separate scratch and content-addressed cache are implemented. Accepted QM
  data are not fitted FF parameters. RESP/cap constraints, bonded/torsional
  validation and production use remain unfinished. Existing unrelated tmux
  sessions were left untouched; no other workstation was needed.
- Completed ligand audit in tmux:
  `output_nyl12_resume_20261020/confidence_audit_v2/report.json`.
  Literal sequence is **MAXA|MAXA**, producing two MAXA fragments. Stable target
  ID remains `poly_xama_maxa`. Ligand mean pLDDT is 29.64 (59/69 atoms <50),
  terminal monads mean 13.36; reactive C/O/N 83.84/54.43/60.16. PA6 ligand mean
  60.84; reactive C/O/N 92.53/85.86/87.64. Retain full MAXA|MAXA as the
  product-consistent reference while sampling uncertain tails; no truncation
  or coordinate change occurred. Initial audit failed on rounded JSON versus
  full-precision CSV ranking scores and is preserved under `confidence_audit/`.
- Completed both descriptive metric baselines in tmux under
  `output_nyl12_resume_20261020/active_site_baseline/`:
  PA6 (d1,d2,d3,a_BD,phi) = (1.3320,2.7468,2.9584,101.9127,-99.1807);
  poly = (1.2864,2.9455,2.9865,103.3262,87.9888), distances A/angles degrees.
  Signed projected Herschlag convention matches the original reference math.
  These are measurements of original unrefined starts, not acceptance results.
- Recommended downstream plan: validated TI1 -> TS1 -> ES; focused LigandMPNN
  initially, optional RFD3 local backbone diversification; independent covalent
  AF3 and unforced noncovalent ES predictions; cheap multistate checks including
  acyl enzyme/water and TI2 geometry; reaction barriers for a shortlist.
  Map the reference cycle coarsely so TS1 is not assumed rate-limiting. This is
  a documented plan, not completed redesign or a reaction profile.
- All **149 tests pass**, Ruff passes. New-code coverage: geometry 97%, metric
  CLI 95%, template preparation 94%, Gaussian acceptance 100%, Gaussian runner
  98%, confidence audit 82%. Evidence:
  `provenance/qm_preparation_20261020_coverage.json` and
  `provenance/qm_design_checkpoint_20261020.json`. Routine scripting delegated
  to GPT-6 Luna and reviewed. No dependencies added; changes remain uncommitted.

TI1 water-contact failures remain unresolved. No new Sander optimization,
accepted TI1 minimum, transition state, ES model or validated parameter set
has been produced. Next: inspect the running Gaussian jobs without duplicating
them; fit/validate parameters after scientific acceptance; implement local
TI1 geometry repair before starting reaction paths or redesign.

## Resumed checkpoint — 2026-10-20, 17:50 America/New_York

This section supersedes the earlier pause/runner status below. The previous
bwrap-blocked resume attempt made no changes. This session also encountered
`bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted`; approved
escalated execution worked. No infrastructure settings were changed.

- Integrated the expanded audited mask into the TI1 input writer. Added
  `scripts/run_nyl12_ti1_qmmm.py`, `ti1_attempt.py`, and `ti1_checks.py` for
  explicit fresh attempts and terminal-checkpoint continuation. Legacy implicit
  build-and-run is rejected before writing. No dependencies were added.
- New attempts preserve actual inputs, code snapshots, environment/version
  metadata and hashes, original restraint references, and separate process and
  scientific acceptance records. Continuation rejects mismatched topology,
  reference, target/method, executable, DFTB data, or checkpoint digest.
- Both real-model preflights ran in tmux and have terminal wrapper exit 1:
  `/data/jpc/covalent-enzyme-design/output_nyl12_resume_20261020/jobs/preflight_pa6_l2/`
  and `jobs/preflight_poly_xama_maxa/`. Attempts are under the same output root
  at `{pa6_l2,poly_xama_maxa}/starting_model_preflight/`.
  Both passed expanded-QM and proton/connectivity inventories, then failed the
  original water-placement contact gate: maximum shortfalls 0.686619 and
  0.661263 A, respectively, against 0.10 A. These are the original LEaP seeds,
  not the improved but still failed rigid-search candidates from earlier jobs.
  Sander did not start; there is no new optimized geometry or restart.
- Routine CLI/test scripting was delegated to GPT-6 Luna and reviewed. All
  89 tests pass; Ruff lint passes. New-module line coverage: attempt 84%,
  chemistry/convergence checks 92%, CLI 96%. Full coverage evidence is
  `provenance/ti1_resume_20261020_coverage.json`. The new files are uncommitted.
- The concrete next work is local conformer/geometry repair as described in
  `docs/ti1_attempt_protocol_20261020.md`. No repair success is claimed.
  Gaussian/RESP fitting, separate neutral-NH2 ES, and parameter validation
  remain unfinished. Do not lower the water gate or repeat an unchanged job.
- A completed tmux export produced both full hydrogenated, canonical-chain
  starting models for PyMOL:
  `/data/jpc/covalent-enzyme-design/output_nyl12_resume_20261020/nyl12_pymol_starting_models.tar.gz`.
  It includes the viewer, residue maps, original audits, provenance, and hashes.
  The viewer now ends with the user-requested commands
  `show sticks, resi 113 +186 and chain A` and `util.cnc`.
  The original archive/viewer/manifest remain in `inspection_revisions/v1/`
  under that output root; the revision script and source export snapshot are
  in `inspection_source/`. No enzyme coordinates changed during this revision.

No jobs from this resumed checkpoint remain active. The two unrelated tmux
sessions listed in the historical pause note remain untouched. Code tests and
passing inventory checks do not certify physical models or MD readiness.


## Paused at the user's request — 2026-10-20

Ready to resume; no calculations from this work remain active. All tracked
job attempts have terminal exit records. Failed water searches and superseded
diagnostic attempts remain preserved alongside successful later diagnostics.
All delegated agents finished. Existing unrelated tmux sessions
`nyl12_effective_strain_20261015` and
`nyl12_focused_redesign_preflight_20261016` were left untouched.

Implementation checkpoint: `ab85811962f40f8988f9113c3c9e7b9edcc0e13c` on
`codex/initial-import`, pushed to the private remote. This pause update changes
documentation only; the last checks remain 51 passing tests and passing Ruff.
No new jobs were launched for the pause.

To resume, open `/data/jpc/covalent-enzyme-design` and read this note first.
Start with the **Next scientific work** section below. Do not rerun the old
TI1 runner unchanged: it still uses the earlier QM boundary and lacks safe
checkpoint continuation. Geometry repair and charge fitting remain unfinished.

Suggested resume request:

> Read `/data/jpc/covalent-enzyme-design/docs/standalone_handoff_20261020.md`
> and resume Nyl12 TI1 work. Delegate routine scripting to less expensive
> models, run scientific jobs in tmux, and preserve manuscript provenance.


## Authoritative locations

- New development checkout: `/data/jpc/covalent-enzyme-design`.
- Private remote: https://github.com/parksjm/covalent-enzyme-design .
- Original source and scientific artifacts: `/data/jpc/enzyme-substrate-modeling`.
- Current calculation artifacts: `/data/jpc/enzyme-substrate-modeling/output_nyl12_ti1_refinement_20261020`.
- Full scientific history: `docs/nyl12_covalent_resume_20261019.md` (historical relative output paths refer to the original workspace).

Continue code development in this standalone checkout. Use explicit absolute
input/output/tool paths while workstation-specific defaults are being removed.
Existing jobs and results were not relocated. The import namespace and CLI
names remain unchanged. The older workspace is retained for data and provenance.

## Extraction provenance

`provenance/source_inventory.json` and `provenance/import.json` record the
42-file source snapshot and archive digest. `provenance/import_snapshot.py`
verified the complete archive inventory and each file digest before writing.
The archive, original scripts, and tmux export/import records remain under
`output_nyl12_ti1_refinement_20261020/repository_export` and `jobs/repository_*`
in the original workspace. The imported inventory generator is a migration
utility with the original explicit allowlist, not a manifest of every future
file added to this repository. Git commits and the uv lock now version further
source changes.

Root README/AGENTS/pyproject/gitignore started from the recorded templates.
Post-import changes include the uv lock, pytest-cov as a development-only
dependency, and pytest's root import path for the scripts namespace. Production
dependencies remain NumPy 1.26.4, RDKit 2025.3.6, and ParmEd 4.3.0.

## Scientific status

Both TI1 topology inventories pass; neither geometry is an optimized minimum.
The rigid-pocket relay-water search has no accepted candidates. Contact
reports diagnose heavy-atom crowding as well as hydrogen contacts. Avoid
interpreting this finite rigid search as evidence against the reaction mechanism.

The expanded QM selection passed structural audits: 109/181 atoms, formal
charge 0, selected partial sum -0.0294, link-H electron counts 364/628, and a
single CA930–C930 cut. **The old QM/MM writer still needs this new selection
integrated.** This audit does not validate electrostatic embedding or DFTB3.

Gaussian 16 C.02 on violet passed a real water HF/STO-3G single-point preflight.
Use `/usr/local/g16`, `g16root=/usr/local`, and scratch
`/data/scratch/gaussian_scratch/jpc`. The current long-lived process needs
`sg g16` to acquire the already granted group; a fresh login should acquire it
normally. The preflight is not a charge-fitting result.

## Next scientific work

1. Integrate the audited QM selection and explicit new-attempt/checkpoint paths
   into the TI1 runner; require convergence and post-run chemistry checks.
2. Develop a documented local geometry-repair strategy for crowded relay-water
   sites, with subsequent contact and protonation checks. Retain failed seeds
   and distinguish restrained starting models from stationary points.
3. Prepare a documented Gaussian/RESP neutral beta-N template protocol and
   cache validated artifacts by chemical state, conformer, method, and versions.
4. Build/audit the requested neutral-NH2 ES state separately; TI1 has NH3+ at
   reacting Thr233 and an intact substrate C–N bond.
5. Complete full bonded/charge parameter validation before classical MD.
   Reaction paths and enzyme design follow reliable TI1 structures.

The measured-product gates remain PA6 L2 and poly(X,AMA) MAXA. New covalent AF3
models preserve custom atom names; no affinity calculation or naming increment
is part of this workflow.

## Skills

`scientific-tmux-jobs` is installed at
`/home/jpc/.codex-hosts/violet/skills/scientific-tmux-jobs` and versioned under
`tools/skills/`. Its installed pipeline adapter already points to this checkout.
The covalent-modeling skill is authorized but awaits a validated refinement
workflow. Consider parameter-fitting/cache, canonical inspection, and manuscript
provenance skills when repeated use justifies separate routing.

## Final import checks

51 tests pass in the standalone uv environment; Ruff passes. Overall source
line coverage is 38%, with older chemistry/topology/prediction code untested.
New workflow modules/wrappers each have at least 80% line coverage. This
coverage does not validate the physical models. The installed tmux skill now
points to this checkout. Full machine-readable coverage, source differences,
and external artifact hashes are retained under `provenance/`.
