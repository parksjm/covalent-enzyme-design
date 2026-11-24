## CURRENT — three converged releases; window033 curvature running, 2026-11-23

- **GPU follow-up now active:** two independent 100 ns replicas each of WT,
  broader mask_07, focused mask_03 and V111C, four PA6 L4 ligands per tetramer.
  Violet runs WT/broader first; Azalea runs focused/V111C first; every sequence
  has its second replica queued on the other host. All starting sites are
  flipped (negative signed projected phi). WT seed04 has four Tier A sites
  but zero native-Thr passes; do not discard that known-active control or tune
  the gate. Follow `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/README.md`.
- **Mandatory Azalea runtime:** read `AMBER_WORKSTATION_RUNTIME.md` and root
  `AGENTS.md`. The generic pmemd26 binary lacks Blackwell kernels; the already
  validated CUDA12.8 build is required. This turn's two failed generic starts
  are preserved in remote campaign clutter; corrected jobs advance native MD.
  The reproducible host launcher now rejects the wrong binary/hash/configs.
- Canonical masks, all substitutions and the testable flipped-pose hypothesis
  are in `ligandmpnn_mask_site_rationale_20261123.md`, manuscript Methods/Results
  and slides36–38 of the updated43-slide deck. Current WT negative-phi Nalpha–
  leaving-N distances are longer on average than positive-phi ones; sign alone
  is not evidence for a competent proton relay. No new positive-pose MD launched.
- All three native optimizations force-converged. Window033 retained distorted
  tetrahedral-like geometry (C–OG1 1.740 Å, C–N 1.490 Å, C–O 1.309 Å;
  HG1–Nalpha 1.050 Å). Window032 is noncovalent but N-terminal-protonated;
  window028 reverted to E–S-like geometry with HG1 on OG1.
- Window033 repeated-gradient and parity checks passed. Coarse .001/.0005 Å
  curvature probes failed; an additional72-evaluation fine-step diagnostic
  found a stable plateau without changing thresholds or the Hamiltonian.
- Full12495-DOF minimum Hessian now runs in tmux
  `nyl12_relay_w033_full_hessian_20261123`, h=.000125 Å, eight nice workers,
  unchanged12 Å movable region and fixed exterior. Follow
  `minimum_verification/06_full_hessian/progress.json` and
  `minimum_verification/verification_summary.json` under window033.
  **Minimum order is pending**, not established by force convergence.
- Canonical VMD exports now include the exact final accepted restart and
  the established QM/substrate display conventions; final-coordinate checks pass.
  No new reaction-path or umbrella campaign was launched.

## Earlier checkpoint — calibrated gate and relay releases, 2026-11-23

- Unchanged native-Thr gate: all ten original 9Z0A/9Z0B sites and four
  derived n12c sites pass; crystal OG1/N maxima 0.1232/0.0841 Å. Known-active
  rejection still prevents activity classification. Full calibration:
  `native_thr_gate_calibration_20261123.md`.
- Three unbiased DFTB3/MM releases run in Violet tmux:
  `nyl12_relay_{ti1_window032,ti1_window033,early_window028}_20261123`.
  Existing direct-path `umbrella_sampling/nterminal_relay_releases_20261123/`
  holds exact frames, manifests and verified progress VMD bundles. Initial
  SCC/energy/gradient checks pass; final convergence is not yet claimed.
  Original 12 Å/4165 movable atoms, fixed exterior, charge −2 and atom map
  are retained, without positional or reaction-coordinate potentials.
- 35.81 kcal/mol is existing direct NEB **ΔE**, not ΔG‡. Its TS candidate
  and 47.05-kcal/mol apparent NEB-US PMF already exist; no duplicate campaign
  was launched. Existing first-order TS verification remains open.
- NylC MCCE 5XYG/7YU0/7YU1 reruns completed natively but remain pKa-
  unqualified; no states assigned. Audit curves and templates before use.
- Progress PowerPoint/PDF: `docs/presentations/nyl12_progress_20261123/`.
  The pull script includes it, calibration and relay VMD. Muse means
  Glimmer 30B; no new model server launched.
- Follow `NYL12_MECHANISM_ROADMAP.md`. A TI1 minimum is optional; proton
  occupancy and connectivity, not labels alone, define chemical states.
- PI offline instruction: leave all releases to their convergence checks.
  `nyl12_relay_minimum_w033_20261123` is actually running and waiting after
  all 40 prerequisite hashes passed. It tests window033 for a minimum only after
  independent force convergence and tetrahedral-like retention. Repeat/parity
  and derivative-step checks precede a full 12495-DOF Hessian with expected
  index zero; keep the original 12 Å movable shell and fixed exterior. No new
  reaction path is authorized automatically; an inconclusive curvature test
  remains unverified. See `relay_minimum_queue_20261123.md` and the roadmap.

## Earlier checkpoint — restored pKa controls and TS gate correction, 2026-11-23

- NylC5XYG/7YU0/7YU1 reruns passed native MCCE step1/2 and are in DelPhi
  step3 on Violet. Each case has `rerun_20261123/` and a current pointer in
  the existing `repaired_comparison_20261120/` campaign. Accepted heavy
  coordinates are unchanged; the REMARK parser failure is fixed.
- The old TS electronic gate falsely rejected recovered SCC inner trials.
  New parsing retains retry evidence and rejects actual fatal/unrecovered
  events and temperature changes. Fresh native gradients accept both
  saved endpoints as E-S/acylenzyme. TS repeats are identical and meet the
  original1e-4 Eh/bohr force criterion. Candidate-to-E-S energy is46.41
  kcal/mol; it remains a high-energy route, not a certified viable barrier.
- First-order TS/IRC remain unverified because of Hessian soft-mode error.
  Next recommended step is12 targeted HB-column evaluations at the same
  candidate; no full-Hessian redo, optimizer restart or shell reduction is queued.
- Use `NYL12_MECHANISM_ROADMAP.md` for the concise forward queue and
  `nyl12_stationary_validation_loop_20261121.md` for detailed diagnosis.
  The five-point evidence is at `.../validation_loop_20261121/diagnostic_repeat_20261123/`.
  Mechanism-agent next has a bounded human-reviewed NylC pKa QA task; servers
  stay off between requests. Original failures/inputs/acceptance are retained.

## CURRENT — missing laptop movie resolved, 2026-11-22

PI uses `vmd -e view_umbrella.vmd`; the laptop lacked `path.dcd`. The pull
script's generic `--max-size=100m` silently omitted the 453,002,956-byte
(433 MiB) movie. The NEB-US bundle now uses a dedicated 500 MiB cap and verifies
`vmd_dynamic/MANIFEST.sha256` after transfer. A real rsync regression with a
file exceeding 100 MiB passes. Copy the updated
`/data/jpc/covalent-enzyme-design/scripts/pull_nyl12_inspection.sh` to the laptop
and rerun it, then launch VMD from the movie directory. No additional path
fallback or GUI chooser was introduced. Full trajectory/topology are retained
on Violet; all previous viewer fixes still apply.

## CURRENT — matched-mask structural campaign finished, 2026-11-22

All 52 new Boltz predictions and Amber minimizations finished; including three
qualified reused WT controls, all 55 endpoints are force-converged. Geometry
tables are complete. The last plotting step rejected new cohort names;
reporting-only recovery is in `ensemble_mask_comparison/validation/reporting_repair/`,
with original failure receipts preserved and no science rerun.

Joint Tier A/B + native-Thr passes: broader 12/80 sites in 7/20 model seeds,
focused 8/80 in 4/20, WT 1/20 in 1/5. Restricting to negative signed Herschlag
phi gives broader 11/80, focused 3/80, WT 0/20. Broader mask_07 has 6/20 joint
sites in 3/5 seeds, all negative phi. Active controls K59D and V111C have
0/20 and 7/20 joint sites, respectively, so this is not a validated activity
classifier. The WT baseline now uses matched Violet runtime for all five
seeds; two previous mixed-runtime seeds were regenerated. Counts are
structural observations, not evidence of lower barriers or higher activity.
See `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/README.md` and
`validation/analysis/comparison_summary.json` for source-bound tables.

## CURRENT — umbrella VMD display repair, 2026-11-22

The NEB-US `vmd_dynamic/view_umbrella.vmd` was regenerated in place after the
PI reported a blank VMD display. Fixed invalid three-argument VDW syntax,
explicitly activated frame zero before labels/camera measurements, and wrapped
loading in a fail-fast procedure. Original viewer/checksums are retained in
`viewer_repair/original`; the repair script/source and receipt are alongside.
The trajectory and topology hashes are unchanged (1,890 frames, 19,967 atoms).
The existing pull script includes the repaired bundle. Tcl loader regression
and checksum tests pass; GUI VMD is unavailable here, so laptop display remains
to be confirmed. Scientific jobs were not changed.

## EARLIER CHECKPOINT — structural pilot and mask comparison launch, 2026-11-22

Live matched comparison: `nyl12_matched_masks_validation_20261122` on Violet.
Native sequence sampling completed (48 audited, 39 unique; broader24/24unique,
focused15/24). Eight representatives and K59D/V111C controls are selected.
`ensemble_mask_comparison/validation/status.json` reports the actual phase: 52
new Boltz predictions, three eligible reused WT models, then automatic Amber
relaxation and four-site analysis. Input/protection/hash/native-runtime checks
passed; postprocessing/startup failures are retained separately.

The 45-model Boltz/Amber campaign finished: all 45 predictions and minimization
processes completed; **44/45** reached the force tolerance. Variant 05 seed 04
has RMS gradient 0.1455 kcal/mol/Å and remains diagnostic. After relaxation,
WT has 3/20 joint Tier A/B+native-Thr sites in 2/5 models. Variants 03/04 each
have 7/20 in 4/5 models; variant 06 has 8/20 in 3/5 models. Seven of eight
variants have a joint pass. Across all endpoints there are 40/180 joint sites;
restricting to force-converged cases gives **38/176 sites in 20 models**.
Raw models have no Tier A/B sites. Full counts, source hashes and model paths
are in `design/nyl12_ligandmpnn_pilot/boltz_validation/analysis/comparison_summary.*`.
These five-seed structural observations are not catalytic-activity estimates.

PI authorized a matched 19-versus-29-position mask test on two accepted WT
backbones, 48 raw LigandMPNN samples, unchanged catalytic/maturation protections
and explicit A–D/E–H sequence ties. New mutable positions passed SASA and
interface checks. At most four distinct sequences per arm plus K59D/V111C
controls proceed through five Boltz seeds and matched Amber relaxation. New
work is on Violet; WT reuse requires runtime/input identity. See
`design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/README.md` and its
actual job receipts. No TS geometry or activity-trained filter is introduced.

Future laptop pulls exclude Hessian column/displacement trees. The original
workstation data and already downloaded laptop copies are preserved.

## CURRENT — NEB-seeded US closeout, 2026-11-22

Only the historical neutral-Thr, ASP272/ASP274 direct-route NEB-seeded US at
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/`
is closed out here. All 63 windows completed six accepted 5 ps blocks
(30 ps/window): eight Violet workers exited 0 on November 22 at 10:43–10:47
UTC and 24 Azalea workers exited 0 on November 21 at 11:55–20:49 UTC. Only the
local analysis PID 4100242 and sync PID 4098079 were stopped with
SIGTERM, both absent afterward; the receipt is at
`umbrella_sampling/dispatch/stop_20261122/violet_monitor_stop_receipt.json`;
no campaign Sander remained active on either host.

The apparent PMF maximum is 47.0475 kcal/mol versus 46.9414 in the same E-S
density bin. Connected overlap and target ESS 7.29 do not establish
convergence; windows 028/030/032 have drift flags. This is not a validated
barrier. Stop additional investment in this forced direct route, preserve raw
data/restarts/sources and the negative finding, and queue no extension. Other
jobs, design work and separate ASH274 NEBs remain outside this decision.

The verified dynamic movie is at
`umbrella_sampling/vmd_dynamic/view_umbrella.vmd`: 30 frames/window,
1,890 frames and 19,967 atoms. All selected frames match original reaction
distances within 0.00005 Å; bundle checksums pass. The laptop pull includes
the final bundle. Independent biased windows are concatenated for inspection.

## CURRENT — 45-model structural validation snapshot, 2026-11-22

PI approved WT plus eight unique LigandMPNN variants, four PA6 L4 molecules per
tetramer and five matched seeds per sequence. The latest campaign `status.json`
snapshot is `running`: **21/45 normalized; 15/45 relaxed**. These are receipt
counts, not a fresh host-wide liveness audit. See
`design/nyl12_ligandmpnn_pilot/boltz_validation/README.md` and
`design/nyl12_ligandmpnn_pilot/boltz_validation/verification/checkpoint_receipt.json`
for detailed records.

The first WT seed completed MM minimization with force convergence (final RMS
gradient 0.048767; target 0.05 kcal/mol/Å). Its single-seed geometry preflight
reports raw 0/4 Tier A/B, 4/4 native-pass and 0/4 joint-pass; relaxed 1 B, 2/4
native-pass and 0/4 joint-pass. Four sites are nested observations in
one model; these results do not rank variants or imply activity. The preflight
tables and plots are under `boltz_validation/checks/WT_seed_01_analysis/`.

The existing modular RESP/GAFF2/ff19SB/TIP3P preparation uses neutral Thr NH2
and anionic ASP272/ASP274, with no restraints during minimization. Prediction
hosts use the same frozen Boltz source and checkpoint SHA256
`090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1`; Python,
Torch and RDKit builds differ and are recorded in
`boltz_validation/dispatch/runtime_versions.json`. PI authorized Tulip for
Boltz and Geranium for Amber relaxation after the campaign availability check;
both were confirmed same-architecture as Violet. The tools were already
available, so no dependency installation was needed.

The campaign README records the coordinator session
`nyl12_design_validation_checked_20261122`, Violet prediction sessions
`nyl12_design_boltz_violet_gpu0_20261122` and
`nyl12_design_boltz_violet_gpu1_20261122`, Tulip session
`nyl12_design_boltz_tulip_20261122`, and Geranium relaxation sessions
`nyl12_design_relax_checked_gpu0_20261122` and
`nyl12_design_relax_checked_gpu1_20261122`. These are recorded session names,
not a host-wide liveness report. Preserve per-model receipts and failed
preparation history; do not infer completion or silently rerun failures.

Parallel PI request: recompute MBAR with the latest NEB-US sampling and show the
PMF. That scientific analysis is handled separately; its current process state
is not audited in this structural checkpoint. No new QM/MM or mechanism branch
is authorized by this structure-validation task.

## CURRENT — Model23 LigandMPNN pilot complete, 2026-11-21

The complete-runtime run finished with exit status 0: 24 raw samples, eight unique
sequences and 16 duplicate samples, with 4-6 substitutions per precursor. All
chain-length, symmetry, fixed/protected-identity and tensor/FASTA checks passed.
Native runtime was approximately 4 minutes 14 seconds for the four six-sample jobs. The tmux
session ended normally. Read `design/nyl12_ligandmpnn_pilot/README.md` for the
eight variant mutation sets and `accepted_parent/run/results/` for FASTAs and
full lineage. No candidate structure prediction, relaxation, MD or activity
validation was launched. Existing QM/MM, US and WT MD remain untouched.

The PI accepted the recommended 24-sample single-parent Model23 pilot and
instructed work to continue while signing off. Only violet and azalea may be
used. Preserve all existing QM/MM, US and WT MD jobs and the queued matched-WT
analysis; do not begin new mechanism branches or model-server benchmarks.

Campaign: `design/nyl12_ligandmpnn_pilot/`; accepted inputs and live statuses:
`accepted_parent/run/preparation_manifest.json`, `accepted_parent/run/job/status.json`,
`accepted_parent/run/results/status.json` and campaign `status.json`. Read those
receipts for launch/completion; absence of a tmux session alone is not success.

Original WT40/41 failed occupied-site native-Thr checks (OG1 shifts1.261/1.480 A,
N shifts2.149/2.530 A, psi deltas-56.03/-38.17 degrees). Only Model23 passed the
seven-model audit: OG1 shift0.188 A, N0.705 A, delta-psi-30.09 degrees and
chi1+4.52 degrees. Cutoffs remain unchanged. The accepted source identity and
coordinates match the physically split input exactly. Model23 is static Tier A:
d1/d2/d3=2.885/2.886/2.873 A; prior Model23 QM/MM deterioration is retained as a
negative result. One-parent sampling does not test opposite-pose diversity.

Native parsing passed all1312 residues, alpha A-D lengths232 and beta E-H96,
33 ligand heavy atoms C24N4O5, deterministic A-H ordering and complete backbone
coordinates. The runtime is existing LigandMPNN Python3.11/Torch2.2.1+cu121,
NumPy1.23.5/ProDy2.4.1, GPU hidden and one CPU thread, low priority on violet.
No dependencies were added. Four jobs use temperatures0.1/0.2 (dimensionless),
two distinct RNG seeds per temperature and six samples/job. The mask has17
alpha positions plus beta39/41, explicitly tied across copies with weight0.25.
Catalytic/maturation protections are unchanged; lengths stay232/96, no indels.

After completion inspect `accepted_parent/run/results/unique_sequences.fasta`,
`samples.csv`, `report.json`, `scientific_acceptance.json` and independent
`tensor_fasta_acceptance.json`. Duplicate sequences are reported, not replaced.
Candidate folding, geometry after relaxation, binding and activity remain
unvalidated. No candidate prediction/relaxation/MD campaign is authorized here.
The no-argument laptop pull includes compact pilot artifacts in the same layout.

The immutable old pair preparation and failed delegated setup attempts remain
under `prepared/` and `backbone_audit/`. Authoritative native measurements are
from the separately frozen lead audit under `native_preflight/`; the CIF has313
observed author-numbered residues, not a distinct313-residue construct. Missing
coordinates are positions1-7,228-232,326-328. Other copies' native poses are
diagnostics; all-copy identity/overlap integrity is required.

The initial native launch failed before sampling because the frozen source lacked
LigandMPNN's bundled OpenFold package. All failed receipts and source remain
under `accepted_parent/job/` and `results/`; no sequence was generated there.
Installed support files were copied without adding dependencies. Full native
imports then passed. The active complete-runtime attempt is under
`accepted_parent/run/`, tmux
`nyl12_model23_ligandmpnn_24sample_complete_runtime_20261121`.

## CURRENT — first filter audit complete; noisy assay evidence retained, 2026-11-21

PI authorized WT/known-active calibration and planning a small symmetry-aware
LigandMPNN campaign while all existing QM/MM/US jobs continue. The controlling
master plan v0.3 now permits exploratory E-S-based sequence engineering before
verified TSs; TS-specific claims still need verified stationary points and TI1
only if present. No new design generation or prediction/MD campaign launched.

`analysis/nyl12_filter_calibration/` completed in tmux: 900 prescribed static
models (WT+eight mutants,50 each,PA6/PA66 separate); historical counts reproduced,
zero Tier A/B gate mismatches. Every variant retains an A/B model; PA6 WT7/50,
V111C1/50,K59D18/50. Angle failures dominate. This does not validate activity
ranking or the newer MD filters. Thresholds unchanged. Fourteen focused analysis
tests passed (92% combined coverage; each module>80%), two laptop tests passed.

PI emphasizes crude-lysate/MS noise and discordant **K59E** experiments.
The primary cb150 XLSX is in the older repo's `config/` (Git-ignored):
`cb150_nyl12_ssl_sites_1_4_pa6_pa66_pmaxa.xlsx`.
Its 99 cells are retained with sheet/cell identity in `experimental_workbook_audit/`.
PA6 L2 B5/B6=0/1.1465832774; legacy0.5732916387 is their mean, not a selected
experiment. Do not average silently or use that mean as clean ground truth.
WT is absent from the workbook and supplied by downstream normalization.
General manuscript methods describe independent re-screen n=3 and triplicate MS
with standards, but these cannot yet be mapped to the two K59E source rows.
The uploaded DOCX text and latest Figure7 image are retained in
`manuscript_context/`. PI believes image/workbook use the same data, with
experimentalist confirmation deferred a few days. No stronger linkage assumed.

Plan: `docs/nyl12_ligandmpnn_pilot.md` and
`configs/nyl12_ligandmpnn_pilot_plan.json`:24 raw samples across two opposite-phi
candidate WT backbones40/41, pending identity/native-Thr/assembly preflight.
Explicit alpha ABCD / beta EFGH ties; protect catalyst/maturation, no indels.
Noisy folds are not candidate ranking targets. The native-backbone theozyme
idea and NylB/poly(MAXA)->MAXA are documented future hypotheses, not new queues.
Mechanism-agent remains human-reviewed; current local models are not qualified
for autonomous scientific planning. No new local-model deployment was started.
The normal pull script now includes compact calibration artifacts and current
planning docs, using only violet. Source/workbook/failed cosmetic checks retained.

## Current job inventory — 2026-11-21 evening

Read `docs/CURRENT_JOBS.md` and `provenance/current_jobs.json` for the latest
live/queued audit. Ten of twelve matched WT replicas are accepted; the last two
are running. Eight violet US workers and five violet NEB engines are active.
Azalea has model servers but no native scientific calculation running.
TS directional follow-up completed but failed electronic acceptance; no verified
TS. All three newer NylC MCCE comparison controllers stopped on an unsupported
PDB REMARK record; they are not still running or queued. Earlier dated entries
below retain the launch history and must not be treated as current liveness.

## CURRENT — automatic final matched WT analysis queued, 2026-11-21

PI requested that remaining WT MD finish, then compare Tier A/B occupancy across
all three matched states and screen the newly completed trajectories for relays.
Tmux: `nyl12_final_matched_wt_analysis_20261121`.
Directory: `qmmm/pa6_flipped_model30/10_wt_protonation_md/analysis_after_md/`.
The live receipt gate accepted 8/12 at launch and waits for all twelve 100 ns
replicas, including the two currently running and two still queued. It validates
both production stages, successful terminal job, explicit 100 ns marker and
topology/source hashes. Upstream failure stops analysis without changing MD.
See `waiting/status.json` and `job/{status.json,job.log}`. It polls every 60 s
with a 48-hour waiting budget. All MD/US jobs are left running unchanged.

After completion, all-frame occupancy/native analysis and four new relay screens
run on two low-priority single-thread CPU processes. Eight completed matched
relay results are reused, with retained original topology paths and explicit
byte-identical canonical aliases. Fresh occupancy trajectory hashes must match
cached relay hashes. Legacy unmatched cohorts are excluded. State labels include
the Thr N-terminal charge so NH2 and NH3+ cannot be merged accidentally.

Compare Tier A, exclusive B, inclusive A-or-B, signed phi, native-Thr gates and
their intersections across four independent replicas/state, plus 10 ns stability
blocks. Relay counts are per unique sampled frame, joined to all-frame native
metrics. Planned products: `plots/tier_occupancy.png`, `tables/`,
`occupancy/summary.json`, `relay/comparison/summary.json` and candidate records.
No new dynamics, protonation, minimization or reaction-path launch is queued.
Frozen source/config/environment retained; root focused suite: 40 passed.

Additional hosts: geranium and tulip are PI-authorized. The initial default-key
failure was resolved using the user-specified `~/.ssh/id_ed25519_aster`; this key
is authorized for all workstations. Explicit identity SSH as `jpc` succeeded on
both. Geranium: 48 logical CPUs; tulip: 24. Both: about 503 GiB RAM and two
idle RTX 6000 Ada GPUs (49,140 MiB each at inspection). Geranium Amber paths:
`/data/jpc/code/ambertools26/venv/bin/sander` and
`/data/jpc/code/pmemd26/bin/pmemd.cuda_SPFP`; not on the noninteractive PATH.
Tulip's Amber installation was not located in the bounded paths checked.
See `provenance/workstation_resources_20261121.json`; no key material was read.
No jobs/transfers/installations were started. Existing MD/US and analysis waiter
are unchanged; verify native environment/compatibility before future launches.
PI restriction: use geranium and tulip simultaneously only for jobs expected
to finish within a few hours. Before each launch, verify nobody else is using
the proposed host, including sessions, CPU/GPU workloads and queued work.
Ask the PI when resource availability or expected duration is uncertain.

## CURRENT — relay screen/refinement complete; diagnostic candidates selected, 2026-11-21

PI approved the expert review's next step: screen existing WT MD for a complete
proton-delivery network to leaving amide N, then select a bounded candidate set.
Campaign: `qmmm/pa6_flipped_model30/11_leaving_group_relay_screen/`.
Fourteen accepted local 100 ns replicas cover historical NH2/ASP/ASP, ASH272,
ASH274 and completed newer matched NH2/ASP/ASP, NH2/ASP/ASH274 and NH3/ASP/ASP.
Two completed original-pose WT replicas add the opposite signed-phi orientation.
All sixteen coarse cases completed: 1.6 μs, 32,000 screened frames at 50 ps.
Actively written trajectories are excluded; historical charged replicas only on
azalea are explicitly outside this local first pass. Original source MD retained.

The original two low-priority workers and the added original-pose queue exited
zero. Baseline refinement/selection exited zero, but its `selected/` files
mostly identify known zero-intermediate direct-route controls. They are not
newly discovered assisted relays. All raw records are retained.

The authoritative assisted-only queue completed exit zero in
`nyl12_assisted_relay_candidates_20261121`. It requires ≥1 water/internal
residue before ranking, evaluates native geometry, and refines ±4 saved frames.
Final report: `assisted_review_refined/{selection.json,summary.md,candidate_readiness.json}`.
One unique strict Tier B frame is matched NH2/ASP/ASP replica 03, frame 9386:
Lys156→water29282→water4517→leaving-N. OG1/Nalpha native shifts 2.396/2.958 Å
fail current gates. Best single-water nearpass is ASH272 replica 1, frame 7316:
OG1→water9397→leaving-N, Tier A; terminal H–N 2.50964 Å misses the 2.5 Å
screen by 0.00964 Å, but Thr displacement/psi gates also fail.

Three exact source restarts were exported as **diagnostics**, not accepted new
QM/MM seeds. None of 20 bounded native-audit shortlist frames passed every native
gate. "Joint qualified" in automated selection refers only to HB+Tier; readiness
adds the remaining criteria. Coarse/refined counts repeat source frames; frame 9386
is one unique candidate. All original records remain; no claim that unexamined
MD frames or reorganized networks cannot work. No new MM/QM/MM or altered
protonation was launched. Inspect these candidate structures/native criteria
before deciding the next chemical test. US inputs and jobs remain unchanged.

Real-frame pilot passed. Root ran 29 focused screen/geometry tests, the
selector's seven tests after correcting a ranking-order bug, and four assisted
review tests. New-code coverage
is above 90% in delegated receipts. Original source snapshot retained; corrected
selection snapshot is `selection_revision_02/source_snapshot/`. The original
orchestration parent was intentionally paused before its unexecuted selection
stage. The replacement waited for both successful worker exits, then retired
only that parent and completed baseline refinement/selection.
Read `jobs/postprocessing_revision_{guard,outcome}.json`; expected parent143
does not indicate scientific failure. Current US jobs/inputs are unchanged.

Read `docs/NYL12_MECHANISM_ROADMAP.md` for the controlling protocol. A geometric
relay is not a proton-transfer barrier or validated mechanism. Inspect actual
per-job status and scientific_acceptance.json before claiming completion.

## CURRENT — preliminary NEB-seeded PMF; sampling continues, 2026-11-21

Updated requested refresh frozen at 19:08 UTC: all 63 windows, 1,575 ps total,
15–30 ps/window (8×15, 8×20, 23×25, 24×30 ps). Native two-CV MBAR gives
q_acyl marginal maximum **46.94136 kcal/mol** above the unchanged E–S density
reference bin; fixed-peak-bin 95% block-bootstrap interval 45.24714–47.92707.
All 50 bootstrap refits connected; target ESS 6.88, drift flags 028/030/032.
Mixed proton branches and sampling uncertainty remain; no activation-barrier
or sampling-convergence claim. Earlier snapshot below is retained for comparison.
Source/data freeze: `analysis/pmf_refresh_20261121_190812/` within the US campaign;
`analysis/pmf_latest/` now exports its verified compact results for laptop pulling.
No MD/US settings or jobs changed. Fifty focused bias/MBAR tests passed;
source/input identities, original reference bin and independent 2D marginal
verified. Initial publication gate rejected regenerated copied Python bytecode;
the retained correction excludes derived caches only, preserving all source/data
checks. No estimator or sampling data were altered for that correction.

PI requested continuation, documentation and a scoped Git checkpoint.
See `docs/neb_us_preliminary_findings_20261121.md` and versioned plot/table
evidence. Frozen MBAR: 63 connected windows, 1,190 ps, target ESS 28.28,
q_acyl marginal maximum 46.49584 kcal/mol above the fixed E–S density bin
(pointwise 95% block-bootstrap interval 44.53939–47.55736 at that bin).
Preliminary mixed-proton-branch projection; not a converged activation barrier.
Source regular NEB: 35.80895 kcal/mol at image 8, endpoint +4.92846,
native convergence cycle 477; separate TS verification still required.

HG1-on-Thr-N occurs in windows 028–034 while both imposed CVs closely track
targets; the distance differences do not uniquely specify proton ownership.
Do not merely strengthen k, silently filter branches, or change active inputs.
All samples retained; initial 13.96 supplemental plot is explicitly INVALID
from weighted-histogram cancellation, with original artifacts preserved.
Use `umbrella_sampling/analysis/pmf_latest/primary_q_acyl_pmf.png`.
The laptop pull now includes that compact export.

At ~14:32 UTC, 8 violet and 23 azalea workers were active, with one remote
worker finished; completed coverage was 8×10 ps, 54×20 ps and 1×30 ps.
Keep the unchanged campaign running toward its minimum 30 ps/window.
Current controllers: `dispatch/coverage_first_20261120/{violet,azalea}/`.
No simulation was stopped, restarted or reconfigured for this checkpoint.

## CURRENT — Hessian complete; semantic model comparison reviewed, 2026-11-21

Full 12,495-DOF Hessian completed. Stored symmetric matrix: one negative and
12,494 positive eigenvalues; chemical proton-transfer/acylation mode about
1307i cm^-1. NOT YET a fully verified TS: global heuristic antisymmetry gate
left 4,535 soft-mode signs unresolved and stopped downstream descents. Audit
localizes most asymmetry to a small coordinate subset; do not repeat the full
matrix or lower the gate simply to claim success. Read
`docs/nyl12_stationary_validation_loop_20261121.md` and retained
`validation_loop_20261121/completed_hessian_audit_20261121.{json,md}`.
Targeted native follow-up is actually running in tmux
`nyl12_candidate_directional_followup_20261121`: 27 bounded native probes
with two nice workers and conditional independent candidate-mode downhill
minimizations. Stable pointer: validation_loop_20261121/diagnostic_current.json.
Original acceptance remains unchanged; these descents are diagnostics, not a
claim of a verified TS or formal IRC. No full matrix redo was launched.

Both model inference runs completed 14 cases. Actual native responses were
semantically reviewed independently of format compliance. Rubric-based
lead-agent grades: GPT-OSS 2 critical-error cases and 3 major-error cases;
Mistral 1 critical-error case and 6 major-error cases. Mistral corrected the
acylenzyme/deacylation trap but both failed H+ charge conservation. This small
single-sample pilot does not show convincing overall Mistral superiority.
Strict response-contract compliance was GPT-OSS 2/14 and Mistral 0/14;
format failures were not counted as chemistry errors. Neither is qualified
for autonomous scientific planning/dispatch. Native Mistral serving did work;
old controller_status.json belongs to a failed earlier attempt. Use current
phase and per-model completed status plus the paired semantic review in
private `benchmarks/mistral_small4_20261120/LEAD_SEMANTIC_REVIEW_20261121.{json,md}`.
Original responses, source freezes and scientific acceptance remain intact.

## CURRENT — model comparison, stationary-point validation and substrate proxy

PI clarification: PA6 L4 is a tractable local proxy for a likely much longer
polymer. L4-to-two-L2 remains a model reaction, not an established experimental
precursor/product pair. Binding direction and proximity of either physical
terminus are unknown. Both-capped interior versus near-N/near-C one-ended
models are future explicit sensitivity tests; no current chemistry, charges,
parameters or frozen model-benchmark evidence were changed.

The official Mistral Small 4 FP8 setup/benchmark runs independently in tmux:
private `/data/jpc/mechanism-agent/benchmarks/mistral_small4_20261120/`.
Read `current_phase.json`, `prepare_controller_status.json`, `controller_status.json`
and `BENCHMARK.md`. Frozen 14 cases comprise 6 development and 8 held-out-in-workflow
cases with auditor-only rubrics. GPT-OSS completed 14: 12 response-contract failures
and 2 contract-valid responses. Lead review of those 2 found critical acylation/
deacylation and proton/charge/topology errors; the other 12 are chemically
ungraded. No overall chemical-error rate or Mistral superiority is claimed.
Official Mistral files downloaded; use current receipts to establish whether
hash/preflight/server/inference stages have completed. Existing science preserved.

Nyl12 stationary validation remains UNVERIFIED. Stable pointer:
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/ts_image08_attempt1/validation_loop_20261121/current.json`.
An isolated DL-FIND precision-only output patch preserves installed engines.
Energy parity/old rounding bounds passed. Tight SCC/density convergence did NOT
remove coarse-column discrepancies; the nonmonotonic step sequence does not
prove a clean truncation explanation. Exact requested-to-native evaluated-point
identity and numeric warning-token checks are part of the current gate repair.
No complete Hessian or first-order index has been claimed. Follow the exact
pointer/retained failed attempts; a session name or successful optimizer exit
is not TS verification. Keep all 4,165 movable atoms and original 0 K Hamiltonian.

Public strategy/resource/calibration docs pushed through 621dcb6. Private model
benchmark/tool/strategy software pushed through 7090d58 (read Git for latest).
Automatic review blocked a broader public scientific/manuscript payload;
those files remain local. Tests: main 1,753 passed / 13 skipped; private 148 passed;
benchmark 28 passed / 92% coverage; stationary focused coverage reported separately
in its receipts. Development-only Matplotlib/MDAnalysis declarations fixed
pre-existing full-suite collection errors; production pins unchanged.

Additional verified live checkpoint before PI logout:

- Nyl12 full 12,495-DOF Hessian is actually running with 8 workers in
  `nyl12_stationary_verified_frame_attempt05_20261121`. The stable current.json
  above points to the authoritative source-frozen attempt. Actual native
  coordinate/force identity, old/new serialization parity and a fine-step
  derivative plateau at 0.00025 Å passed. SCC tightening alone did not
  remove coarse-step differences. Full inertia and both-sign downhill/minimum
  Hessians remain conditional; no first-order saddle claim yet. New 24 focused
  tests passed at 95% coverage. All older owned attempts are terminal/retained.
- Model benchmark continuation is `mechanism_mistral_benchmark_attempt04_20261120`
  on violet. HF-format expert tensors proved incompatible with the installed
  loader. Official SAME-revision native FP8 consolidated artifacts are now
  downloaded separately by the detached chain, then hash/key validation,
  native-format preflight/server and unchanged 14 Mistral cases run automatically.
  Mistral has not yet been benchmarked at this checkpoint. Private continuation
  software pushed 439b025; 29 benchmark tests passed. Read current_phase.json.
- Existing capped/register-model source map is
  `provenance/existing_capped_pa6_models_20261121.json` plus readableCSV.
  It indexes 15 representative rows / 8 campaigns, with 220 retained relaxationdirs
  observed in that bounded scope. Both neutral capped PA6 trimer registers have
  BCC-labeled GAFF1 MOL2/LIB/FRCMOD/PRMTOP assets, correcting the earlier missing-
  parameter statement. Native AM1 execution is unverified from conflicting
  retained SQM files; small ±0.004 e topology drift is recorded, not silently
  normalized or blamed for mechanism failure. The two graphs are identical;
  register labels do not establish binding direction. No qualified capped RESP
  release is assumed.
- A fresh one-component AM1-BCC/GAFF2 baseline completed in
  `capped_pa6_am1bcc_gaff2_20261121`, under
  `parameters/capped_pa6_am1bcc_gaff2_baseline_20261121/job/status.json`.
  Native SQM confirms explicit AM1 optimization, charge 0 / singlet. All 17
  isolated-ligand technical gates passed; native charge residual -0.001999 e
  is retained within documented three-decimal charge-output precision. Actual
  script coverage is 99%. The shared component is ready for explicit enzyme
  preparation/atom mapping, not certified enzyme MD. No enzyme MD or capped
  QM/MM system was launched; all original assets remain unchanged.

Existing MD/US/MCCE jobs were not stopped or reconfigured. Earlier checkpoints
below are history; their current status must come from native per-job receipts.

## CURRENT — three independent NylC preparation/MCCE queues active, 2026-11-20

Authoritative campaign:
`analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/`.
Mature5XYG queue waits exact running Amber19SB/GBn2 repair. Precursor7YU0/7YU1
queues wait verified PDBFixer builds, then run the same full-Amber refinement
and strict sequence/bond/clash/true-terminus gates before native MCCE. Residue267
stays internal in both precursors; all deposited mutations retained. No fresh
MCCE pKas are claimed before the engines and interpretation gates complete.

Current sessions: `nylc_mcce_5xyg_20261121`,
`nylc_mcce_7yu0_attempt03_20261121`, `nylc_mcce_7yu1_attempt03_20261121`.
Earlier precursor queue failures/status-schema fixes are preserved. The first
builder source was formatted after launch without a retained pre-format copy;
that source-history gap was disclosed, and authoritative builds were replayed
once from immutable source snapshots (job_verified/build_verified). Do not use
the initial diagnostic-only candidates. Native MCCE source/binary/parameter
archive and all profiles/helper sources are retained.

The two Violet GPU jobs are WT Nyl12 protonation MD, not PA6 FF followups:
nh3_asp272_asp274 replica01 and nh2_asp272_asp274 replica02 under
qmmm/pa6_flipped_model30/10_wt_protonation_md. No jobs were stopped or moved.
Gpt-oss120b's scientific review failures block autonomous chemistry/job dispatch;
local bounded extraction/reporting remains useful with validation. No model
replacement/download or new permission was enabled.

## CURRENT — independent mature/precursor NylC controls requested, 2026-11-20

PI explicitly requested repair/rerun5XYG plus uncleaved7YU0/7YU1 now. Native
soft-PDBFixer candidates failed peptide/ring geometry. Technical diagnosis
found no atom-copy bug; weak soft geometry is insufficient here. Full
Amber19SB/GBn2 new-atom/H refinement of5XYG is running (observed heavy atoms
fixed; temporary building H stripped before predictors). Both precursor
PDBFixer builds are launched independently, retaining their deposited mutations,
fullA/B1–355 chains and internal267. Their follow-on full-Amber refinement and
MCCE queues are under repaired_comparison_20261120 in the original NylC MCCE
analysis campaign. All geometric/terminal gates remain; results diagnostic,
not accepted protonation assignments. Read per-job status before claiming
any MCCE engine started or completed.

The two Violet pmemd.cuda_SPFP jobs are authorized Nyl12 WTprotonationMD,
not PA6forcefield followups: nh3_asp272_asp274/replica_01 and
nh2_asp272_asp274/replica_02 under qmmm/pa6_flipped_model30/10_wt_protonation_md.
No processes were stopped or moved. gpt-oss120b failed critical scientific
review criteria; retain bounded/read-only utility roles, no autonomous
mechanism selection/job submission. Original failed reviews remain evidence.

## CURRENT — 5XYG geometry continuation running, 2026-11-20

### Controlled reconstruction continuation launched

At 21:41 EST Nov20, `nylc_5xyg_softmin_continue_20261121` started the
new-atoms-only refinement from the retained failed candidate. Work directory:
`analysis/negoro_2023_protonation/pdbfixer_repair_20261120/repair_run_20261121/`.
Read `softmin_continuation.log` and `softmin_continuation_status.json`.
The separate `_softmin.pdb` requires its validation receipt before use.
No valid repaired structure or new NylC pKa result is claimed at launch.

## CURRENT — 5XYG candidate failed geometry gate; controlled refinement authorized, 2026-11-20

PDBFixer native repair wrote a candidate at
`analysis/negoro_2023_protonation/pdbfixer_repair_20261120/repair_run_20261121/result/5XYG_repaired_protein.pdb`
but exited 1 at 21:35 EST: modeled alpha-A tail peptide A9–A10 C–N=1.967 Å.
No scientific acceptance was issued. Original source and failed candidate are
preserved. The PI-authorized repair continues with one explicitly documented
same-PDBFixer-soft-potential minimization of newly modeled atoms only, freezing
all original heavy coordinates, tighter force tolerance and unchanged bond/
terminal/cleavage gates. Check the continuation's native status, not this
candidate filename, before using coordinates for any pKa calculation.
No MD/QM or US job was altered. Detailed trust/nearby-pKa findings and
manuscript notes are in `docs/pka_trust_and_active_site_audit_20261120.md`.

## CURRENT — pKa trust audit; PDBFixer repair authorized, 2026-11-20

The PI requested investigation of suspicious MCCE values and then explicitly
authorized PDBFixer installation/5XYG repair. Read
`docs/pka_trust_and_active_site_audit_20261120.md`. Nyl12 native processing and
structural gates pass, but quantitative pKas are not qualified: saturated
curves, a poor Asp272-F fit and an inconsistent His257-E >14 fit label; the
latter conflicts with positive occupancy <=0.04 throughout pH0–14. The pinned
FTPL zero PRO-backbone charges differ from bundled legacy parameters; source
provenance alone does not validate them. No parameters or MD/QM states changed.

Matched HIS-normalized PROPKA predicts Lys156 8.46–9.91, His257 4.85–5.22,
Tyr267 14.00–14.57 and Asp45 -1.97 to -0.06. These are follow-up flags,
not catalytic assignments. Actual MD site pocket RMSDs ~1 A mean equivalent
chains are not identical static structures. Source-hashed comparison/nearby
evidence and scripted curves are in the MCCE pilot's titration_audit_20261120.

Historical mature5XYG MCCE was rejected; uncleaved7YU0/7YU1 were staged but
never run with MCCE. PROPKA was run for all, without backmutation. PDBFixer
1.12.0/OpenMM8.6.1 are installed in tools/pka-preparation using an isolated uv
lock. Repair staging/status is analysis/negoro_2023_protonation/
pdbfixer_repair_20261120/. Split mature alpha1–266 and beta267–355 before gap
reconstruction; preserve source atoms, record modeled residues and gate all
true termini. Verify native repair/acceptance before claiming repaired input.

## CURRENT — MCCE pilot complete; NEB-US coverage still incomplete, 2026-11-20

Verified at about 21:10 EST. Corrected MCCE/DelPhi pilot completed at 21:02 EST
(exit 0): all 2,527 native conformer entries have completed electrostatics and
pair-energy files; no failed-conformer report, native steps 3/4 complete and
all four reported artifact hashes match. Results are preliminary, not a
validated microstate assignment. Read the completion section of
`docs/mcce_preparation_repair_20261120.md` and canonical `mcce_pka_summary.json`.
Asp272 fits are <0 except chain F 1.326; Asp274 <0 in all four copies;
Thr233 N-terminal fit 6.513 in E and <0 in F/G/H. Strong chain variation,
protein-only preparation and limited sampling/dielectric/ensemble qualification
preclude adopting those numbers as definitive catalytic pKas.

NEB-US: 32/63 windows each have 20 ps accepted production (640 ps pooled).
The 31 remaining windows are advancing under the authorized coverage-first
schedule (8 violet, 23 azalea); the remaining remote worker waits at the
all-63-first-block gate. Latest analysis update_00042 contains those same
32 production windows and two disconnected overlap components. No global PMF
or barrier is supported yet. Earlier all-window first-5ps estimate remains
approximately 03:00 EST Nov21, contingent on throughput; coverage alone will
not demonstrate convergence. No jobs or scientific inputs changed in this audit.

## CURRENT — optional design and interaction resources recorded, 2026-11-20

The PI proposed RFD4-Proteina over RFD3 and tmol, PLIP/ProLIF, the existing
PA6 parameter release and AF3 contact-restraint fork as expert resources.
They are options, not a mandatory stack. Read
`docs/redesign_tool_inventory_20261120.md` and the private mechanism-agent
`docs/TOOL_CATALOG.md`. Foundry/RFD3 checkpoints and AF3 checkouts already
exist; RFD4/tmol/PLIP/ProLIF were not found in bounded environment checks.
The local AF3 checkout differs from current GitHub source. PA6 L4 release
scope and compact-conformer limitations are recorded; do not invent catalytic
state parameters. No packages, models or jobs were installed/launched/changed.
The fixed-length LigandMPNN baseline and current mechanism priorities remain.

## CURRENT — MCCE structural repair validated; complementary pilot active, 2026-11-20

The PI requested an honest mechanism-agent tool list and MCCE preparation
repair. Private catalog: `/data/jpc/mechanism-agent/docs/TOOL_CATALOG.md`.
Implemented agent operations remain evidence/review and bounded PROPKA;
general MD/QM/tool adapters are planned, not yet enabled.

MCCE failure was diagnosed and fixed at the input level: the old writer
right-justified atom names, while MCCE matches exact four-column names.
The earlier global 3 Å explanation was wrong: ddvv is squared distance, so
its unchanged cutoff 3 Å² is about 1.732 Å. No MCCE source/cutoff patch was made.
The corrected complete Nyl12 input retains all 9,928 heavy atoms and 8 OXT atoms.
Native step1 now passes 8 NTR + 8 CTR, zero error messages and 0.000 Å displacement.
Both new preparation/audit modules pass 19 focused tests (86–87% coverage), Ruff.
Read `docs/mcce_preparation_repair_20261120.md` before resuming.

Canonical experiment directory:
`analysis/nyl12_mcce_step1_preflight_20261116/canonical_atom_names_repair_20261120/`.
`scientific_acceptance.json` is true. Conformer generation completed. All 227
native zero-charge notices were matched to explicit official FTPL zeros;
audit retained, no charges changed. One native-generated zero ammonium torsion
remains an explicit model assumption. Initial PBE launch had a PATH error,
now archived; no method change or conformer rerun. Corrected continuation:
`nyl12_mcce_pbe_titration_20261120`, `pbe_retry_job/status.json` and
`pilot_status.json`. Four lower-priority one-thread DelPhi workers, dielectric
4/80, salt 0.15 M, pH 0–14, six native MC runs at 298.15 K, fixed seed 20261120.
Outputs and stage budgets are scripted. Bounds remain bounds, and completion
will not automatically assign protonation or validate microscopic pKas.
All US and unrelated jobs remain unchanged. Email authorization is pending.

## CURRENT — mechanism-agent pKa capability verified, 2026-11-20

The PI requested pKa tools. Private `/data/jpc/mechanism-agent` now has an
explicit typed local PROPKA 3.5.1 adapter and optional review enrichment using
registered structure IDs or a checked existing report. Read `docs/PKA_TOOLS.md`
there. No new dependencies, model-issued shell, automatic protonation changes,
MD/QM dispatch or new model inference were enabled. Full suite: 148 passed,
95% coverage, Ruff lint/format passed. Native default/coupled Nyl12 validation
completed in tmux; 12 target sites/mode, all 16 physical-chain termini and
28 HIS groups passed. The resulting three-role evidence packet validates.
Accepted artifacts: `runs/pka_tools_validation_20261120/attempt2_nyl12/`;
source/commands/failure history: `runtime/pka_tools_validation_20261120/`;
current attempt is identified by `current.json`.

The historical Nyl12 predictor PDB used Amber HIE names and had no HIS summary
groups. A separate explicit HIE-to-HIS copy changed names on 280 atom records
(28 residues), retaining coordinates and all originals. Do not silently replace
historical results with this different preparation. A stricter native trial
rejected the previous mature5XYG NylC input: C-266 groups on alpha chains A/B
were missing from the prediction, with coordinate-start N+ groups at 19A/3B.
Its original output/registration failure is retained, and that input is disabled
pending terminal/missing-tail preparation validation. MCCE remains blocked by
its previously failed preparation. None of these failures establishes a
protonation state. Private code/data have not been published to the public repo.

US timing update: the slowest local workers project first 5 ps production
across all 63 windows around 02:51 EST Nov 21; communicated approximately
03:00 EST (provisional 02:30–03:30). Remote handoff succeeded: 21 azalea workers
had started queued-window heating. Exact native-step-based timing evidence is
in `dispatch/coverage_first_20261120/first_block_eta_20261120.json` under the
canonical US campaign. A PMF still requires connected overlap and convergence
assessment. Native sampling and sync/MBAR/gate monitors continue after logout.

## CURRENT — NEB-US coverage-first scheduling active, 2026-11-20

The PI authorized finishing active 5 ps blocks, then prioritizing all missing
windows before further extensions, preserving existing data. At 17:18 EST,
all 32 original Sander PID/start identities were retained and integrating:
8 violet and 24 azalea. Each was in production block 04; 32 windows already had
15 ps accepted, 31 were queued (including gap windows 045–052).
The canonical campaign remains the historical 35.80895 kcal/mol NEB's
`umbrella_sampling/`; no new pathway or scientific inputs were created.

Active control is `dispatch/coverage_first_20261120/plan.json` with matching
source snapshot on both hosts. Native Sander was never signaled. Only old
Python queue controllers were retired and window drivers suspended; new
tmux workers wait for exact native block completion, record exit status via
Linux pidfds/procfs, preserve original process metadata, validate the native
checkpoint, and run queued windows through their first 5 ps production block.
An independent tmux coordinator `nyl12_us_coverage_first_gate_20261120`
releases extensions only after all 63 windows have accepted diagnostic-ready
first blocks. Each worker then extends its least-sampled assigned window to
the original minimum 30 ps. The original source snapshots/plans/receipts and
native outputs remain. Retired controller exit 143 is intentional; consult
new per-host `worker*/status.json`, postflight and activation records.

All 39 focused tests pass, with 86% combined new runtime statement coverage;
Ruff lint/format pass. Existing dispatcher/sync/analysis checks passed; one
plotting regression required the installed analysis environment because uv
lacks matplotlib and passed there. A dummy Linux child/parent handoff and
real accepted-checkpoint-prefix validation on both hosts passed. Global
sampled overlap remains disconnected; a preliminary PMF is still gated by
overlap, not merely window count. Existing sync/MBAR/pull monitors continue.

## CURRENT — mechanism-agent email delivery blocked by ORNL relay, 2026-11-20

The authorized setup test was accepted locally but bounced upstream. The PI's
Postfix journal trace identifies queue `3241E43E00D2`: `relay.ornl.gov` rejected
parksjm@ornl.gov with `554 5.7.1 Relay access denied` at 16:35:21 EST. Nothing
relevant arrived in Junk/quarantine. End-to-end notifications are unavailable.
ORNL IT must authorize violet's relay access or provide an approved submission
route. The exact relay policy reason remains unknown; the undefined HOSTNAME
warning is not established as its cause. No system configuration changed and
no duplicate test was sent. Private setup, delivery-trace evidence and current
status are in `/data/jpc/mechanism-agent/docs/EMAIL_NOTIFICATIONS.md` and
`provenance/email_delivery_failure_20261120.json`; original receipts are intact.
Notification software remains tested (101 tests); no scientific jobs changed.
Email replies are not automatically processed. This local resume entry and
private notification evidence have not been published to the public repository.

## CURRENT — local pathway review completed; scientific corrections required, 2026-11-20

The PI requested mechanism-agent pathway ranking and one next-calculation
proposal using the corrected experimental baseline, in tmux. The independent
local-model review completed all five stages in307.3seconds under
`mechanism_pathway_review_20261120_extended`.
Authoritative private project: `/data/jpc/mechanism-agent`; current output:
`runs/pathway_review_20261120_extended_budget/`; process metadata:
`runtime/violet/pathway_review_extended_budget_job/status.json`.
Read its `docs/ROADMAP.md` and `docs/STANDALONE_OPERATION.md` before resuming.
No GPT API is used by local inference; chat/laptop may disconnect. The final
files are `PATHWAY_REVIEW.md` and `PI_DECISION.json`. Read the separate
`review_results/pathway_review_20261120/LEAD_SCIENTIFIC_AUDIT.md` before using
the model recommendation: critical nucleophile/protonation/source/energetics
errors remain and the proposed TS job was not launched. All science proposals
await PI review; no model-directed scientific tools exist.

Three earlier incomplete/invalid requests are retained with exact identities;
request schemas now bind stages and source IDs, require mechanistic comparisons,
and include explicit inference budgets. Full custom suite66passed, with review
modules99-100% coverage; source extraction tests91%. No native science changed.
The user-supplied current commercial-polyamides DOCX and papers are indexed in
the private project. Three original PI intake questions are answered; a new
manuscript caption/Methods assay inconsistency is pending PI clarification.
This local resume entry has not been published to the public repository.

## CURRENT — read-only mechanism-agent pilot completed, 2026-11-20

The PI approved the user-owned azalea Ollama runtime, official gpt-oss:120b
(~65 GB), official ollama/Pydantic uv dependencies, and one idle GPU with two
CPU inference threads. Runtime and full single-GPU residency are verified.
The five-case read-only pilot completed; all five outputs passed structural
checks, but scientific review remains pending. Unsupported model claims were
identified; no scientific execution permission was granted. It has no scientific execution tools. This supersedes the
older documentation-only status below. The private repository is
`/data/jpc/mechanism-agent`; its current roadmap, source/input hashes and
retained setup failures are authoritative for that pilot. The next PI task
is review of the generated five-case report packet in the private repository.

US remains 8 violet + 24 azalea native jobs. Eight of 63 windows have completed
production blocks (now two blocks each); overlap remains disconnected, so no
full-path PMF is defensible yet. Production-based timing evidence is in the
canonical `umbrella_sampling/timing_update_20261120.json`: first expanded
assessment around 12:30 EST today; all-window first-block coverage around
10:05 EST Nov 21. These are throughput extrapolations. A preliminary PMF still
requires connected overlap, and converged barriers need further diagnostics.
All native scientific jobs and persistent analysis/sync monitors continue.

## CURRENT — approved two-host US running, 2026-11-20

This supersedes the older pending-transfer checkpoint. The PI explicitly
approved transfer, native compatibility testing and up to 24 azalea CPU
workers. There are now **32 actual native Sander processes**: 8 on violet,
24 on azalea; positive native steps were verified for every remote window.
Eight old local queue-controller PIDs alone were retired after exact argv/
PID-start-time checks. All eight native violet jobs remained unchanged.

The same canonical US directory exists on both hosts beneath the historical
35.80895 kcal/mol NEB. `status.json` points to current execution records.
`dispatch/plan.json` assigns 16 windows to violet and 47 to azalea, each once.
`distributed_launch_acceptance_20261120.json` records verified native PIDs,
progress, prepared hashes and the 100-step native probe. The probe initial
EPtot agreed exactly at printed precision; tolerance was 0.001 kcal/mol.
All 1,046 copied inputs/source files were hash-verified. Azalea uses installed
native Sander/Python; no binaries or environments were copied. `pmemd.cuda`
is for classical MD and was not involved in this CPU QM/MM launch.

A remote tmux/libtinfo conflict initially prevented all 24 window processes
from starting. Those failed controller/window-job artifacts remain under
`dispatch/preparation_history/ncurses_attempt1/` on azalea. The tested fix
isolates only tmux's LD_LIBRARY_PATH, retaining inherited/explicit scientific
library paths in job scripts and metadata. Remote retry source/plan lives in
`dispatch/azalea_retry02/`; canonical `dispatch/azalea` points to its workers.
The original plan/local controllers and all frozen scientific inputs remain
unchanged. No numerical MD samples were lost to that pre-MD failure.

Persistent monitors: `nyl12_neb_us_result_sync_20261120` mirrors compact,
validated completed-block outputs and endpoint views onto violet every 60 s;
`nyl12_distributed_neb_us_analysis_20261120` updates CPU-only MBAR/plots.
Both are configured for seven days. No raw trajectory/restart transfer is
required for monitoring/laptop inspection. Source/transfer receipts are kept.
Analysis source remains the qualified CPU/overlap-guard snapshot.

At the latest completed-block refresh, 8/63 windows had one 5 ps production
block each. Sampled overlap is disconnected, so global PMF generation is
correctly withheld. `first_block_summary_20261120.json` records competing
HG1 transfer to Thr Nalpha in middle windows; all valid samples are retained.
The current `analysis/latest/` contains coverage and overlap plots, not a
converged barrier. The remaining azalea windows are heating. Initial measured
remote speed ~2.9 ps/hour is encouraging but does not establish production
throughput. Provisional minimum-campaign ETA is around two days; revise it
from complete production blocks and extend if convergence requires it.

The updated pull script includes current diagnostics/portable views with
workstation-identical layout and includes small active frozen Python snapshots and excludes nested staging,
compiled caches, job duplicates and historical attempts. Recopy the script to the laptop.
The user requested roadmap+pull+PMF status+mechanism-agent status+commit/push
before logging off. Native jobs/monitors continue independently of chat.
Mechanism-agent remains a private documentation scaffold, with enzyme-general
profiles/NylC/NylB comparisons and five-case read-only model qualification;
no local LLM runtime/model service was installed or dispatched.

Verification: 1,643 collectable tests passed and 13 existing tests skipped.
Full collection remains blocked by missing MDAnalysis in the available test
environment; standard uv also lacks matplotlib. No dependencies were added.
Focused launcher/dispatch/sync regressions and actual-rsync pull test pass.
The native architecture probe and actual 32-process progress checks are
separate from software tests and sampling convergence.

## CURRENT — US capacity review and explicit transfer approval pending, 2026-11-19

Eight violet US windows remain in heating, with no production block complete.
Initial measured speed is ~1.32 ps/hour/window: at least 30 ps production plus
6 ps preparation across 63 windows extrapolates to ~8.9 days with eight workers.
This is heating-only, contention-dependent timing, not a reliable completion
promise. Full evidence is in the canonical US `timing_estimate.json`.
Violet has24 physical cores and is busy; azalea has28 physical cores nearly idle,
working installed Sander/analysis Python and two 96 GB GPUs. Its pmemd26 binary
matches violet's hash; no GPU compatibility simulation was run.

A tested static dispatcher is prepared in `scripts/dispatch_neb_umbrella.py`:
eight controllers adopt the currently running local window jobs and each then
run one new local window;24 remote controllers cover the other 47 windows.
Every window has exactly one owner. This could reduce wall time to about two
to three days if azalea throughput is comparable; benchmark first. It does
not change Hamiltonian, inputs, seeds or minimum sampling. It has NOT launched.
Stop only the verified old worker Python controllers during takeover; native
window jobs/drivers must stay alive. Frozen sources and explicit host hashes
are required.19 focused dispatcher/jobs tests pass; core coverage97.33%.

Automatic review rejected the azalea staging transfer of topology, 24 seed
restarts and campaign inputs/scripts, explicitly requiring confirmation.
No transfer/probe/remote US launch occurred. An explicit async approval asks
the PI to approve that transfer, a short native compatibility check and up to
24 remote CPU workers. Do not retry the rejected transfer without that reply.
The original 8 jobs continue. A scalar/observable/results synchronization plan
is still needed before remote production analysis; raw trajectories should
remain at their execution host and compact results be available on violet.

Other jobs: current 2 WT GPU replicas are ~16 ns at 270–275 ns/day, ~7.5 h left
each. Both remain on violet. Their cached parent queues still own all pending
replicas; do not launch duplicates on azalea. Current PA6 Gaussian water
probes have roughly 3 h left; the L4 intermediate geometry has only 5 optimization
steps and no credible convergence ETA, with further workflow stages queued.
No Gaussian or MD process was signaled during this capacity audit.

`mechanism-agent` now explicitly models enzyme profiles and Nyl12/NylC/NylB
comparative reasoning, with a five-case read-only qualification pilot. The PI
reports prior dissatisfaction with Gemma 4 31B and possibly gpt-oss-120b; neither
is qualified by assumption. Azalea's default Ollama service is not running;
no local model/runtime was installed. See the separate repo's updated design.

## CURRENT — direct NEB umbrella sampling launched, 2026-11-19

The user authorized US from the most favorable NEB path: the historical
force-converged regular band (cycle 477, 35.80895 kcal/mol; no verified TS).
Canonical directory is its existing `umbrella_sampling/` child, beneath
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/`.
Read that child's `README.md`, `prepared.json`, `workers/*/status.json` and
`native_retry02_startup_check.json` rather than interpreting submission as MD
completion. Eight single-thread CPU workers have initialized native Sander
with finite initial energies, the temperature ramp and both bias terms.
All eight workers subsequently reached positive timesteps; see
`launch_acceptance.json` for the immutable verified launch/progress receipt.
No production block was complete at this checkpoint.

There are 63 windows (all 24 original image centers retained, CV-only inserted
centers, nearest actual image restart seeds), two independent chemical
biases and at least 30 ps production per window after 1 ps heating/5 ps
equilibration at 300 K. Exact source Hamiltonian, 98-real-QM/charge −2,
neutral Thr N, ASP272/ASP274, and 4,165-atom movable boundary are retained.
PMF interpretation is conditional on the frozen exterior; no NEB springs
or Cartesian image restraints are carried into US. Native RK=100 for each
CV uses E=kΔq². Competing proton branches are retained and flagged, with
permutation-invariant proton observables alongside HG1-labeled distances.

Initial workers failed before any MD because native Sander truncates long
filenames. Their logs are preserved in preparation history. Retry02 uses
short local filenames/symlinks and local DISANG without changing physics or
seeds. Prepared SHA is
`3188c7f0f807e19530476ed8a9f6f97a934cdac1f6ed52a40784e1d898b87c9d`.
The GROUP-card boundary/thermostat probe and first bad card-order probe are
also preserved. An explicit heating-only temperature gate differs from the
300 K equilibration/production gate; neither establishes sampling convergence.

CPU-only analysis watcher `nyl12_direct_neb_us_analysis_20261120` publishes
`analysis/latest/` after completed blocks; disconnected overlap suppresses
global surfaces. Laptop pull now includes current diagnostics and compact
window endpoint views, without full raw trajectories. Endpoint views are
not full MD movies. Existing MD, five NEB engines and Gaussian jobs continue.
Focused new-suite checks: 74 tests passed across existing test/plotting
environments; Ruff, shell syntax and whitespace checks passed. Core new
modules measured 90–99% coverage; CLI parser branches remain less covered.
No new dependencies installed.

The requested local-agent proposal is updated in
`docs/local_agent_mechanism_project_20261024.md`: bounded tool-driven planner,
deterministic executor and independent reviewer, official non-Chinese-origin
model candidates, two reported 96 GB GPUs on azalea, and qualification on
archived cases before autonomous submission. The redesigner-agent README was
read for lessons; its implementation was not adopted. No LLM installation,
service, model download or new autonomous framework has been launched.

Separate human-in-the-loop repository created at the user's request:
`/data/jpc/mechanism-agent`, PRIVATE
https://github.com/parksjm/mechanism-agent (initial commit `fca54b8`,
branch `codex/initial-scaffold`). Read `docs/HUMAN_WORKFLOW.md` there for
PI decisions and `docs/DESIGN.md` for six specialist roles plus coordinator.
This is a documentation-only scaffold, not a deployed agent service.
Scientific jobs/data remain here; no model/dependency installation occurred.

## CURRENT — three-condition WT MD launched, 2026-11-19

The authorized matched campaign is
`qmmm/pa6_flipped_model30/10_wt_protonation_md/README.md`. Three accepted
common-box preparations cover neutral Thr-N/ASP272/ASP274, neutral Thr-N/
ASP272/ASH274 and cationic Thr-N/ASP272/ASP274 in all four copies. The
controller `wt_three_protonation_20261119` performs fully unrestrained GPU
MM then queues four independent 100 ns replicas per state on two Violet
GPUs. Frozen sources, explicit charge normalization and counterion changes
are recorded. No azalea or preexisting NEB/Gaussian jobs were changed.

The historical ASP272/ASP274 direct acylation NEB is confirmed converged
(cycle 477, maximum 35.80895 kcal/mol) but has no validated TS. Current
ASH274+QM-water bands differ in charge, QM selection, starting conformations
and boundaries; they do not isolate a protonation effect. The subsequently
authorized historical-band US launch is recorded in the newer checkpoint above.

## CURRENT — public repository and NEB/protonation checkpoint, 2026-11-19

The user explicitly requested public visibility; GitHub confirms
`parksjm/covalent-enzyme-design` is PUBLIC. Existing tracked history is
included. Credential-pattern screening found no high-confidence matches
(scope/limitations retained in the action journal). The parameter repository
was not changed.

The roadmap's new cationic-Thr comparison section records the completed
charged-N MD evidence and current five-run NEB status. Extended CG is at
cycle 260, maximum 43.642 kcal/mol and max-force/target 8.093; it is still
unconverged. Updated plots are in the existing NEB `plots/status_refresh_20261120_021148/`.
No new MD/QM campaign was launched; all active NEB engines remain untouched.

## CURRENT — student MD handoff, 2026-11-19

A transferable handoff is in `student_handoff/nyl12_ash274_md/START_HERE.md`.
The immediate hypothesis is a matched WT protonation-state effect, followed by
rebuilding the prior eight mutants with the same PA6 RESP/GAFF2 release.
ASH274 has 8.385% pooled Tier A occupancy versus 1.375% for ASP/ASP, but its
replicas are 16.77% and 0%; joint native-Thr gates do not show a corresponding
advantage. This is a testable hypothesis, not an established catalytic benefit.
Use independent 100 ns replicas, matched 300 K protocols and replica-level
uncertainty. Preserve the historical site-AE comparison and implement validated
four-site tracking before ranking LigandMPNN designs. The old mutant panel
uses legacy AM1-BCC/GAFF substrate parameters and cannot isolate protonation.
No new MD was launched while preparing the handoff. Three WT microstates are
portable prepared starts; mutant coordinate references still require rebuilding.
This MD subproject does not replace the current mechanism/TS priorities.


## CURRENT — long CG same-start replay launched, 2026-11-19

User requested a supported restart with a considerably higher budget and
explicitly allowed starting over. No usable DL-FIND checkpoint exists: the
installed Sander adapter disables dump/restart. A supported identical-start
replay is now running in the existing CG trial's `long_run/` child:
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/04_neb/spacing_tests/04_uniform_24_cg/long_run/`.
Native session: `nyl12_neb_cg_long_20261119`. Only maxcycle 250→32,000 and
Amber maxcyc 100,000→300,000 changed, yielding 224,000 evaluations (128×pilot).
All controls, chemistry, force tolerances, CG settings and 12 Å active region
are unchanged; inverse substitution verifies the rest of mdin byte-for-byte.
The completed pilot and other running NEBs are preserved.

The first 24-image profile and Max/RMS gradients match the pilot exactly.
Replay checking is scripted in `nyl12_cg_replay_check_20261119`; it compares
cycles 1–79 and then records advancement, without claiming science convergence.
Allow roughly 65 minutes to replay the old pilot. VMD and numeric comparison
sessions `nyl12_cg_long_view_20261119` and `nyl12_cg_long_comparison_20261119`
publish child `progress/latest/vmd/view_path.vmd` and `plots/comparison.png`.
The updated pull supports nested paths; recopy its laptop script. New entry
points are `scripts/prepare_long_cg_replay.py` and
`scripts/verify_neb_cg_replay.py`; run records/source hashes are retained.

## CURRENT — parallel NEB spacing tests running on violet, 2026-11-19

Four native paths now run concurrently in the existing
`09_water_relay_acylation/04_neb/spacing_tests/` directory:
01 uniform24 L-BFGS, 02 uniform48 L-BFGS, 03 adaptive24 L-BFGS and
04 uniform24 CG. All share the retained cycle268 source band, the same
DFTB3/3ob charge−1 QM102/H54/eight-link setup and 4,662-atom 12 Å movable
region. Endpoints/exterior and tolerances are unchanged. The three L-BFGS
trials are primary; CG has only a 1,750-evaluation diagnostic budget.
See the comparison README for exact numerical changes and explicit shared-H
review of one interpolated48-image seed. No chemical restraints were added.

Each native engine and VMD exporter runs in its own recorded tmux session.
Comparison monitor `nyl12_neb_spacing_comparison_20261119` writes
`spacing_tests/plots/comparison.png` and current/history JSON/CSV. Portable
movies publish beneath each trial's `progress/latest/vmd/` after a
coherent native band is captured. Laptop pull filters now include current
trial artifacts and exclude accumulated raw/history data; recopy the script.
The original NEB remains running. New profiles are very early and can be
high-energy; improved convergence and a reliable barrier are not established.

Root verification: 50 focused tests passed, Ruff/shell/whitespace checks pass.
Measured new-code coverage is85–92%; preparation recipe90%. No dependencies
installed. Initial combineduv collection lacked matplotlib; existing analysis
Python plus existinguv test tools resolved that environment mismatch.
Next compare forces/reset frequency/chemical gaps and work at50/100/250cycles,
not maxima alone. Actual native startup mapping agrees within 2.29e-7 Å for all four trials.
Separate `startup_precision_review.json` records accept mapping at 1e-6 Å;
original stricter 1e-7 diagnostic failures remain preserved. All four portable
movies have published. This is separate from force/chemical acceptance. No azalea work is authorized
while its outage is handled by the sysadmin.

## CURRENT — NEB energy-profile evolution inspected, 2026-11-19

The native ASP272/ASH274 relay-acyl NEB remains running on violet. A frozen
151-cycle snapshot and reproducible plotting script are in
`09_water_relay_acylation/04_neb/plots/` (under the existing ASH274 campaign).
`neb_evolution.png` overlays selected profiles and plots every completed
iteration's maximum and native gradient/target ratios. The maximum changes
from 52.05 to 40.97 kcal/mol; the last 20 cycles decrease by 6.31 kcal/mol.
A large intervening excursion approaches 78 kcal/mol. At cycle 151, maximum
gradient is 172.73 times target and RMS gradient is 1.181 times target.
The band is still evolving; convergence and a TS barrier are not established.

Subsequent diagnosis found repeated automatic DL-FIND kink repairs, not SCC
failure. In the installed Cartesian NEB implementation, bends above 90 degrees
trigger interpolation between neighboring images. After resets, maxima jump
35.08→82.74 kcal/mol (244→245) and 32.83→75.15 (269→270). Saved geometry also
shows a stretched transferring proton at the high-energy repaired image.
At cycle 295 the maximum is 54.67 kcal/mol, with maximum gradient 238.53 times
target. Local relaxation is repeatedly interrupted; sustained convergence is
not established. Evidence and installed-source hashes are saved in
`04_neb/plots/convergence_diagnosis.json` and the recent-band geometry diagnosis.
The run remains unchanged. Next priority is a bounded test of improved band
handling near images 9–10, preserving endpoints, Hamiltonian and the 12 Å shell.
The plot job completed in tmux; engine and monitoring settings are unchanged.
Existing laptop pull filters include the compact plots and scripts.

Azalea is back online; its outage investigation is assigned to the sysadmin.
No further azalea diagnostics or system changes are requested. Laptop pull now
uses violet's existing mirrored files without triggering an azalea refresh.

## CURRENT — exact legacy07 VMD viewer repaired, 2026-11-19

The user was opening `07_ash274_comparison/07_water_to_asp272_ti1_like/`
`vmd_progress/view_path.vmd`, not the newer09viewer. That exactfile and its
exporter are now corrected. Rootcauses: PRO1312in MMcarbonselection, partitioned
representations dropping C–H/C–N connections, oldfixedcolors, invalid shared
labelargument/return assumptions, and improperly nested camera matrix values.
Earlier headless tests were too permissive; new checks follow actual VMD API
contracts. Full visible Element Licorice plus MMcarbonorangeoverlay preserves
all146retained displayed chemicalbonds. PRO1312rendered0;54Hdynamicbonds use
resolution12. Seven07-specific distances plus BDangle have individual offsets.
Current09viewer also received a view-only refresh; no historical coordinates
were rebuilt or changed. Movie/topology hashes remain unchanged.

Exact07flatpull verified all10checksums; current09portablepull verified25root+
10VMDchecksums. Rootfocusedregression17tests,Ruff/shell/whitespace passed.
Originalfiles/sources are retained in07/history/viewer_*; backups excludedfrom
ordinarypulls. Existing07pull block delivers the correctedfile at the same
laptop path. Source in a fresh VMDsession to clear old duplicate views.
ActualGUIrendering remains unavailable here; tests check realcommandcontracts,
selections, bonded-endpoint coverage, camera math and delivery, not pixels.

Only the owned NEBinspectionmonitor was paused during sharedcode changes and
restarted at its canonicalprogress_job. ScientificSander was never signaled.
The ASP272/ASH274NEB still uses the protocol below. OfficialAPIevidence is
recorded in viewer manifests and the actionjournal (UIUCVMDmanual/source).

## CURRENT — VMD display corrected; ASP272/ASH274 NEB launched, 2026-11-19

User requested repair of the VMD script and NEB refinement of the latest path.
The stray fragment was protein PRO1312 incorrectly included by an off-by-one
substrate display selection. The renderer now derives the bonded P6 substrate
component (actual residues1313–1316), uses explicit stripped-atom mapping,
centers only the visible active-site atoms, and defaults to six distances plus
the OG1–C–O Bürgi–Dunitz angle. Dynamic bonds retain the agreed thresholds.
The transferring H is20318, not originalHG1 (currently on QM water). Re-sourcing
safely replaces only the same-bundle molecule/labels; unrelated objects remain.

Corrected coordinate-driving inspection is `progress/latest/` → update000089
in the same `09_water_relay_acylation` campaign. Movie/topology coordinates were
reused unchanged; all25root+10VMD portable manifest entries verified in a real
local rsync copy. Headless Tcl checked reload/camera/selection behavior; VMD GUI
is unavailable here, so no pixel-render claim is made.

**Native NEB launched:** tmux `nyl12_ash274_acyl_neb_20261119`, calculation
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/04_neb/calculation/`.
Runnable launcher `04_neb/launch.sh`; input `calculation/minimize.mdin`; job
metadata/logs `04_neb/job/`. Same periodic0K DFTB3/3ob,102realQM/54realH,
eightlinks,charge−1,ASP272/ASH274 and neutral betaN termini; all4,662atoms
in the original complete-residue12Å region remain movable for interior images.
QM-weighted image distances do not freeze MM forces. Endpoints are frozen
standard NEB endpoints; the original exterior remains frozen. No extra reaction
or positional restraints are included.

The band has24regularimages plus one disabled climbing slot (`neb_nimage=25`),
13trace-derived controls, `maxstep=0.05`bohr, `maxcycle=16000` and the installed
adapter's112,000energy-evaluation budget (7×maxcycle). Source50frames are38
force-converged endpoint structures plus12accepted checkpoints bridging
window025. Exact unbiased E–S and acyl endpoint ancestry, all H identities,
native slot placement, predicted24images and23midpoints passed preflight.
No control was dropped. Native alignment changes the QM atoms by ≤0.00390Å.
An earlier draft trace-order error was caught before any engine launch; that
rejected preparation is under `preparation_history/neb_seed_order_rejected_20261119/`.

**Live monitor launched:** `nyl12_ash274_neb_progress_20261119`, records in
`04_neb/progress_job/`. It pairs complete native energy reports with verified
`nebpath.xyz` bands and publishes `04_neb/progress/latest/` movie/energy/distance/
coordination artifacts. Numbered native NetCDF restarts remain unusable and
are not used. Native Sander input was accepted and energy evaluation started;
launch/startup is not NEB convergence or a verified transition state.

The laptop pull script now includes this latest NEB inspection and current
mdin while omitting raw native/history files. Refresh the laptop script copy.
Focused root verification:73tests pass; Ruff, shellsyntax and whitespace pass.
Separate preparation tests32pass, newhelpers92%coverage; VMD tests14pass,
97%coverage; monitor/helper coverage93–96%. No new dependencies installed.

## Previous checkpoint — acyl-enzyme retained; reverse scan complete, 2026-11-19

Both overnight jobs exited 0. The independent acyl-enzyme survived zero-bias
DFTB3/MM optimization (106 steps, max active gradient component 0.11697
kcal/mol/Å; gate 0.12). OG1–C=1.34456 Å, C–Nleave=3.28712 Å and C–O=1.24229 Å.
Thr-NH2, neutral departing amine, intact QM water and ASH274 proton are retained;
no unassigned/shared QM H or changed spectator-H counts were detected.
This is a force-converged endpoint candidate; curvature is not verified.

All 36 reverse windows individually force-converged, followed by a zero-bias
E–S release (38 steps; OG1–C=3.17092 Å, C–N=1.35691 Å). The acyl endpoint is
+11.06 kcal/mol relative to this returned E–S. Maximum sampled physical energy
is window_020, +50.89 kcal/mol relative to returned E–S; its bias is only
0.000108 kcal/mol. This is NOT a verified TS barrier, MEP or PMF. Window24→25
has an18.59kcal/mol energy drop with a proton-occupancy change and needs local
path refinement. Individually converged windows do not prove path continuity.

The tracked reacting H is original water H20318; original HG1 remains on water
at both endpoints and the peak. Thr Nalpha remains near NH2. Thus this scan
has not demonstrated the proposed Nalpha–water–leaving-N relay. The high
sampled maximum remains a concern for this branch.

No campaign job remains running. Next: inspect/refine window24→25 and assess
window020 as a saddle candidate with curvature/downhill tests. Preserve the
new product endpoint; a stable TI1 is not required. All evidence and the final
VMD movie/plots are in the same `09_water_relay_acylation/` campaign described
below; compact results are `overnight_summary_20261119.json`. The laptop pull
includes `progress/latest/vmd/view_path.vmd` and the final accepted-state plots.

## Previous checkpoint — acyl-enzyme endpoint workflow launched, 2026-11-18

User authorized the revised product-connected mechanism plan. A stable TI1
is optional; success requires credible reactant/product connectivity and
verified transition state(s). See the top of `docs/NYL12_MECHANISM_ROADMAP.md`.

The existing-MD relay screen is complete: 60,000 frames, 18,958 Tier A/B with
positive signed phi, 51 also passing the unchanged native-Thr gates, and zero
preorganized one-water Nalpha–water–leaving-amide-N bridges. No candidate or
protonation state was silently promoted. Full counts, controls and limitations
are under `09_water_relay_acylation/00_screen/` in the ASH274 campaign and in
manuscript Methods/Results. This is not proof that a reorganized/multiwater
relay cannot form.

**Actual launched work:** tmux `nyl12_acyl_endpoint_20261118`, launcher
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/02_acyl/job/`.
Entry point `scripts/run_relay_acyl_campaign.py`; explicit inputs in the
campaign `case.json` and `protocol.json`. It starts an independent balanced
acyl-enzyme guess from the already accepted water-QM E–S source. QM102,
54 real QM H, eight link H, charge −1, ASP272/ASH274 and neutral beta-N termini,
periodic DFTB3/3ob, electronic temperature 0 K and the original complete-residue
12 Å/4,662-atom movable selection are unchanged. No new classical acyl/TI1
parameters or E–S MM minimization were used. No positional restraints are added;
the original frozen exterior is explicit.

The initial ideal planar-ester guess clashed. The retained bounded 285-position
leaving-fragment search plus 24-value acyl torsion sweep yielded a guess with
**zero severe contacts**, maximum heavy/H displacements 2.172/2.366 Å, and
unchanged protein/water-oxygen coordinates. The selected acyl torsion is +15°
about actual bonded atoms C19920–C19919. Its acute source-like ester angle
83.97° must relax; the guess is not an optimized endpoint. All failed guesses,
sources, scores and the final construction audit are retained in
`01_endpoint_inputs/`. TIP3P's H–H constraint is excluded only from the temporary
chemical graph; the topology itself is unchanged.

The persistent script first allows 1,500 steps of product-seed preoptimization
with three explicit harmonic terms (OG1–C, C–Nleave, leaving-N H coordination).
It then removes all terms for up to 2,500 steps. Only a force-converged, fully
audited acyl-enzyme/amine endpoint triggers 36 reverse asynchronous windows
(1,500 steps each), then a zero-bias E–S release. The single scalar couples
attack, scission and label-invariant OH/Nleave-H occupancy; Thr Nalpha is free.
There is no fixed TI1 proton inventory. An unintended endpoint or force/process
failure stops the queue explicitly in `campaign_status.json`; no TS or barrier
is certified automatically. Curvature/downhill verification remains a next gate.

Live inspection session: `nyl12_acyl_progress_20261118`, launcher `progress_job/`.
It exports accepted-only plots and the water-retaining VMD movie to
`progress/latest/` in the same campaign. `scripts/pull_nyl12_inspection.sh`
includes compact reports and latest inspection files, excluding raw native
restarts/evaluations and full MD. Check both launcher and scientific records:
launch alone is not convergence or chemical acceptance.

## CURRENT — D returned to E–S; mechanistic objective revised, 2026-11-18

Read the controlling revision at the top of `docs/NYL12_MECHANISM_ROADMAP.md`.
The user asked why/whether TI1 is needed and requested a literature review.
The revised success criterion is a credible product-connected hydrolysis path
and verified transition state(s), **not mandatory identification of a stable
tetrahedral minimum**. No new calculations were launched during this review.

Both D jobs completed normally. The zero-bias release took 249 accepted steps
and returned to E–S (OG1–C 3.211 Å, C–O 1.264 Å, C–N 1.356 Å). Thr N is NH2;
water H1 binds OG1 and original HG1 binds water. ASH274 keeps its proton.
`pilot/point_D/unbiased_release/release_assessment.json` and `process.json`
record the force/chemical outcome. The finished DFT/xTB reference is under
`d_reference/{results.json,gradient_projections.json}`: both electronic methods
locally oppose further attack at fixed occupancy in this geometry. This
neither proves absence of TI1 nor validates total periodic DFTB3/MM globally.

New primary references R44–R48 are saved in manuscript references/BibTeX.
FAAH calculations explicitly released an early proton-transfer structure back
to reactants yet found complete acylation without metastable intermediates;
hASNaseIII instead found a stepwise path. Nitrogen inversion and proton access
to the leaving amide N must be included. NylC maturation and mature hydrolysis
remain separate. The current fixed high-q_occ attack pilot did not test a full
cycle in which Thr-NH3+ donates a proton onward and returns toward NH2.

Coordinate audit (`analysis/nyl12_mechanism_evidence_20261118/acid_reach.json`)
finds Nalpha–leaving-N 5.80 Å in E–S and 5.23 Å at D, incompatible with a
ready direct transfer. D's QM water is closer to the leaving N (O–N 3.35 Å,
nearest H–N 2.60 Å), motivating evaluation of a reorganized relay. A different
water role does not reinstate the earlier fixed-pocket pseudo-ring assumption.
Nyl12 F113 lacks NylC Y146's side-chain OH while retaining the backbone NH;
do not assume the same proton-delivery mechanism or an automatic F113Y rescue.

Next work is **proposed, not launched**: screen existing equilibrated models
for a leaving-group proton network, prepare one balanced acyl-enzyme endpoint,
and test an asynchronously coupled attack/cleavage/proton-transfer path without
a required TI1 minimum. Check deacylation/TS connectivity before design/PMFs.

## Previous checkpoint — pilot complete; D release and sparse reference launched

User authorized the next steps. All four occupancy-pilot points reached the
unchanged maximum biased-gradient criterion ≤0.12 kcal/mol/Å. Actual OG1–C
distances A–D are 3.1403/2.6623/2.2266/1.8537 Å; physical energies relative
to original E–S are +17.052/+20.633/+28.416/+42.763 kcal/mol. All retain
N-terminal H1/H2/HG1, no OG1 proton, intact WAT1576 and the ASH274 proton.
Point D is pyramidalized (carbonyl out-of-original-neighbor-plane height
0.278 Å), but fails the original TI1 gates: OG1–C >1.75 Å and C–O
1.292 Å <1.30 Å. It is an attack-like endpoint, not an identified TI1.
Phe113 N–substrate O worsens from 3.161 to 4.110 Å; Asn186 ND2–O remains
short (2.837 Å). The physical profile is uphill; it is not a barrier or PMF.
Reproducible audited evidence: `occupancy_callback_pilot/endpoint_audit/`.

One unbiased release is launched from D in tmux `occupancy_D_unbiased_release`.
Entry point: `occupancy_callback_pilot/release_endpoint.py`; output:
`pilot/point_D/unbiased_release/`; launcher: `job_D_release/`. Both harmonic
springs are exactly zero, with total=physical energy checked at every saved
state. Same QM102/−1, 54 real QM H, 4,662 movable atoms, periodic DFTB3/MM,
0 K electronic temperature and force criterion. Original D is retained.
`release_assessment.json` will classify the final geometry and proton inventory;
Hessian verification remains separate. Its provisional geometric classifier
is not a replacement for the original stricter TI1 gates quoted above.
No follow-up path or sampling launches automatically.

Live release export is launched in tmux `occupancy_D_release_progress`
(`job_D_release_progress/`), using `monitor_release.py --watch`. It writes
the same accepted-state conventions under
`pilot/point_D/unbiased_release/progress/`: energy/coordinate plot and
`vmd/view_path.vmd`. All 54 QM H and WAT1576 are retained, with dynamic
bonds, visible labels and standard substrate colors. The pull script includes
these small products and endpoint/reference reports, excluding native call
restarts, checkpoints, scratch and point-charge arrays.

One parallel sparse method check is launched in tmux
`occupancy_D_method_reference`, launcher `job_D_reference/`, results
`occupancy_callback_pilot/d_reference/`. It uses this exact N3/OG0 D geometry,
110 QM/link atoms and 104,824 adjusted MM charges, for
ωB97X-D/def2-SVP Force (8 CPU/32 GB), then GFN2-xTB gradient (1 thread).
It is a finite nonperiodic embedded-QM diagnostic, **not full periodic DFT/MM**.
Exact atom/charge/coordinate checks precede the electronic calculations.
Saved-gradient projections include attack/proton occupancy and restoration
of Phe113/Asn186 oxyanion contacts at fixed occupancy/attack. No new DFT
optimization or broad campaign was started. Await the release chemistry and
these matched local derivatives before prescribing another reaction path.

## Previous checkpoint — cutoff discontinuity diagnosed; pilot resumed in place

The user authorized diagnosis and repair. Retained final trials establish an
ordinary **MM 10 Å cutoff** discontinuity, rather than an SCC failure or an
occupancy-derivative defect. Pair 11272 OE2(GLU754)–11330 O(GLN757) crosses
10 Å between native calls 283/284. Observed ΔE is −0.0010040 kcal/mol for a
4.77×10⁻⁹ Å maximum move; analytic gradient work is +3.10×10⁻⁸ kcal/mol.
The screened-Coulomb plus LJ pair-removal prediction (−0.0009792 kcal/mol)
and equal/opposite force jumps agree within printed energy/force precision.
QM cutoff membership is unchanged. Evidence/scripts:
`occupancy_callback_pilot/line_search_diagnosis/`.

The repair replaces SciPy's smooth-energy Wolfe line search with the already
installed DL-FIND gradient-based L-BFGS step routine through a guarded ctypes
binding (`src/nyl12_covalent/dlfind_lbfgs.py`). Native units are bohr and
hartree/bohr; per-atom displacement is capped at 0.05 Å. Nonpositive secant
curvature triggers an explicitly logged history reset. Chemistry, QM102/−1,
12 Å/4,662-atom movable region, electronic settings, both springs, cutoff,
Hamiltonian and 0.12 kcal/mol/Å force criterion are unchanged. Energy remains
monitored; this change does not smooth or silently alter the physical model.
The installed-library quadratic test passed; 43 focused tests passed overall,
with 95% wrapper line coverage (continuation runtime coverage unmeasured).

Continuation launched in tmux `occupancy_attack_pilot_continuation`; entry point
`occupancy_callback_pilot/resume_pilot.py`, launcher `job_pilot_continuation/`.
It starts `point_A/accepted_00158_final.ncrst`, appends evaluation numbering
from 285 and accepted iteration numbering from 158 in the **same point_A**,
then gates the existing B–D sequence on force convergence. Original terminal
records move to `point_A/history/line_search_stop/` and
`pilot/history/line_search_stop/`; immutable continuation source/protocol
records are in `pilot/continuation_001/`. Check actual launcher/status and
scientific records before reporting convergence; launch alone is not proof.

**Verified repair outcome:** A passed the unchanged force gate after three
additional steps (accepted iteration 161; maximum 0.109355 kcal/mol/Å).
Its N retains H1/H2/HG1, OG1 is deprotonated, WAT1576 retains both H, and
ASH274 retains HD2. These assignments describe the **biased** endpoint, not an
unbiased minimum. B has started and shortened OG1–C toward its 2.60 Å target.
Monitor `occupancy_pilot_progress_continuation` now refreshes the same
`pilot/progress/` plots/movie; records are in `job_progress_continuation/`.
The pull script includes diagnosis, continuation protocol and archived stop
records while excluding full native-evaluation/snapshot binaries.

## Previous checkpoint — pilot stopped at point A; diagnosis required

ETA inspection found the pilot had stopped after 158 accepted iterations and
284 native evaluations (about 20 min). L-BFGS-B reported status 2, `ABNORMAL`;
the maximum active biased-gradient component was 0.277869 kcal/mol/Å versus
the 0.12 gate. Native process/SCC checks passed. The outer launcher exited 0
because the runtime records scientifically incomplete termination explicitly:
`occupancy_callback_pilot/pilot/completion.json` says `incomplete`, stopped A.
Points B–D did not run. Final accepted coordinates and all native trial calls
are retained. Do not mistake launcher exit 0 for convergence or restart the
whole pilot unchanged. The next task is a bounded line-search/energy-gradient
consistency diagnosis using retained evaluations, then an explicitly recorded
continuation from the final accepted A snapshot. No force tolerance, chemistry
or Hamiltonian has been changed. A completion ETA is not yet defensible.

## Previous checkpoint — analytic callback verified; diagnostic pilot launched

User authorized steps 1–2 immediately and step 3 after readiness. Steps 1–2
are complete: the real optimizer receives the occupancy bias energy and
analytic force, and saved three-state reference gradients have been analyzed
without new QM calculations. Read the current decision at the top of
`docs/NYL12_MECHANISM_ROADMAP.md` before continuing.

Callback proof: `qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid/occupancy_callback_pilot/callback_proof/scientific_acceptance.json`.
All seven gates passed; expected/applied initial bias 0.20 kcal/mol, four native
evaluations, two intentionally bounded optimizer iterations. This verifies
coupling, not convergence of a molecular stationary point. Native PMD inputs
remain unsupported in this DL-FIND minimization mode; no Amber installation
was changed. The adapter and coordination tests passed (31 tests; library
coverage reported 100%; runtime/proof coverage unmeasured).

Step 3 is launched in tmux `occupancy_attack_pilot`, with launcher records
`occupancy_callback_pilot/job_pilot/`. Runtime entry point:
`occupancy_callback_pilot/run_pilot.py`; serial output:
`occupancy_callback_pilot/pilot/point_A` through `point_D`.
Targets: q_occ=2.10; independent OG1–C distances 3.1251784/2.60/2.10/1.60 Å.
Same QM102/−1, 54 real QM H, 4,662 movable atoms, periodic Hamiltonian and
0 K electronic temperature. Each point stops at force/SCC failure or budget;
force acceptance is max active biased-gradient component ≤0.12 kcal/mol/Å.
All source copies, hashes, accepted iterations and raw native calls are retained.
Inspect actual `job_pilot/status.json`, stage acceptance and completion records;
launch success alone does not establish chemistry or convergence.

Accepted-step inspection products are in `occupancy_callback_pilot/pilot/progress/`:
`accepted_progress.png` (physical/total energies relative to the first pilot
physical energy; absolute bias also shown), `accepted_iterations.csv`, and
`vmd/view_path.vmd` with `path.parm7`/`path.dcd`. Frames retain WAT1576 and
come only from accepted initial/every-25/final snapshots, not line-search
trials. The pull script includes these named small products and gradient
analysis outputs, preserving workstation/laptop layout and excluding native
evaluation restarts/scratch. Exporter entry point: `monitor_progress.py --watch`.
It is launched in tmux `occupancy_pilot_progress`, with `job_progress/` records,
and refreshes every 120 s until the pilot terminates (48 h maximum). The real
two-frame export smoke check and all ten bundle checksums passed. Duplicate
label additions were removed; bonds/angles are explicitly visible. Initial
exporter sources were retained before that presentation-only correction.

Saved-gradient evidence: `reference_checks/gradient_analysis/acceptance.json`,
`results.json`, `projections.csv`, `directional_derivatives.png`. DFT/xTB
agree that attack is locally uphill along the tested directions, but disagree
on some proton-transfer slopes. No fast Hamiltonian is newly validated; native
total QM/MM slopes and embedded-QM-only slopes are kept separate.

## Previous checkpoint — occupancy tasks complete; native CV integration issue

All 25 grid-point proton inventories are audited (54 real QM H each), in
`08_activation_attack_grid/proton_occupancy_audit/`. The single authorized
unbiased release finished with force/SCC convergence and returned to E–S by
Thr/water proton exchange; no activated minimum was retained. Its native
energy is 4.7757 kcal/mol below the original E–S reference and geometry is
Tier B. Source, optimized coordinates, assessment and 13-frame VMD movie are
under `08_activation_attack_grid/diagnostic_activation_release/` (`vmd/view_optimization.vmd`).

The label-invariant coordinate and analytic derivative are prepared in
`src/nyl12_covalent/proton_coordination.py` (12 tests, 100% line coverage;
agrees with all 25 independent audit values). **Do not use native NFE/PMD
inputs with this Sander/DL-FIND minimization:** the bounded preflight found
zero applied bias instead of 0.20 kcal/mol and no PMD initialization. The
scientifically rejected integration test and corrected half-k assessment are
under `occupancy_coordinate_preflight/attempt_2/`; the coordinate definition
and analytic implementation metadata are in its parent. The installation was
not changed, and no new reaction path was launched. Next integrate the analytic
coordinate into a supported energy/gradient callback and check the saved DFT
reference gradients, following the current roadmap.

The two reference continuations have now finished normally (three total
geometries including E–S). Raw and relative energies are in
`08_activation_attack_grid/reference_checks/results.json`; gradient analysis
and scientific method acceptance remain pending. The earlier gated release
queue was not promoted; the new diagnostic was separately authorized.

## Previous checkpoint — grid complete; references resumed

The five rows and 25 reported grid points finished; no accepted TI1 was found.
Rows 1–4 ended at OG1–C 1.887/1.886/1.860/1.783 Å and physical relative
energies +28.34/+32.71/+36.62/+54.58 kcal/mol. Row 5 instead exchanges Thr
HG1 with water H1, leaving Thr–OH intact; its lower-energy points are
rearranged reactant states. Do not interpret the nominal high-q coordinate
as complete Thr-N protonation. No unrestrained release ran.

The E–S benchmarks completed, but selected-state reference exports failed on
Amber's 20-character executable field. The short-wrapper/template repair
passed an actual Sander export preflight. `ash274_reference_continuation`
runs only the activated and attack DFT/xTB comparisons in parallel, retaining
all completed E–S and DFTB3 data. Entry point:
`08_activation_attack_grid/reference_checks/continue_references.py`;
launcher: `reference_checks/job_export_repair/`. The final 25-point plot/movie
refresh was also submitted (`job_final_export/`). Method interpretation and
candidate releases remain pending. Read the dated grid result at the top of
the mechanism roadmap before further calculations.

## Previous checkpoint — activation grid and reference launches

No calculation is active in the completed 25-window water-relay path. It
transferred both intended protons but did not meet the tetrahedral endpoint
gate: OG1–C = 1.852 Å, relative physical energy +54.0235 kcal/mol, monotonically
increasing. The conditional unbiased release did not run. This branch leaves
both Asp residues protonated and Thr N neutral; it does not rule out a
Thr-NH3+/oxyanion intermediate or a route without a persistent TI1 minimum.
Read the controlling dated decision at the top of
`docs/NYL12_MECHANISM_ROADMAP.md`. The bounded 2D activation/attack grid and
sparse reference checks are authorized and started. Five sessions
`ash274_grid_row_1` through `ash274_grid_row_5` run independent q rows under
`qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid/`.
The protocol contains 25 reported grid points and 135 bounded preparation/
continuation stages, with QM102/charge −1, 4662 movable atoms, periodic
electrostatics and two independent restraints per stage. See each
`job_row_N/` and `row_NNN/row_status.json` for actual process state; the
immutable prepared protocol is intentionally not rewritten during execution.
`ash274_grid_progress` refreshes `progress/grid_heatmap.png` and per-row
`progress/row_NNN/view_path.vmd` every 120 s. At least one initial frame was
successfully exported; these rows are separate scans, not one reaction movie.
`ash274_grid_reference_queue` now runs the corrected bounded comparison
chain: E–S first, then at most three selected grid geometries, with up to two
eight-core DFT workers. Its finite nonperiodic full-PC embedding is distinct
from the periodic grid Hamiltonian. Source copies and hashes are under
`reference_checks/source_snapshot/`; launcher records are in
`reference_checks/job_queue/`; new results are under
`reference_checks/evaluations/`. The E–S DFTB3 reference finished and direct
Gaussian Force started. Analytic-gradient projection and scientific method
acceptance remain pending.

The earlier Sander/EXTERN E–S job returned process exit 0 despite Gaussian
failing at its field-query centre limit; its zero-atom FCHK and EXTERNESCF
are **rejected**. Do not use that job as a successful DFT reference.
The corrected chain retains the exact 110 QM/link atoms and 104,824 adjusted
MM charges (all non-QM atoms except eight link parents), runs Gaussian Force
directly without the excessive `Prop=(Field,Read)` query block, and requires
normal termination plus 110 atoms and 330 finite gradient values. Gaussian
energies are corrected by its own logged point-charge self energy. Original
setup failures remain preserved under `reference_checks/attempts/` and
`reference_checks/es_seed/nonperiodic_full_pc/`.

`ash274_grid_release_queue` waits for terminal grid rows and a root-reviewed
qualitative method-comparison pass, then permits at most two candidate
unbiased releases. An inconclusive comparison leaves releases pending;
it is not silently promoted. No new intermediate or barrier is claimed.
The complete water-relay result and provenance are recorded at
the end of this handoff and in `manuscript/RESULTS_DRAFT.md`.

## Previous checkpoint — water-in-QM E–S reoptimization, 2026-11-17

The corrected-shell ASH274/ASP272 E–S DFTB3/MM optimization converged and
retained Tier A geometry. Its mobile water 1576 forms a local
Thr233 OG1–HG1···water O–H2···Asp272 OD1 hydrogen-bond arrangement in the
endpoint, but not a bridge to Thr233 N or the substrate leaving nitrogen.
The full geometry and input hashes are in
`qmmm/pa6_flipped_model30/07_ash274_comparison/06_qmmm_optimization/replica_1_tier_A_seed/corrected_shell/water_1576_full_hbond_audit_final.json`.
This single geometry does not establish occupancy or catalytic proton transfer.

One water-in-QM reoptimization completed on violet in tmux
`ash274_water_qm_reopt`, under sibling `water_qm_reoptimization/`. It started
from the force-converged corrected-shell endpoint (SHA-256
`f09eac265cc74f4105080ac3bf0900f19f257942c9e407272b22082e37464bc8`),
with 102 QM atoms, charge −1, 378 electrons, eight unchanged link cuts and
4,662 movable atoms in a complete-residue 12 Å region. The four additional
boundary waters and zero starting-shell omissions are recorded in
`shell_audit_water_qm.json`; the executable input, script and hashes are in
`preparation.json`. Sander exited 0; DL-FIND converged at cycle 248. The
endpoint retained Tier A E–S geometry (d1 3.100 Å, d2 3.063 Å, d3 2.959 Å,
attack angle 116.44°, signed phi +98.75°). Thr OG1–HG1, both water O–H
bonds, and ASH274 OD2–HD2 remained intact. The local Thr HG1···water O and
water H2···Asp272 OD1 contacts remain oriented (1.932 Å/165.5° and
1.730 Å/172.7°, respectively), while water O···substrate leaving N is
4.107 Å. `water_1576_full_hbond_audit_final.json` and
`scientific_water_assessment.json` give hashed geometry and separate gates.
The generic `acceptance.json` contains a stale inherited initial water O···OG1
distance; use `scientific_water_assessment.json` for the actual start/end
comparison. This result supports one bounded water-assisted Thr activation
and OG1 attack test toward TI1, with explicit proton inventory and an
unrestrained endpoint release. It does not establish water occupancy,
proton transfer or a complete relay to the leaving group.

## CURRENT — bounded mechanism workstreams, 2026-11-16

PROPKA 3.5.1 was found already installed; a three-frame Tier A protein-only
screen completed in tmux without adding dependencies. See
`docs/propka_asp_pair_pilot_20261116.md` and
`qmmm/pa6_flipped_model30/09_propka_static_tier_a/screen_attempt1/` for
methods, all four beta-site predictions, input hashes and limitations.
The substrate-facing site E estimates vary with structure and do not
consistently favor Asp272-COOH. This is a hypothesis screen, not a
protonation assignment.

The matched six-replica late Tier A persistence and three-replica candidate
water-oxygen screen are in `docs/asp_pair_md_water_screen_20261116.md`.
The water scan finished in tmux; per-frame event CSVs and 10 ps identity
streaks are under `09_propka_static_tier_a/water_screen_attempt1/`.
Oxygen proximity alone did not establish a hydrogen-bond relay; the
inspected ASH274 Tier A frame has an OG1-HG1-to-water hydrogen bond but
the water hydrogens do not point toward the Thr N terminus. Next inspect
explicit hydrogen-bond geometry across candidate frames before selecting
a follow-up QM/MM E–S snapshot.
The two-leg hydrogen-bond scan subsequently finished in `hbond_job/` and
`hbond_screen_attempt1/`. Among Tier A frames, complete
OG1-HG1→water→Thr-N geometry occurred in only 1/238, 6/3886 and
32/1677 frames across the three screened replicas; the longest
same-water streak was three 10 ps samples. No Tier A frame had a
complete two-leg bridge to the substrate leaving-group N. See the
dated water-screen note for criteria and limitations.

The two ASH274/ASP272 300 K MD replicas completed 100 ns each with process
exit 0. A matched six-replica Tier A and Asp-pair analysis completed in
`qmmm/pa6_flipped_model30/08_asp_pair_md_comparison/analysis_attempt2/`;
its immutable inputs and tmux log are in `inputs/manifest.json` and
`job_retry/`. The first analysis launch failed because the project `.venv`
lacked MDAnalysis; the successful retry used the previously validated
`enzyme-substrate-modeling/.pixi/envs/analysis` interpreter without adding
a dependency. Tier A counts per 10,000 frames are ASP272/ASP274 37 and 238;
ASH272/ASP274 3886 and 1724; ASP272/ASH274 1677 and 0. Site-1 median
closest Asp-pair O–O distances are 4.89/4.90, 2.56/2.56, and 2.73/2.65 Å,
respectively. These are replica-specific associations, not pKa estimates or
proof that either Asp is protonated in the enzyme. Inspect time blocks,
substrate register, proton geometry, and water before selecting follow-up
snapshots. The Asp-pair comparison is the current E–S decision gate.

The direct-route DFTB3/MM dimer converged at cycle 99. Its block-LOBPCG
low-mode calculation also finished, but approximate eigenvalues −0.227 and
+0.033 Eh/(Å² amu) have residual norms 0.051 and 0.061, versus the 10⁻⁴
target. Thus Hessian index one is **not verified**. The near-identical prior
candidate has downhill paths toward E–S and acylenzyme; this is supportive
but does not establish the new candidate's Hessian index. The bounded
operator check below was performed before choosing another mode verifier.
The 35.81 kcal/mol NEB peak remains an approximate
PES image energy, not a verified barrier.

A bounded ten-gradient operator diagnostic then completed in tmux; see
`docs/ts_hessian_operator_diagnosis_20261116.md` and the TS folder's
`hess_op_20261116/diagnosis.json`. Repeatability was exact at the tested
step, chemical-mode step sensitivity 0.93%, and bilinear asymmetry 0.37%.
The native force operator appears numerically stable enough for this
test; full Hessian index remains open. Do not rerun the same stalled
unpreconditioned low-mode calculation unchanged.

## HISTORICAL — two bounded mechanism workstreams, 2026-11-14

**Protonation comparison.** Violet tmux `nyl12_ash274_md_20261114` is
running the four-beta-chain Asp274–COOH / Asp272–COO⁻, neutral-Thr-N
alternative at 300 K, with two independent GPU replicas targeted to 100 ns.
Run root: `qmmm/pa6_flipped_model30/07_ash274_comparison/`. The staged
Asp274 proton is on OD2 in each beta chain. Dry charge changed by exactly
+4 e from the ASP/ASP parent, the substrate parameters and neutral Thr
charges were retained, and the final ion count was audited (110 Na⁺,
106 Cl⁻). Initial MM preparation passed readability and process gates; MD
sampling and scientific acceptance are pending. The existing ASP/ASP and
ASH272/ASP274 campaigns both originated from the same static model30
protein/substrate pose and used the same 300 K MD protocol; they are
independent solvent/seed realizations, not paired solvent configurations.
After ASH274 completion, apply one signed-phi/Tier and Asp-pair/H-bond
analysis to all six 100 ns replicas. Compare persistence by replica and
time blocks, not by treating correlated frames as independent experiments.

**Direct-route TS.** Violet tmux `nyl12_direct_ts_image08_20261114`
completed a 0 K DFTB3/MM DL-FIND dimer refinement from image 8 of the
converged 24-regular-image NEB. Run root:
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/ts_image08_attempt1/`.
The valid seed restart was reconstructed from native NEB XYZ; the known
invalid native image restart writer was not used. The mode is the image
9-minus-7 tangent displaced 0.01 Å on the rotation mask. The original 12 Å
region (4,165 atoms), 98-atom QM region, charge −2, and zero electronic
smearing are retained. Peak-image local clash audit found zero severe
contacts. The scripted run allows up to 1,000 cycles and stops on first
DFTB SCC failure. It converged at cycle 99 without SCC failure, meeting the
native maximum-gradient criterion (9.2423×10⁻⁵ versus 1×10⁻⁴). Final
OG1–C = 2.0805 Å, C–leaving-N = 1.5403 Å, and HG1–leaving-N =
1.2184 Å; Thr N retains H1/H2. The new candidate differs from the earlier
downhill-tested dimer structure by only 0.0058 Å RMSD over the 4,165
movable atoms (0.0084 Å over mapped reactive atoms). This supports the
same local chemical saddle hypothesis; it does **not** verify Hessian index
one. A distinct block-LOBPCG finite-difference low-mode diagnostic is
running in tmux `nyl12_ts_block_modes_20261114`, with 0.05 Å√amu
displacements and a hard 160-HVP cap under
`ts_image08_attempt1/mode_probe_block_0p05/`. Require one reproducible
negative mode, no additional negative low modes, and downhill connectivity
before calling this a true first-order TS.

The task priority is to complete these two gates before opening another
mechanism or temperature campaign. Assay temperature remains 348.15 K;
current MD stays at 300 K as requested. Job records and
`output_agent_actions/20261021.jsonl` retain commands, hashes, and status.

## CURRENT — endpoint chemistry, movies, and path audit, 2026-11-14

The matched ASH272 MD SI panel already exists at
`qmmm/pa6_flipped_model30/06_ash272_comparison/04_md_validation_20261113/results_attempt1/si_matched_panels/ash272_si_trajectory_panel.png`
(PDF alongside it). All four converged ASH272 QM/MM optimizations now have
stripped VMD bundles in their respective `vmd_optimization/` directories;
the laptop pull script includes them. Their audited endpoints all pass the
geometric Tier A filter. Replica 1 Tier A has the most favorable *geometry*
of this set (d1 2.778 Å, d2 2.961 Å, d3 2.981 Å, alpha_BD 103.84°,
signed Herschlag phi +85.72°, CB–OG1–C 113.45°). However, **none retains
the assigned ASH272 state**: in each endpoint the initial Asp272 HD2 is
1.02–1.05 Å from Asp274 OD2 and 1.49–1.65 Å from Asp272 OD2. The
Thr233 HG1 and substrate scissile C–N bond remain intact. See
`06_qmmm_optimizations/endpoint_chemistry_geometry.json`. Treat these as
geometrically promising, Asp274-protonated *outcomes*, not accepted
ASH272-protonated E–S structures or proof of a catalytic protonation state.

The converged direct-route NEB's 24-image VMD movie, exact relative PES
profile, and reactive-distance plot are under its `vmd_final_24img/final_24img/`
bundle, with laptop pull support. Maximum image energy is 35.80895 kcal/mol
relative to image 1 (image 8, cycle 477). Identity and endpoint checks pass;
the largest adjacent changes are OG1–C 0.388 Å, C–N 0.290 Å, and HG1–N
0.515 Å. A separate full nonbonded-clash sweep is outstanding. The older
~46 kcal/mol number is a dimer candidate relative to a separately minimized
E–S structure, so the two values are not comparable barriers. Neither is a
verified activation free energy or a confirmed first-order saddle point.

The bounded dimer precision diagnosis found energy and force curvature agree
within 0.73% at ±0.05 Å along the tested mode. This argues against a gross
force-sign or derivative-precision error in that direction; the global
low-mode eigensolver still produced no converged eigenpairs after 800 Hessian
vector products. Next: refine a chemically accepted NEB maximum with a TS
optimizer and independently verify exactly one negative mode and downhill
connectivity. See `provenance/ts_verification_diagnosis_20261114.md`.

Continue MD at 300 K for now. The assay was 75 °C (348.15 K), 6 h. A future
348.15 K PMF/rate comparison would be more directly assay relevant. There is
no universal viable barrier cutoff: for kappa=1, Eyring estimates at 348.15 K
give roughly 0.11 s⁻¹ for ΔG‡=22 kcal/mol, 0.0062 s⁻¹ for 24, and
~1.1×10⁻⁶ s⁻¹ for 30. These are illustrative unimolecular rates, not
predicted assay conversion, and the NEB PES maximum cannot be substituted
for ΔG‡.

## PRIOR — completed processes requiring scientific acceptance, 2026-11-14

No tmux sessions remain at this status check. The four ASH272 DFTB3/MM
optimizations exited with process status 0 and each native mdout reports
`Optimisation converged` with maximum-gradient criterion met. Their existing
scientific-acceptance JSON files still say pending/running; **do not interpret
them as chemically accepted** until proton retention, scissile-amide bonding,
Thr233 orientation, substrate register, and final Tier geometry are measured
from each final restart and the records are updated. Inputs and results are in
`qmmm/pa6_flipped_model30/06_ash272_comparison/06_qmmm_optimizations/`.

The 25-image trace-seeded 0 K NEB also exited 0, and native DL-FIND reports
`Optimisation converged` after cycle 477. Its last complete 24-regular-image
energy report has a maximum at image 8, approximately 35.81 kcal/mol relative
to image 1. This is a **converged optimizer report, not yet an accepted MEP or
barrier**: validate image-wise proton/connectivity, spacing, clashes and
endpoint chemistry, then consider CI refinement/TS search. The new
80-vector-subspace Hessian low-mode probe exited 0 but exhausted its 800-HVP
(1,600-gradient) budget without two eigenpairs; first-order saddle character
remains unverified. Do not raise the same numerical budget again without a
different eigensolver/precision diagnosis.

Assay-temperature record: the user confirmed the experimental protocol was
**75 °C (348.15 K) for 6 h**; the later “65 K” message was not the assay
temperature. The pH 8 Tris measurement temperature remains unknown. MD was
at approximately 300 K. No simulation temperature was changed by this note.

## PRIOR — ASH272 charge-aware QM/MM optimizations launched, 2026-11-14

The active ASH272 QM region was remapped by residue position and atom name from
the prior ASP topology. Its audited inventory is 99 real atoms + 8 link H,
formal QM charge -1, 368 electrons, and singlet multiplicity. The source ASP
region is 98 real atoms + 8 links, charge -2, and also 368 electrons; the
added catalytic Asp272 HD2 is included. This corrects the earlier mix-up of
QM nuclei with electrons. The remapped QM mask is recorded in each case's
`preparation.json`.

The first two native Sander starts stopped before gradients because DL-FIND's
AmberTools26 `MASKLEN=2560` truncated a 2,913-character active atom mask. Those
logs and metadata are retained in each Tier A case's
`attempt_01_parser_mask_truncation/`. The 4,166 movable atoms were found to be
exactly 441 complete residues. A 1,438-character residue mask was checked with
ParmEd to select exactly the same 4,166 atoms; the preparation code now fails
if it encounters a partial residue or if the compressed mask changes
membership. The focused tests pass (8); Ruff check and format checks pass.

Four unrestrained native DFTB3/MM optimizations are now active, initialized
from the four accepted, completed MM endpoints. Tier A retry sessions are
`ash272_qmmm_replica_1_tier_a_retry2_20261114` and
`ash272_qmmm_replica_2_tier_a_retry2_20261114`; each mdout reports QM charge
-1, 99 QM atoms, and 4,166 DL-FIND active atoms. Tier B sessions are
`ash272_qmmm_replica_1_tier_b_20261114` and
`ash272_qmmm_replica_2_tier_b_20261114`. Their command/input hashes and process
records are alongside the run files under
`qmmm/pa6_flipped_model30/06_ash272_comparison/06_qmmm_optimizations/`;
Tier B attempts use the project launcher under each case's
`attempt_01/job/`. At this handoff all four sessions were active and scientific
acceptance was still pending. Check each process status, then assess force
convergence, retained proton/scissile-amide chemistry, and endpoint geometry
separately. No optimization outcome or chemical acceptance is claimed yet.

## CURRENT — 300 K SCC diagnostic passed cycle 60, not a TS, 2026-11-13

Tmux `model30_dimer_scc300_cycle60_20261113` completed at 15:43:21 UTC with
process exit 0 and **zero SCC-failure warnings**. This passes the original 0 K
first-failure point (cycle 50). The job is a bounded electronic-convergence
diagnostic, not a validated TS. DL-FIND says `WARNING: OPTIMISATION HAS NOT
CONVERGED`; cycle-60 maximum gradient is 0.040982 vs target 0.0001, maximum
step 0.1 vs target 0.0004 (DL-FIND internal/atomic units). Midpoint energy
changed from -718.1684607 Eh at the initial guess to -717.9685621 Eh at
cycle 60, +125.44 kcal/mol. These are *search* energies, not a reaction barrier
or PMF. Final saved geometry is product-like: Thr929 OG1–PA6 carbonyl C
1.401 Å, Thr HG1–leaving-group N 1.024 Å, C–leaving-group N 2.107 Å, and
OG1–HG1 2.131 Å (measured with the saved Amber topology/restart). A true
saddle and downhill connectivity remain unestablished. The run's
`pilot_status.json` field `scientifically_accepted=true` **only encodes process
exit 0 plus no SCC warning**, not geometric or saddle acceptance. Do not
continue the same dimer search blindly; inspect mode, forces, geometry and
search trajectory first. Full artifacts are under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/06_dimer_ts/scc300_pilot_20261113/`.

## PRIOR — dimer TS interrupted for DFTB3 SCC failures, 2026-11-13

User requested immediate interruption. The original dimer session
`model30_true_dimer_ts_20261113` was stopped by tmux Ctrl-C; the two ASH272
GPU MD replicas remain running. Its native `dimer.mdout` contained 127 SCC
nonconvergence warnings. Amber's implementation explicitly continues with
inaccurate energies and forces after these warnings, so neither the interrupted
search nor its apparent energy profile is a validated TS result.

Scripted one-gradient probes and hashed inputs are under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/06_dimer_ts/scc_diagnostic_20261113/`.
The original TS guess converges at `scfconv=1e-10`, but a saved later geometry
fails at `1e-10`, Amber's default `1e-8`, and `1e-8` with Broyden iteration
budget raised from 70 to 200. At the saved geometry, `dftb_telec=300 K`
converged at both `scfconv=1e-8` and the original `1e-10`; the latter is a
single-variable SCC change. On the original TS guess, 300 K versus 0 K changes
the printed QM energy by about 0.0000024 kcal/mol. This same-geometry check
does not establish negligible energy differences everywhere on the path.
The first diagnostic attempt omitted `AMBERHOME` and failed before evaluation;
its logs remain under `scc_diagnostic_20261113/job/` and `results/`. Corrected
attempts are retained separately.

Tmux `model30_dimer_scc300_pilot_20261113` is a bounded 20-cycle pilot from
the original TS guess, keeping the original method, SCC tolerance and dimer
settings while setting only `dftb_telec=300 K` and the pilot cycle limit.
`scripts/run_dftb_scc_pilot.py` watches native output and terminates sander at
the first SCC warning because sander itself treats that warning as nonfatal.
Its results and process status are under `06_dimer_ts/scc300_pilot_20261113/`.
Successful process completion would still not establish a saddle; assess SCC
retention, force convergence, one negative mode and downhill connectivity
separately. Do not restart from the failed original late checkpoint without a
new geometry audit.

**Pilot outcome (2026-11-13):** the 20-cycle job completed with process exit 0
and no SCC-failure marker. It did not converge as a saddle: final reported
maximum gradient 0.27710 versus target 0.00010, and the final DL-FIND output
says `WARNING: OPTIMISATION HAS NOT CONVERGED`. Retain the 300 K smearing as a
promising electronic-convergence setting, but do not treat the pilot endpoint
as a TS or launch a full search from it before a geometry/optimizer audit.
The same-geometry 0 K versus 300 K energy difference at the *initial* TS guess
was ~0.0000024 kcal/mol; this is not a path-wide validation.

Amber keyword audit: in the installed AmberTools 26 source, `itrmax` already
defaults to 1000; `dftb_maxiter` defaults to 70 and controls Broyden resets.
`scfconv` defaults to 1e-8, but using that value did not converge the saved
late dimer geometry. `dftb_disper` toggles a dispersion correction, not SCC
convergence, and `dftb_3rd_order='PA'` selects an additional parameterized
third-order option rather than a solver aid; neither should be added ad hoc.
The TS input has `imin=1,ntc=1`, so QM SHAKE and an MD timestep are irrelevant
to this minimization. A 100 K smearing test remains optional; 300 K already
passed a bounded electronic check. Source paper: Seabra et al., *J. Phys.
Chem. A* 2007, DOI 10.1021/jp070071l (reported 1e-9 Eh SCC tolerance for
alanine-dipeptide QM/MM dynamics, not a prescription for this active site).

Further DL-FIND/source audit: `delta` is the dimer separation, but the installed
Amber default is already 0.005 bohr (~0.00265 Å), exactly the proposed value;
the run used that default. Its `maxstep=0.1` bohr (~0.053 Å) is already at the
low end of the proposed range. `dimer_maxrot=20` is the installed default, and
the initial dimer mode was explicitly derived from NEB images 7→9. The first
SCC failure in the original run was at NSTEP 255, dimer cycle 50, during a
**midpoint** evaluation after the preceding rotated endpoint converged. Thus
the claim that a displaced endpoint first hit a clash is not supported by the
native log; a geometry/clash audit remains necessary. The active region already
contains 4165 atoms in a 12 Å shell, with only the exterior fixed. Reducing
it to 4–6 Å would be a different mechanical approximation, not a generic SCC
fix, and would repeat prior over-restraining concerns. The user explicitly
directed on 2026-11-13 that the 12 Å movable shell must not be restricted
further. In this build, DFTB3/3OB
already loads `3ob-3-1` and sets H–X damping exponent `zeta=4.0`; no
`dftb_hcorr` namelist keyword was found. Do not add the older `PA` third-order
option without an explicit Hamiltonian comparison. Finite-T occupations may
alter energies/forces near a narrow gap, so any 300 K search must be checked
against 0 K at candidate stationary geometries and along the local path.
The proposed `PA` file (`DFTB_3RD_ORDER_PA.DAT`) is absent from the installed
`3ob-3-1` parameter directory; adding that option now would fail at setup.

**Restart provenance correction:** `calculation/optimized.ncrst` was written
again at output lines 5068, 5305, 5582, 5899, 6772, 7090 and 7244, all
after the first SCC failure at line 4879. It is *not* the clean cycle-49
checkpoint, even though line 4872 says a checkpoint was saved then. A bounded
contact screen on the current restart found no obvious short nonbonded clash,
but does not diagnose the first failure because that file was overwritten.
The action journal retains this correction to the delegated read-only audit.
Tmux `model30_dimer_scc300_cycle60_20261113` now reruns from the original TS
guess, guarded to stop at the first SCC warning, to test the setting through
the original cycle-50 failure point. Its job and calculation directories are
under `06_dimer_ts/scc300_pilot_20261113/` as `job_cycle60/` and
`calculation_cycle60/`. Do not infer the result until native output completes.

## CURRENT — ASH272 MD queue and real saddle search launched, 2026-11-13

User explicitly authorized both launches. ASH comparison root:
`qmmm/pa6_flipped_model30/06_ash272_comparison/`.
Tmux `ash272_2rep_100ns_20261113` completes unrestrained MM preparation then
runs two independent GPU0/1 replicas using the original heating/equilibration
and 10+90 ns protocol. Native GPU minimization verified past cycle6300.
All four Asp272 residues are ASH with original OD2 protonated; Asp274 remains
ASP and all four Thr233 N-termini remain NH2. Dry charge changes by
4.0000000066 e; inherited -0.00400495 e residual is preserved. Salt recount:
106 pairs plus four neutralizing Na, total110Na/106Cl. Substrate released
charges/types/nonbonded parameters match the original topology exactly.
Inputs and acceptance: `00_system/preparation_acceptance.json`; queue/source
snapshots under `03_md/`; process record under `overnight_job/`.
The queue waits for each replica's own successful10ns before its90ns extension.
Do not describe a queued replica as actively integrating before native output.

Tmux `model30_true_dimer_ts_20261113` launched a real native DL-FIND dimer saddle
search in `qmmm/pa6_flipped_model30/05_neutral_reaction_path/06_dimer_ts/`.
Final complete NEB cycle499 maximum is image8, E=-718.1684607Eh, about24.76
kcal/mol above E-S. Terminal XYZ reconstruction was independently checked
against that native energy and matched at printed precision; results retained
under `energy_validation/`. Dimer direction from images7->9, same98QM+8links,
charge-2, same4165movable atoms and fixed exterior; no reaction/position bias.
Native dimer curvature/rotation output verified. Search convergence, one
imaginary mode, and downhill connectivity remain unestablished. This job retains
the original dianionic Asp state; it is separate from the new ASH MD comparison.

## CURRENT — neutral model30 MD inspection movies exported, 2026-11-12

The completed 100 ns neutral replicas now have compact laptop-viewable bundles:
`qmmm/pa6_flipped_model30/03_md/neutral_replica_01/vmd_stripped_100ns/`
and the corresponding `neutral_replica_02` path. Each bundle concatenates the
10 ns and 90 ns production segments in chronological order, retains 401 frames
(approximately 0.25 ns spacing plus the exact final frame), autoimages around
PA6 residues 1313-1316, and retains all 19,967 protein/substrate atoms while
removing bulk water and ions. Each DCD is 91.7 MiB. These downsampled files are
for visualization only; quantitative analyses must use the original NetCDF.
`view_production.vmd`, the matching stripped topology, source hashes, cpptraj
inputs/logs, README and checksums are retained in each bundle. The laptop pull
script includes the bundles but continues to exclude the raw 12+ GiB-per-replica
production trajectories. Attempt 1 failed before reading coordinates because
`libsander.so` was absent from its runtime library path; the preserved attempt-1
record is under `03_md/jobs/neutral_vmd_export_20261112/`. Corrected tmux attempt
2 supplied the explicit Amber library path and completed with exit status zero
under `03_md/jobs/neutral_vmd_export_20261112_attempt2/`.

## CURRENT — matched E-S to acylenzyme NEB running, 2026-11-11

User authorized NEB once endpoint ready. Final unrestrained product accepted:
SHA256 `d4ee10b57fc5844acbfa6e20af3718ba6bf0adac6a3d77b937843cee9dea8d46`;
original optimized E-S `d9479cd85293e4c1ba4a877e0a77326b29d8e814f42e0575a7a6e29b88b222e6`.
Same topology SHA256 `448560852cfe65e86e6ef094f2bfc9551029dd1b3cf9847b17fb1a89116a8f63`.
Root: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/05_es_to_acylenzyme_neb/`.
Tmux `model30_es_to_acylenzyme_neb_20261111`, Sander PID118669 at launch;
native output confirms frozen endpoints, image setting17 (16 regular band
images plus reserved CI slot), no climbing and no image freezing during
optimization. OneCPU DFTB3/3OB periodic QM/MM,98QM+8links,charge-2;
same4165 movable atoms, exterior fixed. `nebk=0.01`,maxstep0.1,maxcycle1000,
tol1e-4,tole1e-5,ntr0,nmropt0,neb_climb_test0,neb_freeze_test0.
Initial path uses original E-S, directed windows02,04,...20, then accepted
product; failed21 excluded. All12 seed images have identical box/atom order,
zero frozen-exterior mismatch, original ThrH1/H2 and substrateHN retained.
Native interpolation redistributes seeds; actual E-S energy -718.2079132Eh
and subsequent image evaluations verified. No converged barrier claimed.

Tmux `model30_es_acyl_neb_monitor_20261111` runs campaign `monitor.py`:
energy exports every 60 s, verified native XYZ reconstruction for first VMD band
then every10cycles. Entries `plots/neb_profile.svg`, `plots/neb_evolution.png`,
`vmd_progress/latest/reaction_path.vmd` appear after complete reports.
Numbered native NetCDF NEB restarts remain INVALID due installed writer defect;
only coherent native `calculation/nebpath.xyz` reconstruction is used for VMD.
Inputs, scripts, source snapshots, exact command/environment and hashes retained.
Preparation/monitor Ruff passed; runtime source/input/proton gates passed;
new artifact-script coverage unmeasured. Monitor launches no candidate releases.

## PREVIOUS — endpoint-first acylenzyme hypothesis optimizing, 2026-11-11

User authorized endpoint-first construction after the direct composite coordinate
became chemically degenerate. Pinned source: direct path window20 SHA256
`d57540b6bd5ede95c360174ea4d1bd92c80625355d3c2ccd13f321443b9674ad`.
New isolated directory: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/04_direct_acylenzyme/acylenzyme_endpoint`.
Fourteen cheap deterministic geometry candidates preserved; selected C-N2.85A,
ester1.40A, carbonyl1.23A, substrate N-HG1 1.03A, neutral ThrNH2 and original
substrate HN retained. Max displacement1.4375A across42 atoms, frozen exterior0.
Zero heavy-heavy overlaps under declared0.65sum-VdW rule; five H contacts retained.
Starting OG1-C-acyl angle136.21deg is a seed limitation, not an optimized result.

Reaction-unrestrained DFTB3/MM running in tmux
`model30_acyl_endpoint_unrestrained_attempt2_20261111`, actual Sander PID103870
verified with cycle7 energy -717.9966818Eh and RESTRAINT=0.0000. Same original
compact4165 movable shell,98QM+8links,charge-2; outer atoms frozen. No proton
holding restraints, restrained preoptimization, artificial guard, or NEB queued.
First startup exit127 (missing explicit Amber shared-library path) is preserved
in `job/`; corrected original Amber environment recorded in `job_attempt2/`.
Calculation: `unrestrained_optimization_attempt2/`. Geometry/provenance and source
snapshots retained. Initial5-frame VMD bundle in `visualization/latest` includes
seed, saved trajectory, restart and all reactive-bond labels; chemistry retained
at initial checkpoint but native force convergence not yet established.
`monitor_direct_acyl_endpoint.py` records independent full-H chemistry and native
convergence gates. Builder tests23 passed, module coverage95%, Ruff passed.
Next: inspect completed reaction-unrestrained endpoint chemistry and convergence;
do not infer a minimum or barrier from process exit alone or start NEB before review.

## CURRENT — direct substrate-N acylation path launched, 2026-11-11

The user authorized a direct four-center E–S-to-acylenzyme path. Tmux
`model30_direct_acylenzyme_20261111` and Sander energy output were verified.
Root: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/04_direct_acylenzyme/`.
The HG1 acceptor was substrate N19929, not Thr N13905. The composite
coordinate had 39 targets from −3.273755808 to +2.5 Å in 0.148044 Å steps.
The planned acyl-like candidate release removed all reaction bias. The path
used the original E–S model, 98 QM atoms plus eight links, charge −2, and
4,165 movable atoms. Original amide HN19937 and Thr H1/H2 were monitored.
VMD exports showed HG1–substrate N and dynamic C–N cleavage; the label watcher
ran in `model30_vmd_substrate_n_labels_20261111`. No other jobs were stopped.
The concurrent TI1 path was still running at verification; stepwise path 02
required a separate check. No completed-window or scientific-success claim
was made at launch.

## CURRENT — stepwise attack and concurrent PT/attack running, 2026-11-11

Both single-CPU DFTB3/MM jobs launched in tmux under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/`, with 98 QM atoms plus
eight links, charge −2, and 4,165 movable atoms. Path
`02_nucleophilic_attack` (`model30_neutral_attack_20261111`) started from the
released proton-transfer endpoint and had 16 OG1–C windows to 1.50 Å, no N–H
restraints, and a planned reaction-unrestrained release. Path
`03_concurrent_pt_attack/work` (`model30_concurrent_pt_attack_20261111`)
started from original E–S and had 29 asynchronous composite targets from
about −4.06 to 0.0 Å with 0.145 Å maximum step. Its planned early release
required a converged TI1-like candidate. Chemistry gates tracked component
distances, all relevant protons, and scissile C–N; no automatic umbrella
sampling was queued. Per-window energy/distance plots and VMD were enabled.
Sander processes were verified running at launch, with no completed-window or
scientific-outcome claim. Focused tests passed (37 stepwise, 25 concurrent);
Ruff passed; new coverage was unmeasured. Protocol:
`docs/neutral_stepwise_concurrent_paths_20261111.md`.

## CURRENT — proton-transferred endpoint retained, 2026-11-11

`05_neutral_reaction_path/01_direct_pt` completed 19 converged windows and a
reaction-unrestrained release. The endpoint had N–HG1 1.1105 Å, OG1–HG1
1.5629 Å, and OG1–C 3.0397 Å; H1/H2 remained on N. DL-FIND converged within
the 4,165-atom movable region; outer atoms remained fixed. Released total
energy was −718.186312526996 Eh versus −718.207913180869 Eh for original
E–S, a difference of +13.5546 kcal/mol. This was a proton-transferred
stationary candidate, not a Hessian-validated minimum. The monotonic driven
profile did not locate a TS. The latest VMD bundle had 20 frames (19 windows
plus release), and its plot referenced PT01 rather than original E–S. The next
step was OG1–C attack from that endpoint while tracking all N–H/OG1–H and
scissile C–N distances, followed by release to test TI1 retention.

## CURRENT — neutral PT path running; charged MM gated out, 2026-11-11

Tmux `model30_neutral_direct_pt_20261111` had verified Sander output in
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/01_direct_pt/`. The plan
was 19 direct OG1-to-Thr-N proton-transfer windows and an unrestrained release,
using DFTB3/MM, 98 real QM atoms plus eight links, charge −2, and the original
4,165-atom movable region. A long atom mask failed before optimization; an
equivalent compact residue mask fixed it. Failure evidence was retained under
`clutter/initial_mask_failure/`. Per-window VMD and plots were enabled and the
laptop pull rules were updated. Charged frame 7013 MM reached force RMS
0.0094117 but failed Tier/native geometry gates, so no charged QM job launched.
See `docs/md_model30_charged_screen_20261111.md`.

## CURRENT — neutral optimized E–S accepted; new paths authorized, 2026-11-11

All four model30 replicas completed 100 ns and were screened. Neutral replica
2 frame 27 (0.27 ns; early transient) completed MM and 12 Å DFTB3/MM
optimization. It passed native convergence and Tier A geometry: d1 3.148 Å,
signed Herschlag phi +100.035°, and Bürgi–Dunitz angle 105.491°. Charged
replica 2 frame 7013 (70.13 ns) was an isolated Tier B/native pass authorized
for MM; tmux `model30_charged_frame7013_mm_20261111` was launched. Evidence
is under `qmmm/pa6_flipped_model30/04_md_screen_20261111/`. At this
checkpoint, the direct OG1-to-Thr-N path was being prepared, not yet launched.
The baseline had 98 QM atoms plus eight links, charge −2, and the original
12 Å movable region. All N–H and OG1–H distances needed monitoring.
Assay conditions were documented as crude lysate, 5 mg polymer powder,
200 mM Tris, reported pH 8, 75 °C for 6 h; pH adjustment temperature was
unknown. MD ran at 300 K; no high-temperature MD had launched.

## CURRENT — matched MD underway; NEB ended, 2026-11-10

The then-current NEB ended at cycle 437 by reaching its energy-evaluation
limit; it had not converged, and its last profile peaked at 86.939 kcal/mol
at image 8. No agent killed it. Two neutral model30 replicas were in
unrestrained 10 ns production on violet under `03_md/neutral_replica_01` and
`_02`. Two all-four-beta-NH3+ replicas were authorized, with charged setup
on azalea underway. Neutral-replica heating (100 ps) and staged NPT (400 ps)
passed. The RMS-force 0.042 endpoint was sufficient for MD preparation; no
extra minimization was required. Job logs, not this historical status,
determine whether any run is now active.

## CURRENT — model30 MM endpoint passes geometry; force target unmet, 2026-11-10

The neutral-all-four-beta model30 build and unrestrained GPU MM run exited 0
after 10,000 cycles. RMS force was 0.042036 versus the 0.001 target, so it
was not force-converged. Endpoint `01_mm_minimization/endpoint_complex.pdb`
retained flipped phi +97.948° and Tier A/native geometry: d1/d2/d3
2.882/2.978/2.871 Å, Bürgi–Dunitz angle 93.082°, and N-terminal-to-leaving-N
distance 3.940 Å. OG1/N shifted 0.263/0.698 Å; psi/chi1 changes were
−21.662/+5.626°. The reproducible audit is in `endpoint_geometry.json`
and `measure_endpoint.py`. At this checkpoint, MD/QM/MM had not launched.

## CURRENT — model30 build/minimization launched; direct-PT US stopped, 2026-11-10

Tmux `nyl12_model30_mm_20261110` launched the model30 Amber build and fully
unrestrained MM minimization under `qmmm/pa6_flipped_model30/`. Inputs,
checks, and logs are in `00_system/`, `01_mm_minimization/`,
`source_snapshot/`, and `job/`. The source was distributed Boltz batch
model 30 from aster, not historical ensemble model 30. All four beta Thr
N-termini were neutral; PA6 L4 was zwitterionic with modular RESP/GAFF2
parameters. The heavy-atom mapping check passed (33 atoms, 32 bonds).
Topology acceptance and force/geometry convergence were separate gates.
The planned pilot comprised two neutral and two all-four-NH3+ replicas after
each starting state passed checks.

The direct OG1-to-Thr-N umbrella controller and 20 children on violet were
intentionally stopped to prioritize MD. All 20 windows retained 22 ps of
production data (2 ps pilot plus 20 ps extension), plus interrupted next
blocks. The latest analyzed 12 ps prefix was unconverged; the stop did not
exclude the pathway. See `docs/direct_pt_umbrella_stop_20261110.md`.
The then-current acylenzyme NEB on azalea remained unchanged.

## CURRENT — flipped-pose screen and donor clarification, 2026-11-10

The then-current NEB was active at cycle 329 and unconverged. A screen of
82 existing WT PA6 L4 structures found three positive-phi models passing
Tier A/B plus native-Thr gates. Current E–S phi was −114.341°; the three
alternatives were historical relaxed model 23 and distributed raw models
30 and 25. Raw model 30 was the new lead: phi +72.24°, Tier B,
terminal-N-to-leaving-N 3.77 Å, and psi/chi1 changes +3.61/−0.60°.
It had not yet been minimized or equilibrated. With fixed atom mapping,
opposite signed phi identifies the flipped local pose; the Tier filter used
absolute phi. For an NH3+ microstate, OG1–HG1 remained the proposed donor.
Reports: `docs/pa6_flipped_pose_screen_20261110.md` and
`docs/nyl12_alternative_proton_routes_20261110.md`; source analysis:
`analysis/pa6_binding_orientation_20261110/`.

## CURRENT — repaired native NEB exports; calculation continues, 2026-11-10

IMPORTANT: numbered native NEB NetCDF restart files are INVALID because the
installed Amber custom writer passes x instead of x(lcrd). Do not use them.
Earlier VMD bands and3failed candidate trials were quarantined under
`reaction_paths/03_acylenzyme_neb/clutter/native_restart_export_defect/`.
Internal NEB energies and native `nebpath.xyz` are usable. Ordinary endpoint
minimization restarts use a different writer; acylenzyme release remains valid.

Corrected live monitor: azalea `nyl12_acyl_candidates_xyz_20261110`, immutable
`source_xyz_verified/`. It reconstructs full coordinates from fixed reference
plus mapped native XYZ, checks elements/backbone/report coherence, exports VMD
every10cycles and tests persistent dips (at most3). Independent cycle67 image10
energy exactly reproduced native -706.7081487Eh (0kcal difference at printed
precision). Old candidate monitor stopped. NEB itself continues unchanged.
Current VMD entry remains `vmd_progress/latest/reaction_path.vmd`; copied to
violet and checksum-verified. Updated laptop pull excludes clutter/raw restarts.
Do not draw scientific conclusions from old failed candidate trials.

Provisional native band: HG1 visits terminal N, returns to OG1, then transfers
to leaving N during acylation. Peak remains about40kcal/mol and NEB is not
converged; no viable catalytic barrier or mechanism established. See expanded
repair/proton-route notes in docs/acylenzyme_backward_path_20261110.md.
User accepted Lys156 plus Gly273-backbone QM sensitivity tests; exact masked
fragment/charge design is being documented separately; no change to currentNEB.

## CURRENT — backward acylenzyme endpoint and NEB queue, 2026-11-10

Azalea tmux `nyl12_acylenzyme_neb_20261110` is now running native NEB.
Reaction-unrestrained acylenzyme release CONVERGED and retained chemistry:
OG1-C1.32384Å, C-O1.25576Å, C-Nleave4.11101Å, Nleave-HG1 1.01974Å.
No Hessian/minimum-order proof yet. `status.json` phase is `03_neb`.
Root: `qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths/03_acylenzyme_neb/`.
Method DFTB3/3OB periodic QM/MM; same 98 QM atoms and eight links; same atom
names/indices. Native DL-FIND FREE-end NEB has 16 regular images plus CI slot,
serial evaluations; 4,620 movable atoms selected by fixed ES reference at12Å.
The construction has six temporary chemical restraints; release and NEB remove
all of them. Outer coordinates are fixed, so do not call this globally free.
`nyl12_acyl_neb_monitor_20261110` updates `plots/`;
`nyl12_acyl_candidates_20261110` publishes coherent VMD bands and launches at
most3 persistent-dip image releases. Entries: `plots/neb_profile.svg` and
`vmd_progress/latest/reaction_path.vmd`. Source snapshots and remote action
journal preserved; laptop pull excludes full-system raw coordinates.

See `docs/acylenzyme_backward_path_20261110.md` for scheme correction, procedure,
measured proton-transfer reach, literature and limitations. ThrN-leavingN is
4.937Å in attack09; Asp272/274 nearest O-to-ThrN2.605/2.743Å, but to leavingN
5.988/6.915Å. Both Asp residues are in QM. A2023 NylC study supports homologous
Asp308-Asp306-Thr267 involvement in substrate hydrolysis, not only maturation.
Existing Boltz MSAs give99.76%/99.95% Asp272/274 conservation (2126 unique
>=80%-coverage beta hits), saved in docs/nyl12_active_site_conservation_20261110.md. No fresh BLAST has been submitted. First seed with
large fragment rotation was rejected; deployed local seed has no heavy-heavy
clashes and retains all atoms including leaving amine. Classical acyl MD still
requires new parameters. An endpoint does not prove how its proton moved.

## CURRENT — DFT/MM TI1 release stopped after reactant reversion, 2026-11-10

User authorized stopping the azalea release. Sent Ctrl-C only to tmux
`nyl12_ti1_unbiased_dft_20261110`; verified no sander/g16 processes remained.
Gaussian abort and Sander external-call error are intentional interruption,
not evidence of a spontaneous software failure. All artifacts retained in
`reaction_paths/02_asynchronous_es_to_ti1/04_ti1_unbiased_dft_mm` under the
frame1603 campaign. Last inspected cycle105 geometry had OG1-C 2.9290 A,
C-O 1.2362 A, OG1-HG1 0.9817 A, N-HG1 1.8904 A, N-H1 1.0212 A,
N-H2 1.0271 A: reactant-like reversion, not an accepted TI1 minimum.
The run advanced beyond that inspection before interruption; the retained
restart must be measured again before quoting its final geometry.
Retained restart SHA256 c86656c20ff2c2b2f7d2f98c03cb645be3bdf7f4d323b71c4796207be47895b6.

Interpretation: collapse of this seed under this embedding/active-region setup
is not proof that TI1 is absent or that the entire fast-QM pathway is correct.
Next priority is a bounded diagnostic before another long trajectory: audit
embedding charge continuity and QM boundaries/protonation; analyze existing
DFT reactive forces and oxyanion/attack geometry during release. Then compare
addition-only versus coupled addition/leaving-group protonation/C-N cleavage
with a small, matched set of constrained DFT/MM geometries. The coupled route
is a hypothesis, not an established Nyl12 mechanism. Do not repeat the same
TI1 release or start more umbrella sampling without a defensible chemical path.
Use independent calculations in parallel within resource limits. A candidate
TS must ultimately pass saddle-order and endpoint-connectivity checks; a
coordinate-driving maximum is not sufficient.

## CURRENT — unbiased TI1-like DFT/MM release active, 2026-11-10

User authorized omegaB97X-D/def2-SVP/MM optimization from corrected attack_09.
Session `nyl12_ti1_unbiased_dft_20261110` is active on azalea; Sander and Gaussian
processes were verified. Run directory:
`qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths/02_asynchronous_es_to_ti1/04_ti1_unbiased_dft_mm/`.
Gaussian uses all 28 physical cores and 64 GB memory. Input retains the audited
QM/link definition, 20 A external embedding and 12 A movable shell; atoms outside
the shell remain fixed. Reaction and positional restraints are absent (ntr=0,
nmropt=0). Budget 500 optimization cycles; ntwx=1 requests optimization frames.
This is an active-region minimum test, not a full-system stationary-point proof.
Monitor proton identities, OG1--C, C--O, leaving-group C--N and convergence.
Finite-cutoff membership changes remain an electrostatic sensitivity concern.

In parallel, useful analysis is projected reactive gradients and updated US
overlap/chemistry audits. Missing TI1 does not rule out acylation: addition,
leaving-group protonation and C--N cleavage may form a coupled route to the
acylenzyme without a metastable TI1 basin. Low percent conversion alone cannot
identify the slow step. Redesign must lower barriers relative to populated
states while avoiding over-stabilization of acylenzyme or product complexes.

## CURRENT — coupled four-method panel complete and analyzed, 2026-11-10

All 24 evaluations completed with required termination checks. Six-state
relative-energy RMSE versus omegaB97X-D/def2-SVP is 3.27 kcal/mol for DFTB3,
10.26 for PM6-DH+, and 3.61 for GFN2-xTB. TI1-like versus E--S errors are
+7.42, -14.17 and -1.57 kcal/mol, respectively. DFTB3's error grows by
9.13 kcal/mol from coupled PT to TI1-like; xTB's direct-PT error is +8.02.
No method passes all promotion gates. Prioritize gradient and unbiased DFT/MM
stability tests; xTB merits attack-region follow-up but is not yet validated
for the full pathway. See `docs/coupled_method_panel_results_20261110.md` and
`qmmm/benchmarks/04_coupled_path_method_panel/analysis/method_errors.png`.
Future independent benchmark states should run concurrently within physical-core
budgets. The laptop pull configuration still needs this new result bundle.

## PREVIOUS — six-geometry four-method panel launched, 2026-11-10

The corrected attack branch completed all nine windows with correct terminal
proton identities. Final OG1--C is 1.6526 A and carbonyl C--O is 1.3350 A.
This is a restrained TI1-like seed; unbiased stability remains untested.

The requested matched panel is active on azalea in tmux session
`nyl12_coupled_four_methods_20261110`, under
`qmmm/benchmarks/04_coupled_path_method_panel/`. Six exact geometries are used:
original unrestrained E--S, direct-PT endpoint, asynchronous window 24, and
corrected attack windows 03, 06 and 09. The methods are DFTB3/3OB, PM6-DH+,
GFN2-xTB and omegaB97X-D/def2-SVP. Gaussian uses 16 cores; states run sequentially.
No optimization or driving bias is applied during these single evaluations.
All source restarts and the topology are hashed, input proton identities checked,
and a complete source snapshot retained. Existing builder/parser tests passed
(17 tests); the new orchestration script passed Ruff and remote import checks,
but its coverage has not yet been measured.

The runner writes `work/energies.json` after each complete four-method state.
It does not yet compute final error statistics or validate projected gradients;
those are the next analysis steps. Compare method-relative double differences;
changing finite-cutoff charge membership forbids interpretation of raw profiles.
The earlier plan to freeze charge membership across states is still pending.
Inspect Gaussian normal termination and embedding inventories before accepting
results. The local laptop export has not yet been updated for this panel.

## CURRENT — corrected attack branch and reaction-specific fast-QM selection, 2026-11-10

The gas-phase xTB advantage did not transfer to the actual enzyme QM region.
Consequently, no production Hamiltonian is being chosen from proton affinities
or the isolated five-window direct-PT test. The next selection uses coupled
proton-transfer, carbonyl-polarization, nucleophilic-attack and oxyanion-forming
geometries, with DFTB3/3OB, native Amber PM6 and PM6-DH+, embedded GFN2-xTB and
omegaB97X-D/def2-SVP/MM evaluated on the same sparse structures.

A corrected attack branch now starts from the exact valid `asynchronous_24`
restart (HG1/H1/H2 all N-bound; OG1--C = 2.104 A). It drives only OG1--C to
1.50 A in increments no larger than 0.075 A while flat-bottom gates retain all
three N--H identities. The accepted 12 A movable shell remains free and no
global positional restraints are present. Biased endpoints are benchmark
seeds, not stationary points. Attempt 1 stopped before Sander because the
minimal source snapshot omitted a transitive visualization import; no science
ran. Attempt 2 passed local and remote source-hash/import checks and accepted
three windows. It stopped after converged `attack_04` only because the 0.115 A
target lag exceeded an arbitrary 0.10 A tracking gate; all three N--H identities
were correct. Tracking is now diagnostic. An immutable acceptance revision
promotes that endpoint, and attempt 3 is active in the same calculation
directory on azalea in tmux session
`nyl12_corrected_attack_window24_attempt3_20261110`.

After coupled-method selection, the productive path will use a 2D coordinate:
`xi_PT = d(OG1,HG1)-d(N,HG1)` and `xi_attack = d(OG1,C)`. This permits stepwise
or concerted chemistry without the proton-identity degeneracy of the old scalar
coordinate. Full gates and stationary-point criteria are in
`docs/fast_qm_method_selection_and_coupled_path_plan_20261110.md`.

## CURRENT — embedded enzyme-path DFTB3/GFN2-xTB/DFT comparison, 2026-11-10

The matched five-window Hamiltonian benchmark completed on azalea in tmux
session `nyl12_direct_pt_embedded_methods_attempt1_20261110`. It evaluated the
accepted direct-proton-transfer geometries `pt_01`, `pt_05`, `pt_10`, `pt_14`,
and `pt_19` with DFTB3/3OB, GFN2-xTB 6.7.1, and
omegaB97X-D/def2-SVP. Unlike the earlier gas-phase model comparison, all three
methods used the exact Amber QM/link geometry and the same adjusted, imaged MM
point charges within 20 A for each individual window. DFTB3 was reevaluated
without QM Ewald/PME so the embedding convention matched the external-Gaussian
calculation.

The retained Gaussian input is converted deterministically for xTB. Direct
interface parsing established 98 real QM atoms plus eight link hydrogens (106
total). This corrects the older window-24 manifest's incorrect statement of
six link atoms/104 total. xTB charges use the legacy
point-charge format in bohr with hardness 999.0; a definition-only xTB 6.7.1
smoke test passed before launch.

All 15 engine calculations terminated normally. A parser-only analysis retry
used Sander `EPtot` and completed without rerunning the calculations. DFTB3's
five-window embedded-QM profile agrees with omegaB97X-D/def2-SVP at 1.49
kcal/mol RMSE and -0.57 kcal/mol endpoint error; GFN2-xTB gives 3.91 kcal/mol
RMSE and +7.93 kcal/mol endpoint error. DFTB3 is therefore the closer method in
the actual enzyme QM region despite xTB's advantage in the small gas-phase
models.

The 20-A charge membership varies from 8,088 to 8,107 atoms among windows, so
the large raw finite-cutoff profiles are rejected as reaction profiles. Applying
only the DFTB3-minus-DFT double differences to the original periodic DFTB3/MM
path gives estimated hybrid relative energies of 0.00, 0.42, 7.27, 12.57 and
15.85 kcal/mol. The corrected sparse profile remains monotonic and shows no
transferred-proton minimum. This does not constitute a DFT/MM optimization or
stationary-point proof, but it makes a DFTB3-only energetic artifact an unlikely
explanation for the missing minimum. Prioritize coupled proton transfer and
nucleophilic attack rather than an xTB rebuild for the isolated PT path.

Attempt 1, the parser retry and final correction are under
`qmmm/benchmarks/03_direct_pt_embedded_method_comparison/`. Full rationale,
tables, limitations and acceptance decisions are in
`docs/qmmm_direct_pt_embedded_method_benchmark_20261110.md`.

## CURRENT — direct proton-transfer interim PyMBAR result, 2026-11-10

The chemically valid direct HG1-transfer campaign now has an interim PyMBAR
analysis combining the complete 2-ps pilot and complete 5-ps extension block 1
(7 ps/window; active partial block 2 excluded). A fresh all-window trajectory
audit of block 1 found H1 and H2 N-bound in all 2,000 frames. HG1 changes from
OG1-bound through ambiguous to N-bound across the intended window sequence,
with no proton-identity swap like the asynchronous pathway.

The combined analysis contains 28,000 raw and 16,268 blockwise-decorrelated
samples. Minimum bidirectional adjacent-window overlap is 0.204; the global
PyMBAR overlap scalar is only 0.0127, so sampling is connected but not
converged. Relative to the first umbrella target (xi = -1.090 A), the final
target (xi = 0.750 A) is 15.49 kcal/mol higher. The PMF remains essentially
uphill and supplies no product-side minimum or downhill branch, so 15.49
kcal/mol is the free-energy cost to the restrained endpoint, not an accepted
activation barrier. Adding the 5-ps block changes the common-domain PMF by
0.55 kcal/mol RMS (1.10 kcal/mol maximum) relative to the 2-ps pilot.

Results and hashes are under
`07_direct_pt_umbrella_30ps_violet/analysis/mbar_7ps_complete/result/`; the
laptop pull script now mirrors that exact campaign-relative path. The active
umbrella campaign should continue to 30 ps/window, with cumulative MBAR and
chemistry audits repeated after each complete block.

## CURRENT — intrinsic gas-phase proton-transfer barrier pilot, 2026-11-10

A second-stage method benchmark is active on azalea under
`qmmm/benchmarks/02_small_model_intrinsic_barriers/`. It reuses the three
validated small-model reaction families from the 81-point scan at the 2.7 A
seed separation. Reactant, transition-state and product guesses are optimized
without geometric restraints at omegaB97X-D/def2-SVP; each receives a
frequency calculation and an omegaB97X-D/aug-cc-pVTZ single point. Accepted
DFT stationary-point geometries then receive matched DFTB3/3OB and standalone
GFN2-xTB single points. These approximate-method values are single-point
barrier estimates, not independently optimized approximate-method saddles.

Scientific acceptance requires zero significant imaginary frequencies for
each endpoint, exactly one for each TS, preservation of the declared endpoint
proton identities, normal forward/reverse IRC termination, and opposite proton
identities at the two IRC endpoints. The analysis reports forward/reverse
electronic barriers, reaction energies, DFT ZPE corrections, and errors of the
two approximate methods relative to the DFT reference. Seven focused tests,
shell syntax and Ruff pass; new Python code has 86% focused coverage.

Attempt 1 exited before scientific work because its new Gaussian scratch
subdirectory did not yet exist. It is retained as process provenance. Attempt
2 is active in tmux session `nyl12_intrinsic_barriers_attempt2_20261110` with
nine four-core Gaussian workers. Automatic validation/analysis is queued in
`nyl12_intrinsic_barriers_analysis_queue_attempt2_20261110` and will write to
`qmmm/benchmarks/02_small_model_intrinsic_barriers/analysis/result/` only if
the calculation job exits successfully.

## CURRENT — asynchronous-US chemistry gate and method benchmark, 2026-11-09

The 25-window asynchronous umbrella campaign on azalea was stopped gracefully
after every window completed its first 5-ps production block. A scripted
all-window trajectory audit showed that windows 0--11 retain the reactant
state (HG1 on OG1; H1/H2 on the terminal N), windows 12--15 make the intended
HG1 transfer, and windows 16--20 retain HG1/H1/H2 on N. Window 21 begins an
unintended H1 transfer from N back to OG1; windows 22--24 have H1 fully on OG1
and the OG1--substrate-C distance retreats. Thus the valid prefix is preserved,
but the composite coordinate becomes chemically degenerate after window 20.
Do not extend or interpret the old asynchronous PMF as an E-S-to-TI1 profile.
Partial block-2 artifacts are retained as interrupted provenance. The separate
direct-proton-transfer 30-ps campaign on violet remains active. A matching
all-window audit of its completed 2-ps pilot found H1 and H2 on N in every one
of 800 saved frames; HG1 alone moves from OG1 to N across the intended window
sequence. That campaign therefore passes the current proton-identity gate.

Tao et al. (J. Chem. Phys. 2024, 160, 224104; DOI
10.1063/5.0211276) establishes that Sander can be compiled with direct DFTB+
and xTB library interfaces, including GFN2-xTB and QM/MM Ewald electrostatics.
The retained local AmberTools26 Sander build still lacks xTB linkage; this is a
local build limitation, not an Amber capability limitation. The current
small-model benchmark therefore uses the installed standalone xTB 6.7.1,
Amber SQM DFTB3/3OB, and Gaussian 16 reference single points. Rebuilding
Sander with xTB is deferred until the benchmark justifies using GFN2-xTB in
the enzyme.

The benchmark is under `qmmm/benchmarks/01_small_model_proton_transfer/` and
contains 81 matched fixed geometries: 2-propanol/methylamine (Thr-to-terminal
amine), 2-propanol/acetate (Thr-to-Asp), and acetic acid/methylamine
(AspH-to-terminal amine), each at 2.5, 2.7 and 3.0 A donor--acceptor
separations with nine proton positions. Attempt 1 failed before scientific
work because AMBERHOME was absent. Attempt 2 completed all DFTB3 and GFN2-xTB
points, but Gaussian rejected `def2TZVPD` as an ambiguous built-in keyword.
Attempt 3 preserves those attempts and uses omegaB97X-D/aug-cc-pVTZ, a
recognized diffuse triple-zeta reference for H/C/N/O. All 81 calculations
from each of DFTB3, GFN2-xTB and Gaussian completed. The calculation wrapper
then exited nonzero because its acceptance parser retained the obsolete
`wb97xd_def2tzvpd` dictionary key; no scientific calculation failed. The key
was corrected, shell syntax/Ruff/three focused tests passed, and the completed
outputs were analyzed without rerunning them.

The analysis is complete under
`qmmm/benchmarks/01_small_model_proton_transfer/analysis/result/`. Across the
nine fixed-geometry profiles, mean per-profile MUE/RMSE values are 3.29/4.71
kcal/mol for DFTB3 and 1.89/2.31 kcal/mol for GFN2-xTB; worst absolute errors
are 17.74 and 8.43 kcal/mol, respectively. GFN2-xTB is closer to the
omegaB97X-D reference for each of the three reaction families in this pilot.
This supports testing GFN2-xTB/MM in the enzyme, but these scans are
fixed-geometry potential profiles, not activation free energies, optimized
proton affinities, or validation of an enzyme reaction path.

The direct-proton-transfer umbrella campaign on violet remains active in
production block 2 across 20 windows. Completion of this block gives about 10
ps/window of production sampling; the chemically contaminated asynchronous
campaign remains stopped.

DeePMD-kit is available through the published Amber infrastructure, but no
Nyl12 delta-MLP is being trained now. A correction would require a curated
DFT QM/MM energy/force training set, uncertainty/transferability tests, and
independent path validation. That is a possible later acceleration layer, not
the fastest route to establish whether DFTB3 or GFN2-xTB is adequate.

## CURRENT — direct-PT umbrella and asynchronous E-S-to-TI1 outcome, 2026-11-09

The original asynchronous construction path stopped at OG1--carbonyl-C =
2.104 A. Its bounded eight-window continuation completed through
`asynchronous_32`. Provenance confirms that window 25 started directly from
window 24's optimized restart. Nevertheless, the system made a sharp downhill
structural jump at window 25 and the continued composite coordinate was
satisfied mainly by proton identity exchange rather than C--O bond formation.
Specifically, Thr929 H1 transferred from the N-terminal N back to OG1 while
the originally tracked HG1 remained on N: N--H1/OG1--H1 changed from
1.021/3.103 A in window 24 to 2.074/0.977 A in window 25. The scalar coordinate
tracked only HG1 and therefore could not distinguish this chemically
reactant-like proton swap from progress toward TI1.
The final window has OG1--C = 2.413 A, N--HG1 = 1.026 A, OG1--HG1 = 3.749 A,
and carbonyl C--O = 1.261 A, so it fails the declared TI1-like gate. No
unrestrained TI1 release or new extension umbrella windows were launched.

At the user's request, continuation windows 25--32 now reside directly beside
windows 01--24 in the canonical `01_coupled_coordinate_drive/` directory. The
obsolete empty extension directories were removed. The canonical
`plots/potential_energy.{json,csv,png}` was regenerated in place for all 32
windows; the prior 24-window plot and the temporary duplicate plot bundle are
preserved under `clutter/asynchronous_path_preconsolidation_20261109/`.

Do not extend the old scalar composite coordinate further. The next TI1
construction should restart from window 24, where N has H1, H2 and transferred
HG1 and OG1 has no proton. Drive OG1--carbonyl-C independently toward a
covalent distance while applying explicit flat-bottom retention gates to all
three N--H bonds and exclusion gates preventing H1/H2 from bonding to OG1.
All other protein/MM degrees of freedom remain unrestrained. Release the
proton-identity and attack restraints in a staged and then fully unrestrained
optimization. Failure to retain TI1 after that corrected construction is a
Hamiltonian/protonation/mechanism result, not evidence that a one-dimensional
composite target was reached.

The VMD display now removes static topology bonds for changing connectivity
and uses separate DynamicBonds cutoffs: 1.40 A for OG1--HG1 and N--HG1, and
1.90 A for OG1--carbonyl-C. Only non-QM substrate carbons are overlaid in
orange; other substrate atoms use element colors. Visible measurements include
OG1--HG1, N--HG1, OG1--C, carbonyl-O--Phe113(N),
carbonyl-O--Asn186(ND2), and the OG1--C--O Burgi--Dunitz angle. The prior
25-frame movie ended with the *converged endpoint* of the unrestrained release
from `asynchronous_24`, which reverted to E-S. During the continuation, that
failed release is omitted so the movie remains a continuous driven path.

The 19-window direct Thr233-HG1-to-N-terminal-N coordinate drive completed.
Its reported Sander potential energy increased monotonically to 15.38 kcal/mol
relative to the first window (15.28 kcal/mol after subtracting the reported
reaction-coordinate restraint). The unbiased endpoint optimization returned
HG1 to OG1, so no transferred-proton stationary point was established.

The 24-window asynchronous composite-coordinate drive also completed. Its
biased endpoint reached OG1-HG1/N-HG1/OG1-C distances of
2.275/1.017/2.104 Å, but unbiased optimization converged back to E-S at
0.977/2.073/3.239 Å. The potential energy rose monotonically to 21.01 kcal/mol
(20.81 kcal/mol with the reported restraint removed), with no downhill TI1
basin. This path is a construction scan and umbrella seed, not a barrier.

Reproducible JSON, CSV and PNG profiles are under each run's `plots/`
directory. The laptop pull script includes both checksum-verified bundles.

The corresponding 20-window direct-proton-transfer umbrella pilot completed
successfully on 2026-11-09. All windows finished 1 ps heating, 1 ps
equilibration and 2 ps production. The preliminary 2-ps histograms show strong
adjacent-window coverage: histogram-intersection overlap ranges from 0.440 to
0.705 (mean 0.558), with no adjacent pair below 0.2. These are encouraging
sampling-connectivity diagnostics, not final PMF convergence evidence. The
checksum-verified plot and numerical summaries are under
`04_direct_pt_umbrella_initial/analysis/03_production_2ps_overlap/`.

The 30-ps extension launched automatically in tmux session
`nyl12_direct_umbrella_30ps_queue_20261109`; all 20 windows are currently in
the additional 5-ps equilibration stage, before six retained 5-ps production
blocks. A general tested exporter now orders windows by reaction-coordinate
target, samples evenly across every completed production block, caps retained
frames independently per window, concatenates them, and writes a frame-level
source manifest plus the established reactive-site VMD display. The completed
pilot bundle contains all 40 available frames from each window (800 total,
19,967 stripped atoms) at
`04_direct_pt_umbrella_initial/visualization/sequential_production_windows/`.
Its `view_umbrella.vmd` uses the same dynamic bonds, element/orange substrate
coloring, oxyanion-hole distances and Burgi--Dunitz angle as the coordinate-
drive movies. The laptop pull script mirrors both the overlap and VMD bundles.

A 25-window asynchronous-coordinate QM/MM umbrella pilot was launched at
14:46 EST in tmux session `nyl12_asynchronous_umbrella_20261109`. It spans
xi = d(OG1-HG1) - d(N-HG1) - d(OG1-C) from -4.317 to -0.800 Å, uses the E-S
structure plus all 24 driven endpoints, and runs 1 ps heating, 1 ps
equilibration and 2 ps production per window at 300 K with a 0.5-fs step and
50 kcal mol-1 Å-2 bias. No positional restraints are present. All 25 seed
restart hashes are recorded in the campaign input manifest. Violet is running
20 direct-PT plus 25 asynchronous Sander processes without CPU oversubscription.
The initial 2-ps production is only an overlap/equilibration pilot; PMF
interpretation requires histogram overlap, autocorrelation/effective-sample
size checks and likely longer production.

The estimator and extension protocol is now specified in
`docs/qmmm_umbrella_sampling_and_pmf_protocol_20261109.md`: MBAR is primary,
WHAM is an independent cross-check, and accepted pilots are extended in 5-ps
blocks to at least 30 ps of production per window, excluding discarded
equilibration. Longer sampling remains mandatory if convergence, overlap,
uncertainty or hidden-coordinate gates fail. `pymbar==4.0.3` is installed and
pinned through `uv`; optional JAX acceleration was not added.

The production campaigns are split cleanly across hosts. Tmux controller
`nyl12_direct_umbrella_30ps_queue_20261109` waits for the direct-PT pilot and
then runs 5 ps equilibration plus six retained 5-ps production blocks on
violet. The short asynchronous violet attempt had not completed even its
heating stage; it and its obsolete dispatcher were stopped and moved intact to
`clutter/obsolete_violet_asynchronous_pilot_20261109/`. A single unified
25-window campaign was launched from the original coordinate-drive seeds on
azalea as `nyl12_async_umbrella_30ps_20261109`: 1 ps heating, 5 ps
equilibration and six retained 5-ps production blocks. Its canonical output is
`06_asynchronous_umbrella_30ps_azalea/`. Azalea has 56 CPUs and 3.2 TB free.
Both hosts use their native AmberTools 26 installations; different executable
hashes are expected because the hosts have different architectures, while the
topology hash is identical. No positional restraints are present.

Completed reaction-path artifacts now use the same canonical relative layout
on the workstations and in the laptop bundle:
`reaction_paths/01_proton_transfer/` and
`reaction_paths/02_asynchronous_es_to_ti1/`. Each PES plot is named
`plots/potential_energy.png`. Historical JSON and engine logs retain their
original absolute path strings as immutable provenance.

Both completed coordinate-drive paths now have checksum-verified stripped VMD
bundles under their matching `vmd_progress/latest/` directories. The direct
proton-transfer bundle contains 20 frames (19 driven endpoints plus the
unrestrained release), and the asynchronous bundle contains 25 frames (24
driven endpoints plus its unrestrained release). The laptop pull script mirrors
both bundles and their PES `plots/` directories under the same relative tree.

## CURRENT — separate polyamide Amber FF project proposed, 2026-10-28

A separate `polyamide-amber-ff` project is recommended for modular PA6/PA66
residues, termini/caps, Boltz/AF3 mapping and enzyme-state extensions. The
Nyl12 repository should consume versioned releases rather than own the general
library. Prior nylon force fields mean novelty must come from the enzyme-ready
state library, reproducible mappings and broad QM/MD validation, not merely
another PA6 charge set. E--S can use the ordinary modular substrate; TI1 and
acyl-enzyme fitting remain gated on validated stationary-point connectivity,
proton ownership and charge. The proposed scope, paper positioning, state
matrix, repository layout and compute allocation are documented in
`docs/polyamide_amber_ff_project_charter.md`. Azalea is reachable and idle with
56 CPU cores, 503 GiB RAM, two 98-GB GPUs, Gaussian access and AmberTools
Sander. Tulip SSH from violet currently fails public-key authentication.

## CURRENT — modular PA6 residue-library proposal, 2026-10-28

The user proposed representing each polyamide monomer as an Amber residue so
arbitrary oligomer lengths can be assembled. The concrete PA6 boundary is the
inter-repeat amide C--N bond: each internal residue contains
`NH-(CH2)5-C(=O)`, connects through head N and tail carbonyl C, and has net
charge zero. Separate charged and neutral terminal templates plus nylon-derived
acetyl-like and N-methylamide-like caps are required. Charges must be fit from
capped/terminal oligomers with RESP equivalence and integer-residue constraints,
not copied from an isolated monomer or directly partitioned from the current
full-L4 fit. Boltz/AF3 conversion must use heavy-atom graph isomorphism and an
immutable atom-map manifest rather than predicted atom names/order. The full
plan and acceptance gates are in
`docs/modular_polyamide_residue_parameterization_plan.md`. No modular
production parameter has been certified, and current QM scans continue.

## CURRENT — PA6 L4 QM bonded-term benchmark active, 2026-10-28

The four-conformer RESP fit, GAFF2 typing, explicit-solvent minimization and
170-ps staged MD sanity test are complete. They establish a usable candidate,
not a certified production parameter set. The remaining gate is a direct QM
benchmark of the transferred torsions and `parmchk2` impropers.

Three small, charge-explicit nylon-6 fragments were optimized normally at
omegaB97X-D/def2-SVP/SMD(water): neutral internal
`CC(=O)NCCCCCC(=O)NC`, N-terminal `[NH3+]CCCCCC(=O)NC`, and C-terminal
`CC(=O)NCCCCCC(=O)[O-]`. None used `CalcFC`; all expected bonds passed the
endpoint gate.

The completed controller `pa6_torsion_scan_controller_20261028` wrote
`05_qm_torsion_benchmarks/02_relaxed_scans/` and launched eight four-core scan
sessions. They cover four internal proper torsions, terminal ammonium rotation,
both flagged amide impropers, and the carboxylate improper. The positive
direction is active in all eight sessions; the negative direction follows
sequentially, limiting demand to 32 cores. Proper scans cover full rotations;
improper scans are complementary +/-60-degree branches from the same minimum,
not full hysteresis cycles. Acceptance requires complete optimized windows,
geometry continuity and agreement wherever angular profiles overlap. Accepted
geometries then receive def2-TZVP
single-point energies and fragment RESP/GAFF2 comparison. No new parameter is
fit automatically, and the full PA6 L4 candidate remains uncertified.

## CURRENT — authoritative plan and Model 23 late-cycle diagnosis, 2026-10-27

The controlling scientific plan is now
`docs/MASTER_PLAN_NYL12_COVALENT_REDESIGN.md`. It defines locked facts,
required chemical states, measurable G0--G7 gates, DFTB3-versus-DFT
validation, first-order-saddle acceptance, multistate symmetric design, data
organization, and explicit pivot conditions. Topic documents remain detailed
evidence; this handoff remains the operational checkpoint.

The Model 23 whole-system DFTB3/MM optimization remains active in tmux session
`nyl12_model23_dftb_es_20261027`. It uses one CPU core, not a GPU; both GPUs
were free when checked. The user requested that it continue as long as it does
not block needed resources. A provenance-preserved cycle-6790 snapshot is at
`qmmm/pa6_model23/02_dftb3_es_optimization/diagnostics/cycle_6790_snapshot/`.
Its PA6 metrics are d1 3.47397 A, d2 3.08211 A, d3 3.03706 A, Burgi-Dunitz
98.8186 degrees, projected Herschlag phi 124.8201 degrees, and
CB--OG1--C54 108.7603 degrees. It is only Tier B and lies extremely near the
d1/phi Tier-B limits.

After fitting the n12c catalytic scaffold, the same snapshot has Thr233
heavy-atom RMSD 1.47992 A, N/OG1 displacements 2.31548/1.91531 A, chi1 change
+25.3189 degrees, and psi change -57.9843 degrees. It is degrading relative to
the accepted starting structure and must be retained as a negative-control
basin, not promoted merely if minimization converges.

The next production branch is a bounded n12c-template-conditioned E-S pilot.
Reuse the validated PA6 L4-to-L2 ligand graph and two existing generation
contacts, preserve the mature tetramer and canonical mapping, generate 12
independent Boltz seeds, then release contacts during staged relaxation and
apply G1 to the complete denominator. The installed custom Boltz2 parser does
support explicit CIF/PDB templates and chain mapping. Boltz reads the four
precursor chains by mmCIF label-asym IDs A--D (corresponding to deposited
auth-asym IDs A/F/H/V). The split-query mapping is therefore A/E->A, B/F->B,
C/G->C and D/H->D. Duplicate template IDs still require an end-to-end parser
preflight before launch. AF3 remains a second arm after the user's separate
contact-restraint implementation is validated.

The corrected mapping passed an installed-parser preflight: 9 query chains,
8 template records and 1 parsed template. The 12-seed Boltz pilot completed
successfully on 2026-10-27 at 19:02:32 UTC with exit status 0 and 12/12 complete
canonical PDB/confidence artifact sets. Inputs and hashes are under
`inputs/template_conditioned_es_pilot/`; the provenance-recorded run is
`models/nyl12_template_conditioned_es/01_boltz_generation/`.

Process success did not produce a G1 scientific promotion. The established
four-site analyzer found the intended PA6 amide and AE site in all 12 models,
but the best geometry per model comprised 0 Tier A, 3 Tier B (models 1, 7 and
8), and 9 rejects. Model 7 is the closest miss: d1/d2/d3 =
2.603/2.887/2.818 A, Burgi-Dunitz = 89.837 degrees and projected phi =
-115.312 degrees. Model 8 has d1/d2/d3 = 2.532/2.800/3.002 A,
Burgi-Dunitz = 87.409 degrees and projected phi = -106.455 degrees. Model 1
has d1/d2/d3 = 3.021/3.115/2.492 A, Burgi-Dunitz = 86.412 degrees and
projected phi = -123.635 degrees. Thus proximity is often good but the attack
orientation is not Tier A. The other nine models fail more strongly.

Protein/assembly confidence is uniformly high (confidence score mean 0.936,
complex pLDDT mean 0.927 and ipTM mean 0.972), but ligand coordinates are much
less certain: ligand atom B-factor/pLDDT means are only 38.3--43.6 on the
0--100 PDB scale despite ligand ipTM values of 0.820--0.898. Confidence does
not rescue the geometry failure. Raw analysis is under
`01_boltz_generation/analysis/raw_geometry/`. Per the predeclared G1 gate,
none proceeds to production relaxation as an accepted seed. Models 7 and 8
may be retained as explicitly labelled near-miss diagnostics when refining the
generation strategy; they are not accepted E-S structures.

## CURRENT — hASNaseIII precedent and replacement E-S seed screen, 2026-10-26

Andjelkovic et al. (JCIM 2023, DOI 10.1021/acs.jcim.3c00900) provide a directly
relevant N-terminal-nucleophile substrate-hydrolysis precedent, distinct from
precursor autocatalytic cleavage. Their selected acylation path first transfers
the catalytic Thr hydroxyl proton directly to the same residue's neutral
N-terminal amine (TS1, 14.0 +/- 0.4 kcal/mol), then performs nucleophilic attack
(TS2, 18.6 +/- 0.7 kcal/mol). They screened alternatives with DFTB3/MM and
refined the selected path with B3LYP-D3/6-31+G*/MM adaptive-string free-energy
calculations. This supports testing a direct Thr233 OG1-to-N transfer before
imposing an Asp272 proton shuttle. It does not establish the Nyl12 mechanism or
transfer hASNaseIII barriers quantitatively. Reference R41 and manuscript notes
preserve the source and mechanistic distinction.

A catalytic-scaffold alignment to the supplied `cif/n12c.cif` confirmed that
the current solvated E-S seed has a suspicious Thr233 backbone/N-terminal pose.
After fitting 16 N/CA/C/O atoms from Phe113, Asn186, Asp272 and Asp274, the
converged E-S Thr233 heavy-atom RMSD is 1.832 A; N and OG1 are displaced 2.730
and 2.257 A, and psi differs by -45.64 degrees. The selected raw MD frame already
showed the defect (Thr RMSD 1.616 A, N/OG1 shifts 2.284/2.146 A, psi difference
-40.85 degrees), so the current minimization did not create it. The raw frame
also shortens OG1--Asp272 OD1 from 4.920 A in the reference to 2.950 A, plausibly
biasing the earlier Asp-relay hypothesis.

The user's proposed replacement strategy is now preferred: use an
OpenMM-relaxed, canonical-chain static Boltz E-S model that first passes the
existing PA6 L4-to-L2 Tier A gate and then passes a separate native-like Thr233
screen. Only WT models 23 and 40 are Tier A. Model 23 is the better catalytic
scaffold seed: fit RMSD 0.568 A, Thr233 heavy-atom RMSD 0.623 A, N/OG1 shifts
0.705/0.188 A, chi1 difference +4.52 degrees and psi difference -30.09 degrees.
Model 40 has a smaller scaffold fit RMSD (0.346 A) but a worse Thr233 pose
(heavy-atom RMSD 1.301 A, N shift 2.149 A, psi difference -56.03 degrees).
Artifacts are under `qmmm/pa6_static_es_seed/`; model 23 is the current leading
seed, subject to visual inspection and topology/protonation rebuilding before
QM/MM. Tier A plus structural similarity is a seed gate, not evidence of a TS.

Parallel DFT/MM was audited but not launched from the defective MD seed. The
installed QUICK backend is CPU-serial (`MPI=FALSE`, `CUDA=FALSE`, OpenMP off).
Amber Sander's external Gaussian interface supports electrostatic embedding,
analytic Gaussian forces and `%NProcShared` through `&gau num_threads`; Gaussian
16 is therefore the available parallel backend. With `qm_ewald=1`, a future
pilot must explicitly use compatible Mulliken charge analysis. Launch only
after the replacement E-S topology, protonation and direct-transfer coordinate
are verified.

## CURRENT — driven path complete; released TI1 reverted, 2026-10-26

All 42 driven windows (`stage1_00`--`stage1_20` and
`stage2_00`--`stage2_20`) converged and passed their declared window gates.
The final driven endpoint is TI1-like: OG1-C 1.70283 A, Thr N-HG1 1.04794 A,
OG1-HG1 2.58677 A, and C-O 1.33213 A. The controller selected this exact
`stage2_20/optimized.ncrst` as the starting geometry for fully unrestrained TI1
validation; its SHA-256 is identical to
`02_ti1_unrestrained_validation/starting.ncrst`
(`956c5667b0cf6f85c2fdd1452e408f26a260f52362279c00555e2c7dfa5b662b`).

The unrestrained DFTB3/MM minimization ran inside the provenance-recorded tmux
controller and converged natively after 416 DL-FIND cycles. It had `ntr=0`,
`nmropt=0`, no DISANG or positional restraints, all atoms active, and
`tole=1e-5` Hartree. It did **not** retain TI1: OG1-C relaxed to 3.05524 A,
HG1 returned to OG1 (0.99217 A), Thr N-HG1 lengthened to 2.39747 A, and the
substrate C-O/C-N distances returned to 1.26529/1.35661 A. The exact acceptance
record is `02_ti1_unrestrained_validation/acceptance.json`, with
`ti1_geometry_retained=false`. Therefore the current DFTB3/MM Hamiltonian and
environment do not establish an unrestrained TI1 minimum, even though the
biased construction profile descends after `stage2_09`. Do not launch an
unchanged duplicate minimization.

A verified stripped VMD trajectory of the released optimization is at
`02_ti1_unrestrained_validation/vmd_optimizer_trajectory/`; load
`view_stage2.vmd`. It contains 16 saved optimizer frames (written every 25
cycles), 19,976 atoms, all 98 physical QM atoms, the complete protein/substrate
display selection, and three complete nearby waters. The DCD coordinates and
fixed boxes were independently verified exactly against the NetCDF source.
These are optimizer snapshots, not MD time; the exact converged endpoint remains
`02_ti1_unrestrained_validation/optimized.ncrst`. The laptop pull script includes
and verifies the bundle.

The energy discontinuity between `stage2_09` and `stage2_10` coincides with
exchange of the terminal Thr ammonium hydrogen donated to Asp970 OD2. At
`stage2_09`, H2--OD2 is 2.03513 A and H1--OD2 is 2.86346 A; at `stage2_10`,
H1--OD2 is 2.04305 A and H2--OD2 is 3.48213 A. This supports the user's visual
interpretation of the discontinuity as a hydrogen-bond rotamer switch rather
than a distinct chemical elementary step.

## CURRENT — stage2_02 reconciled; path continuing to TI1, 2026-10-26

The prior `stage2_02` run exhausted 2000 DL-FIND cycles because its
`tole=1e-7` Hartree energy-change threshold was below the observed DFTB3/MM
numerical fluctuations. The RMS gradient and step criteria were already
satisfied, and the largest residual maximum-component motions were remote MM
atoms rather than relay atoms. The retained failed-attempt artifacts remain as
`stage2_02/*_attempt1`. An in-place continuation from the exact final restart,
with no change to the Hamiltonian or three reaction-coordinate biases and with
`tole=1e-5` Hartree, converged in eight additional cycles. Its accepted endpoint
has OG1-C 1.75170 A, OG1-HG1 1.75522 A, Asp272 OD1-HG1 1.02165 A, and Thr233
N-HG1 2.34463 A. Native convergence, fixed-cell coordinates, and proton
ownership all pass. Acceptance and process provenance are in `stage2_02/`.

The strict reconciliation script `stage2_02/reconcile_stage2_02.py` required
the exact stopped state through `stage2_01`, revalidated the resumed endpoint,
preserved the prior status as
`status_before_stage2_02_reconciliation_20261026.json`, and appended only
`stage2_02`. Its immutable record is `stage2_02/reconciliation.json`; the next
window is `stage2_03`.

Active sequential controller:
`nyl12_coordinate_drive_stage2_03_to_ti1_20261026`; provenance:
`01_coordinate_drive_to_ti1/controller_resume_after_stage2_02_20261026/`.
It skips all accepted windows, runs `stage2_03` through `stage2_20`, and then
performs the fully unrestrained TI1 minimization in
`02_ti1_unrestrained_validation/`. It fails closed on engine, native
convergence, coordinate, proton-ownership, or final TI1 checks. At this
checkpoint `stage2_03` has converged and passed all window gates, `stage2_04`
is active, and no TI1 endpoint is yet claimed.

Active read-only inspection publisher:
`nyl12_path_inspection_stage2_03_to_ti1_20261026`; provenance:
`01_coordinate_drive_to_ti1/inspection_monitor_after_stage2_02_20261026_attempt2/`.
It updates the canonical VMD movie and energy/distance plots after accepted
windows. The first replacement publisher attempt failed without affecting
Sander because the reconciled window lacked canonical `process.json`; that
metadata was regenerated from the preserved attempt and resume records, a
one-shot export passed, and the replacement monitor is healthy.

The `tole=1e-5` generator change and its regression assertion are uncommitted.
Focused relay tests pass (10/10), Ruff passes, and no production dependency was
added.

## CURRENT — retained-coordinate phase II active, 2026-10-26

The first attempted repair of phase II also failed scientifically. A one-sided
1.35-A lower wall on OG1-HG1 allowed a converged shared-proton structure at
`stage2_00`: OG1-HG1 1.22198 A, Asp272 OD1-HG1 1.28484 A, and Thr N-HG1
2.79899 A. DL-FIND converged normally in 393 cycles, so the failure was the
reaction-coordinate model, not the optimizer. The hard proton-ownership gate
correctly rejected it. Its complete artifacts are preserved under
`clutter/failed_phase2_soft_wall_20261026/`; do not treat that endpoint as an
accepted reaction-path window.

Attempt 5 replaces the soft wall with a retained first-transfer coordinate:
`d(OG1-HG1)-d(Asp272_OD1-HG1) = 0.76217766 A`, taken directly from accepted
`stage1_20`. Phase II independently drives
`d(Asp272_OD1-HG1)-d(Thr233_N-HG1)` while retaining the 1.55-A attack target.
There are no MM positional restraints. The exact specification is
`01_coordinate_drive_to_ti1/phase2_retained_coordinate_20261026.json`.

Active calculation session:
`nyl12_coordinate_drive_phase2_retained_20261026`; provenance:
`01_coordinate_drive_to_ti1/resume_job_20261026_attempt5/`. Corrected
`stage2_00` and `stage2_01` both converged and passed proton ownership. At
`stage2_01`, OG1-HG1 is 1.74619 A, Asp272 OD1-HG1 is 1.01772 A, and Thr
N-HG1 is 2.44905 A; `stage2_02` is active. The read-only monitor
`nyl12_retained_path_inspection_20261026` publishes the accepted path after
each window. The live laptop-ready path is
`01_coordinate_drive_to_ti1/vmd_progress/reaction_path.vmd` (23 accepted
endpoint frames through `stage2_01` at this checkpoint).

The failed soft-wall optimizer trajectory is separately inspectable at
`01_coordinate_drive_to_ti1/failed_soft_wall_stage2_00_vmd/view_stage2.vmd`.
It contains 15 saved optimizer frames (the run used `ntwx=25`, so it does not
contain all 393 optimization cycles), a stripped 19,976-atom topology, all 98
QM atoms, and three complete nearby waters. The exporter independently checks
all DCD coordinates and boxes. Its VMD display includes hidden protein cartoon,
orange non-QM substrate carbons, QM licorice, an HG1 sphere, DynamicBonds,
reactive and oxyanion-hole distances, and the OG1-C54-O49 Burgi-Dunitz angle.
The laptop pull script includes this directory and verifies nested manifests.

## CURRENT — phase-II proton-return defect corrected, 2026-10-26

The first phase-II implementation was scientifically invalid. Its sole proton
coordinate, `d(Asp272_OD1-HG1)-d(Thr233_N-HG1)`, could reach its target after
HG1 returned covalently to Thr233 OG1 and Asp272 moved away. At accepted
`stage2_05`, OG1-HG1 was 1.00503 A, Asp272 OD1-HG1 was 1.77026 A, and Thr
N-HG1 was 2.70623 A. The calculation was interrupted during `stage2_06`; no
sander process from that attempt remains. This branch must not be interpreted
as an Asp-to-N proton-transfer path.

All seven invalid `stage2_*` directories were moved without deletion to
`clutter/invalid_phase2_hg1_return_20261026/01_coordinate_drive_to_ti1_stage2_archive/`.
Its `manifest.json` records source/destination mappings and before/after file
hashes; `status_before.json` permits rollback. The live run directory was reset
to the accepted `stage1_00` through `stage1_20` sequence. The repair utility is
`scripts/archive_invalid_phase2.py`, and its provenance job is
`01_coordinate_drive_to_ti1/phase2_repair_job_20261026/`.

Corrected phase II retains the attack and Asp-to-N difference coordinates and
adds a one-sided lower wall at 1.35 A on Thr233 OG1-HG1. This wall prevents
covalent return to OG1 but applies no bias above the lower bound. It is a
reaction-specific internal-coordinate bias, not a positional restraint on MM
atoms. Every phase-II window now fails acceptance if OG1-HG1 is more than
0.10 A below the declared wall. The exact rationale and settings are in
`01_coordinate_drive_to_ti1/phase2_correction_20261026.json`.

Active calculation session:
`nyl12_coordinate_drive_phase2_corrected_20261026`; provenance:
`01_coordinate_drive_to_ti1/resume_job_20261026_attempt4/`. It started from the
original accepted E-S endpoint, reused rather than recomputed accepted stage I,
and began regenerated `stage2_00` from `stage1_20/optimized.ncrst`. Active
read-only inspection session: `nyl12_corrected_path_inspection_20261026`;
provenance: `inspection_monitor_job_20261026_attempt4/`. The monitor updates
`vmd_progress/` after each accepted window.

## CURRENT — converged E-S and unrestrained coordinate path, 2026-10-26

The initial 5000-cycle unrestrained E-S minimization exhausted without native
DL-FIND convergence but retained a valid E-S geometry. The conditional
continuation started from its exact restart and converged after 1121 additional
cycles. Acceptance: native convergence true, fixed finite cell true, no
restraints, OG1-C 3.16529 A, OG1-HG1 0.99215 A, Asp272 OD1-HG1 1.82505 A;
accepted endpoint SHA-256
`93a7f1a4b2508c6ed7bf94b1ec184371da0ae54c1bba4b085584dc6e53d2c80f`.

The first automatic path launch was stopped before a window because the live
publisher consumed the prior terminal status before the resume child changed
it. A tested state-transition guard fixed that race. The next attempt exposed
two non-scientific input/publishing defects: the Amber `&dlfind` namelist must
follow `&wt`/`DISANG` declarations, and a configured new `vmd_progress/`
destination needs an initialized manifest. Both now have regression tests.
Failed zero-step/interrupted artifacts were moved reversibly to
`clutter/failed_coordinate_drive_resume_20261026/`.

Current calculation tmux session: `nyl12_coordinate_drive_resume3_20261026`.
The path uses
the accepted E-S endpoint, no positional restraints, and only the two declared
reaction-coordinate biases. Windows `stage1_00` through `stage1_16` converged
and passed acceptance; `stage1_17` is active. At `stage1_15`, OG1-C is 2.10955
A, OG1-HG1 is 1.42309 A, and Asp272 OD1-HG1 is 1.09954 A. No adjacent-window
distance change exceeds the declared 0.25 A discontinuity threshold. Relative
total energy has nevertheless risen monotonically to 35.27 kcal/mol (32.38
kcal/mol after subtracting the declared driving bias), so this remains a
construction scan and is not a validated barrier. The controller resumes from
accepted windows without recomputing them. Live state:
`qmmm/pa6_asp272/01_coordinate_drive_to_ti1/status.json`; active job provenance:
`qmmm/pa6_asp272/01_coordinate_drive_to_ti1/resume_job_20261026_attempt3/`.

The original live publisher omitted three VMD dependencies and had no live
energy plot. The corrected publisher now atomically includes
`stripped.parm7.gz`, `qm_indices.tcl`, `nyl12_es_vmd_common.tcl`, and the
current `energies.*` and `distances.*` files. All focused publisher, inspection,
energy, and package tests pass (42 tests). Because the active controller has a
frozen source snapshot, a separate read-only tmux monitor,
`nyl12_complete_path_inspection_20261026`, updates the corrected bundle after
each accepted window without altering sander. Its provenance is under
`01_coordinate_drive_to_ti1/inspection_monitor_job_20261026/`; the laptop-ready
bundle is `01_coordinate_drive_to_ti1/vmd_progress/`.

The complete construction schedule is 21 stage-I windows (`stage1_00` through
`stage1_20`), then 21 stage-II windows (`stage2_00` through `stage2_20`), then
a fully unrestrained TI1 minimization. Thus, with `stage1_17` active, 25 driven
windows remain including the active one. The workflow does not stop after the
first proton transfer. Stage II transfers Asp272-H toward the terminal Thr N
while retaining the bounded 1.55-A attack target; the final release tests
whether TI1 persists without either driving bias. If the energy has not turned
downhill or TI1 is not retained after that release, adding extrapolated windows
with shorter bond targets would not be justified; the reaction coordinate or
mechanistic hypothesis must instead be revised.

The live VMD publisher now also installs stage-aware native bond-distance
labels, DynamicBonds representations for the forming/breaking bonds, persistent
O49--Phe113 N and O49--Asn186 ND2 distances, and a small HG1 sphere. The helper
is `scripts/nyl12_reactive_distances.tcl`; the replacement read-only monitor is
`nyl12_vmd_labels_path_inspection_20261026`. It waits 25 seconds after each status
change so its complete annotated bundle is published after the frozen primary
publisher, without modifying or interrupting sander.

Native bond-distance labels are explicitly shown on load rather than relying on
VMD's default label visibility. The two persistent oxyanion-hole measurements,
O49--Phe113 backbone N and O49--Asn186 ND2, are also explicitly shown. A native
angle label reports the Burgi-Dunitz angle as OG1--C54--O49, with C54 as the
central atom. The current monitor provenance is
`01_coordinate_drive_to_ti1/inspection_monitor_job_20261026_attempt3/`.

## CURRENT — top-level output cleanup, 2026-10-25

Fifty superseded top-level output directories were moved reversibly into
`clutter/superseded_exports_and_test_runs_20261025/`: 37 numbered laptop
inspection snapshots, one pre-review snapshot, nine temporary test-job output
directories, and three failed VMD packaging attempts. No file was deleted.
The move verified per-file SHA-256 hashes and totaled 654196829 bytes; see
`provenance/output_directory_cleanup_20261025/move_manifest.json` and its logged
tmux job. Canonical `output_laptop_inspection_20261021/` and
`output_vmd_inspection_20261024/` remain because the laptop pull script uses
them. Action logs, parameter cache, active `qmmm/` calculations, and scientific
`output_nyl12_*` evidence remain in place. Restore any moved directory by using
the manifest's source/destination mapping.

## CURRENT — corrected unrestrained E-S to TI1 workflow, 2026-10-25

User rejected the prior blanket k10 positional restraints on all 103153 MM
atoms. The prior restrained paths remain diagnostic history only. Their bulky,
regenerable XYZ/NetCDF/live-publisher artifacts (283 files; 112478870857 bytes)
and curated inspection bundle were moved reversibly to `clutter/`; see
`provenance/qmmm_cleanup_20261025/cleanup_manifest.json`. Core restarts, mdout,
inputs, acceptance records, source snapshots and append-only provenance remain
in place. The active old raw relay directory is now 380 MB.

New session: `nyl12_unrestrained_es_to_ti1_20261025`. Active calculation:
`qmmm/pa6_asp272/00_es_unrestrained/`. It continues the E-S-like endpoint that
reverted from biased phase I but had not converged after 2000 cycles. Exact
controls: DFTB3/3ob-3-1, unchanged QM98/link8/charge-2 and full solvated system,
fixed cell, all 103251 atoms active, `ntr=0`, `nmropt=0`, no DISANG/reference or
frozen atoms, native DL-FIND maximum-component tolerance 1e-4 Hartree/Bohr, up
to 5000 continuation cycles. `verbosity=1` suppresses redundant multi-gigabyte
DL-FIND XYZ files; Amber NetCDF snapshots remain available.

A verified stripped snapshot of the live optimization is in
`qmmm/pa6_asp272/00_es_unrestrained/vmd/`; load `view_optimization.vmd` in VMD.
The current export contains 59 states (starting structure, all complete saved
optimizer frames, and the latest restart), 19976 atoms, all protein and PA6,
and three complete nearby TIP3P waters (topology residues 3328, 10325, 20216).
Bulk solvent and ions are excluded. The 98 physical QM atoms use Licorice,
non-QM substrate carbons are orange, and the protein NewCartoon representation
is loaded but hidden. `frame_manifest.json` records exact atom remapping,
source-snapshot hashes, and independent DCD coordinate/box validation. The
minimization remains active, so this is a progress snapshot rather than its
final trajectory; regenerate with `scripts/export_active_es_optimization_vmd.py`.
DL-FIND energy progress is extracted reproducibly from the paired
`Energy calculation finished` and `Testing convergence in cycle` records by
`scripts/plot_dlfind_optimization_energy.py`. Current CSV, JSON metadata, and
the absolute/relative energy plot are in
`qmmm/pa6_asp272/00_es_unrestrained/plots/`; rerun the script to update this
snapshot while the minimization continues. These are whole-system optimization
energies and must not be interpreted as a reaction barrier.

Overnight conditional session `nyl12_es_conditional_to_path_20261025` watches
the current 5000-cycle run without interrupting it. If that run converges, the
original controller proceeds to coordinate driving and the watcher exits
without launching anything. If it exhausts without native DL-FIND convergence,
the watcher creates `00_es_unrestrained/continuation_01/`, starts from the exact
latest restart, and permits up to 10000 additional cycles with the same DFTB3
Hamiltonian, QM98/link8/charge−2 definition, cell, tolerances, and no restraints.
Only a natively converged, finite endpoint that passes the declared noncovalent
E-S geometry gate can start coordinate driving. After an accepted continuation,
the watcher runs the existing live path publisher and resumes the 21+21-window
coordinate-driving workflow from that exact hashed endpoint. If either E-S run
fails acceptance, the chain stops explicitly. Live watcher state:
`00_es_unrestrained/conditional_continuation_status.json`; logged tmux metadata:
`00_es_unrestrained/conditional_continuation_job/`.

Only after native convergence plus the noncovalent E-S geometry gate, the same
tmux workflow proceeds to `01_coordinate_drive_to_ti1/`: 21 phase-I windows
approach OG1-C=1.55 A and transfer HG1 to Asp272 OD1, followed by 21 phase-II
windows with the attack target bounded at 1.55 A while transferring HG1 from
Asp272 toward terminal Thr N. No target is extrapolated to a shorter attack or
Asp-H bond. Tracking errors are diagnostics, not arbitrary stop limits. The
environment has no positional restraints; only the two declared reaction
coordinate biases operate in each window. Asp272 OD1 is the actual acceptor
(current HG1 distance 1.83 A); OD2 is tracked but is 3.95 A away.

After phase II, every driving bias is released in
`02_ti1_unrestrained_validation/`. A converged tetrahedral geometry that
persists there is a TI1 minimum candidate. Driven endpoints are construction
and TS-search seeds, not a barrier. The live publisher updates the endpoint VMD
movie and all declared reaction-distance plots after every completed window.
Prepared/executed input hashes are in
`jobs/unrestrained_es_to_ti1_20261025/status.json`.

## CURRENT — two original-interval phase-I windows running

IMPORTANT clarification: the ES reference and current driven path use k10
positional restraints on ALL103153 MM atoms (protein, MM substrate, solvent,
ions), in addition to any driving biases. Only the98 QM atoms lack positional
restraints. Reference is analytic_relaxation_attempt2/starting.ncrst, not the
preceding window. The user requires a truly unrestrained ES reference; original
setup does not meet that requirement and must be described as restrained.
The phase1_unrestrained branch alone removes both positional and driving
restraints. No active path Hamiltonian has been silently changed mid-run.
Native original EAMBER is -397695.4038 kcal/mol without positional bias;
latest reversion evaluation at01:20:59UTC was -443693.5466 (not converged).
Difference -45998.1428 is a whole-system potential-energy change, not a catalytic
barrier or proof of a better active site. Script/evidence:
provenance/es_energy_comparison_20261025/compare_energies.py and comparison.json.

Latest user instruction: keep the phase-I interval unchanged and append only
a few steps. Exactly stage1_17 and stage1_18 were launched, not phase II or
ten finer windows. Session `nyl12_asp272_phase1_extend_20261025`; existing run
`output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/` is unchanged.
Targets (attack Å, q1 Å):17=(1.4505601244,0.9124292155),
18=(1.3511202488,1.0248584310). Start is original stage1_16, not the reverted
unrestrained checkpoint. The window17 and16 minimize.mdin hashes match exactly;
only the harmonic targets advance by the original increments. No rotation,
Hamiltonian, restraint strength, MM reference or convergence tolerance change.

Controller `provenance/asp272_relay_20261025/extend_phase1.py`; executed frozen
copy `run/source/extend_phase1_20261025.py`; plan `phase1_extension_plan.json`.
Original status/source snapshots and an unsuccessful metadata-only preparation
are preserved. The superseded phase-II draft was never prepared/executed and
is archived as text to prevent accidental launch. Root corrected unexecuted
draft/preparation issues before actual execution;36 relevant checks passed.
Publisher `nyl12_phase1_pending_publish_20261025` updates the same
`qmmm/pa6_asp272/plots/energies.png`, `plots/distances.png`, and
`vmd/view_path.vmd` after each completed endpoint. New endpoints appear in
`01_og1_to_asp272/window_17` and `window_18`. Stop after18 for review.
These directories now expose running/queued status.json and available inputs
before optimization finishes. Only finalized endpoints enter plots/movies.
This replaces the earlier completed-endpoint-only export; no new laptop pull
script is needed. Old publisher logs and source snapshots are retained.
The separate unrestrained diagnostic and its VMD publisher remain independent.

## Concurrent diagnostic — fully unrestrained phase1 optimization

User requests observing unbiased rearrangement before considering amine rotation.
tmux: `nyl12_phase1_unrestrained_20261025`; scientific child directory:
`output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/phase1_unrestrained/`.
Original stage1_16 unchanged; copied start/frozen topology, no ntr/nmropt/DISANG,
all103251 atoms active, no SHAKE, same DFTB3/QM98/link8/charge−2/MM model and
fixed cell. Native LBFGS tolerances unchanged, maxcycle2000, ntwx5 snapshots.
No amine rotation or phase2 driving launched. Read process/acceptance records
for live outcome; launch log confirms zero restraint energy.
An intermediate restart (00:42:50 UTC Oct26) returned toward ES: OG1–C3.150 Å,
OG1–HG1 0.992 Å, Asp272–HG1 1.824 Å. Not yet a converged final endpoint;
script/source hash in `phase1_unrestrained/checkpoint_geometry.json`.
Viewer requests: dynamic reactive bonds, small HG1 sphere, O49–Phe113 N and
O49–Asn186 ND2 labels. Routine viewer/export implementation delegated to
GPT-6 Luna; native inspection needs laptop VMD for visual verification.
Stripped optimization viewer/trajectory destination:
`qmmm/pa6_asp272/phase1_unrestrained/vmd/` (publisher status in job records).

## HISTORY — phase1 complete; TI1 construction and benchmark plan

All17 stage1 windows converged; live status says stopped/phase1_complete.
Window16 has OG1–C1.795 Å, Asp272–HG1 1.022 Å, OG1–HG1 1.793 Å,
N–HG1 2.732 Å. It passes the biased stage1 gate, not an unbiased TI1 minimum;
terminal N is still NH2. Relative energy +66.11 kcal/mol including driving
bias (+60.02 after bias subtraction). Neither is a validated barrier.
Phase2 remains held for donor/amine orientation review; no new job was launched
by the literature/documentation update. Next: continue transfer to terminal N,
test bias removal/TI1 retention, reverse path, and refine saddle seeds with
DL-FIND and mode/endpoint checks. Numerical convergence and biological accuracy
are separate; high energies alone do not justify abandoning path construction.

Added manuscript introduction paragraph and R35–R40 references. Small matched
DFTB3/omegaB97X-D proton-transfer benchmark is planned, not executed; see
docs/asp272_relay_pilot_20261025.md. GFN2-xTB has a native sander interface,
but local retained build flags lack XTB support. Current g-xTB has analytic
gradients but no point-charge embedding, so cannot replace our embedded QM/MM
backend unchanged. No new dependencies or MLP training authorized/introduced.
Design references needed: ES, intended TI1 and connecting TS(s); apo control,
then acyl-enzyme/product compatibility before final candidate selection.
Latest plots/movie remain qmmm/pa6_asp272/ (17 endpoint frames).

## HISTORY — phase1 resumed with tracking warnings only

User challenged repeated arbitrary tracking cutoffs. Same run directory,
new frozen controller run_relay_warning_only.py, tmux
`nyl12_phase1_warning_only_20261025`: resume from convergedwindow13 and finish
windows14–16. Tracking warns above0.15 Å and has no hard stop. Convergence,
finite/box and chemical endpoint checks remain unchanged. Phase2 is held for
amine/Asp-OH orientation review. Status/policy history is preserved; active
policy is warning_only_tracking_policy.json. Read status.json for live progress.

Window13 energy rise:+49.67 kcal/mol, comprising+4.46 driving bias,+2.03 MM
position restraints,+46.53 printed DFTB term and−3.35 other terms. Not a barrier.
Both amine H-bonds to Asp residues were confirmed; Asp-OH···N angle95.6° is
poor. Proton transfer is ahead of covalent attack(Asp-H1.094 Å, OG1-C2.053 Å).
Current hypotheses include unfavorable coupling/orientation and proton affinity;
do not attribute the rise solely to one hydrogen bond or infer a Nyl12 pKa.
See latest section of docs/asp272_relay_pilot_20261025.md. Analyses are exposed
in qmmm/pa6_asp272/data/amine_orientation.json and relay_energy_component_audit.json.
Earlier stopped/running statements below are history.

## LATEST — energy plot published; relay stopped at window13

Energy plot: `qmmm/pa6_asp272/plots/energies.png`; numerical data:
`data/energies.csv` and `.json`. Includes all14 displayed endpoints, with
window13 explicitly marked rejected. The optimizer converged there, but
tracking exceeded0.20 Å. No further scientific continuation is queued.
The package publisher now regenerates energy plots whenever a new finalized
inspection snapshot appears. Original calculation files remain unchanged.
Earlier running statements below are historical; read live status.json.

## CURRENT — same-directory Asp272 continuation, 2026-10-25

Preferred inspection entry: `qmmm/pa6_asp272/README.md`. Plots are in
`plots/distances.png`, VMD in `vmd/view_path.vmd`, and the executed source in
`scripts/run_path.py`. Laptop: `nyl12-inspection-20261021/qmmm/pa6_asp272/`.
Numbered elementary stages contain sequential `window_00`, `window_01`, ...
The active calculation stays in its original directory.
The preferred VMD viewer now labels the three reactive distances per stage;
values update with animation and labels switch at the stage boundary.
Publisher session: `nyl12_qmmm_labels_20261025`.

Window10 converged but missed the 0.150 Å tracking gate by 0.003764 Å.
User authorized warning at 0.150 Å and stop at 0.200 Å, explicitly keeping the
same run directory. Continuation began 18:37 EDT at stage1_11 in tmux
`nyl12_asp272_continue_20261025`. Run remains
`output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/`.
Original window10 acceptance and inputs are preserved; reassessment and source
hashes are in continuation_tracking_policy.json; prior status is preserved.
No force constant, Hamiltonian, target schedule, optimizer tolerance or stage-end
chemical gate changed. This is not a TS/barrier claim.

Each new finalized window now synchronously regenerates a verified endpoint
VMD trajectory and 16 reaction-distance traces, with CSV/JSON and jump flags.
Stable files in `output_vmd_inspection_20261024/`: reaction_path.vmd,
reaction_path.dcd and reaction_distances.png/.csv/.json. All are included in
the no-argument laptop pull. Reload the VMD viewer to see newly added frames.
A 0.25 Å adjacent-window jump is a review warning, not automatic rejection or
an instruction to add constraints. No such jump occurred through window10.
Detailed notes, expert references and execution protocol:
`docs/asp272_relay_pilot_20261025.md`. Read live status for actual progress.

## CURRENT — Asp272-mediated relay running, 2026-10-25

User authorized an Asp272 two-stage HG1 shuttle. tmux:
`nyl12_asp272_relay2_20261025`. Run:
`output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/`.
First native L-BFGS energy/gradient evaluations and trajectory writes verified;
no window convergence is claimed at launch. Stage 1 couples OG1–C attack
to HG1 transfer onto Asp272 OD1. Stage 2, conditional on the intermediate
gate, transfers that same proton from Asp272 to terminal N while retaining
the attack geometry. Native maximum-component convergence plus target tracking
replace the earlier permissive global RMS acceptance. Read status.json for
actual progress. No new water, TI1 parameters, reverse scan or MD are queued.
Controls/provenance/acceptance: docs/asp272_relay_pilot_20261025.md and
configs/asp272_relay_pilot_20261025.json. Failed namelist-order attempt 1 is
preserved; use attempt 2. Earlier stopped/running statements below are history.

## PREVIOUS — direct-N forward pilot stopped before reverse, 2026-10-25

Latest: stricter ES stopped at cycle 223 with REPEATED LINMIN FAILURE;
RMS 0.016484 remains above 0.0005. No scientific job is currently active
from this work. User-observed Asp272 contact confirmed: OD1 accepts both
N–H1 and OG1–HG1 hydrogen bonds, whereas initial OG1–HG1···N alignment is
poor (96.3°). Do not assume a blocked amine, a necessary relay water or a
verified mechanism. Compare direct-N and Asp-assisted hypotheses after
repairing optimizer convergence. Notes: docs/asp_thr_proton_network_20261025.md.
VMD non-QM substrate carbons are orange; protein NewCartoon remains available
but is hidden by default. Shared viewer and package checksum were updated.
Running language below is historical; check final process/acceptance records.

User requested an actual barrier and then a movie of the failed scan.
The concatenated forward movie is now in the curated VMD download:
`output_vmd_inspection_20261024/forward_driven_path.vmd` and matching DCD.
It contains 18 actual minimized window endpoints (19976 atoms, three waters,
98 real QM atoms in Licorice), without interpolation or alignment. All windows
finished before their first ntwx=10 snapshot; no within-window intermediate
frames were saved. The frame manifest records target and achieved coordinates.

A separate tighter ES reference is running in tmux
`nyl12_es_barrier_reference_20261025`, directory
`output_nyl12_pa6_dftb_pilot_20261024/es_tight_convergence_attempt1/`.
It starts from the original ES endpoint, retains the same QM/MM model and
original MM restraint reference, and changes drms to 0.0005 with a 3000-cycle
budget. No new driven path or TS search has started. Check its process and
scientific acceptance records, then review chemistry before using its energy.
This stricter minimization is not itself a barrier calculation.

All 18 forward windows completed, but the endpoint failed the TI1 geometry
gate; reverse did not run. No path job remains active. Final OG1–C 2.60843 Å
(target 1.55), proton difference −0.73086 Å (target +0.50). The seven-cycle
final minimization met global RMS 0.04279 but retained a 6.8292 kcal/mol/Å
maximum force on QM atom O49. Global RMS acceptance was too permissive for
the local chemistry. Do not claim a complete reaction path, barrier or unstable
TI1. Review/plots: docs/pa6_overnight_path_review_20261025.md and
output_nyl12_pa6_dftb_pilot_20261024/overnight_review_20261025/.
Next: verify stricter local convergence and target tracking in a few windows
before another full forward/reverse attempt. Older running status below is
historical. Original outputs and all forward movies remain preserved.

## OVERNIGHT LAUNCH — forward/reverse TI1 construction, 2026-10-24

User authorized a sequential overnight coordinate-driving pilot and then
disconnecting. tmux: `nyl12_ti1_forward_reverse3_20261024`.
Run directory:
`output_nyl12_pa6_dftb_pilot_20261024/coordinate_driving_forward_reverse_attempt3/`.
First three forward windows completed and exported; forward_03 was running at
the checkpoint. Read status.json and the job log for current progress.
There are 18 forward targets, followed by 17 reverse targets only if the
forward endpoint passes a candidate TI1 geometry gate. Each window inherits
the preceding converged restart and retains the original MM positional
reference, full system, QM region and Hamiltonian. Native NMR restraints drive
OG1–C and OG1–HG1 minus N–HG1 independently. Stop on engine failure,
nonconvergence, invalid coordinates/box or display-export failure. No MD or
new parameter fit runs tonight. This is a biased construction pilot, not a
validated reaction pathway or barrier calculation. Global RMS convergence
alone does not establish QM stationarity or a TS. Review actual coordinates.

Every completed window retains full optimizer snapshots (ntwx=10), its exact
endpoint, restraint targets, process/geometry records and
optimization_stripped.dcd. A window converging before the first snapshot
contains only its exact endpoint in that display file, explicitly recorded.
The fixed 19976-atom display selection matches the VMD bundle topology.
Attempts 1/2 are retained: DISANG's 80-character filename field required
relative paths; a zero-snapshot immediately converged window required an
explicit endpoint-only display export. Neither changed the chemistry.

VMD bundle: `output_vmd_inspection_20261024/` (17.1 MB, 12 checked payloads).
Default `nyl12_es_optimization.vmd`: 48 ES frames, protein/full substrate,
three fixed waters and QM Licorice. `nyl12_es_initial_final.vmd`: full-solvent
two-frame comparison. Updated `scripts/pull_nyl12_inspection.sh` downloads
both curated packages. Next session: inspect overnight progress and endpoint
chemistry; package forward/reverse path movies for the laptop before expanding
scope. Details: docs/pa6_ti1_coordinate_driving_plan_20261024.md.

## ES REFERENCE CHECKPOINT — 2026-10-24

User inspected the 98-real-atom + 8-link-H QM XYZ and approved the region,
including Asp272/Asp274 and the following Thr234. Keep that region unchanged.
The full 103251-atom solvated tetramer and full PA6 L4 remain in the calculation.
The 25-cycle analytic-gradient pilot completed but did not converge: RMS
2.4646 versus target 0.05 kcal/mol/Å. Thr retained OH/NH2 and the substrate
retained an amide-like geometry; this is ES, not TI1.

A continuation with a 500-cycle budget completed in tmux
`nyl12_es_continuation_20261024` after 107 cycles, exit 0. Its final RMS
gradient is 0.042903 kcal/mol/Å, below the 0.05 target. The fixed box and
finite full-system coordinates were preserved; 21 optimizer snapshots saved.
Directory: output_nyl12_pa6_dftb_pilot_20261024/analytic_relaxation_continuation_attempt1/.
It starts from attempt2/optimized.ncrst and retains attempt2/starting.ncrst as
the original restraint reference. Only cycle budget and output cadence change;
QM atoms remain free, all MM atoms retain the same positional restraints.
Endpoint review retains Thr OH/NH2: OG1–HG1 0.98493 Å, N–H distances
1.02118/1.01706 Å, substrate C–O/C–N 1.25974/1.35877 Å. No TI bond
formed. Maximum QM/MM heavy-atom shifts from the original restraint reference
are 0.28286/0.32613 Å. Accept as a restrained ES reference for subsequent
model construction, not as evidence of a reaction mechanism. Final inspection
export is complete in the continuation directory: inspection/ contains the
QM XYZ, full enzyme/substrate endpoint and 21 saved optimizer states; the
final cycle-107 endpoint is separate from the last saved trajectory frame.
Final d1/d2/d3 = 3.14104/3.02079/2.84060 Å, a_BD = 102.19064°,
herschlag_phi = −89.60256°. Next: construct an enzyme-supported TI1 seed
from this restrained ES reference; native TI1 parameters remain unaccepted.
This reported-gradient convergence does not establish a Hessian minimum or
validate reaction energetics. No MD, reaction path or new parameter fit ran.

The attempt2/inspection/ directory now contains qm_subsystem_endpoint.xyz,
qm_subsystem_initial.xyz, an explicit link-atom map, endpoint_review.json,
full_enzyme_substrate_endpoint.pdb.gz and portable viewers. Pocket PDBs are
display extracts only. Native TI1 parameters remain unvalidated and paused.
Own-code verification: 656 passed/13 skipped in the standard run, then all
13 workstation-fixture cases passed separately; new-code coverage 91.76%.

Keep one active scientific milestone: finish ES minimization and review
geometry, then establish an enzyme-supported TI1 model before a small
symmetry-aware fixed-backbone LigandMPNN pilot. Full pathway/free-energy work,
local LLM agents, NylB and variable-length RFD3 work remain deferred.
Use bounded cheaper delegation for routine tasks and short progress reports.
Active append-only journal: output_agent_actions/20261021.jsonl.
Visualization policy: retain actual optimizer coordinate snapshots, original
starting coordinates and exact final endpoints for each future reactant,
intermediate, TS and product optimization. Set Sander trajectory output
explicitly before launch; preserve frame/source mapping and periodic boxes.
Optimizer frames are not a physical-time reaction trajectory. Current ES
already has 25 + 21 saved snapshots, so no rerun is needed for visualization.
Older sections below are historical and do not override this checkpoint.

## RESUMED: analytic-gradient bounded ES pilot — 2026-10-24

User explicitly resumed next steps, then questioned the numerical-force gate
and asked whether analytic derivatives are available. Sander already supplies
analytic QM/MM gradients; its force implementation was not modified here.
Stop further finite-difference calculations. The earlier failed numerical
check remains failed and its cutoff explanation remains unconfirmed. It is
not a chemically meaningful reaction-energy error or evidence against TI1.
Root narrowed the next calculation to an explicitly declared, exploratory
25-cycle restrained ES relaxation with analytic gradients, retaining the
existing finite analytic baseline, source/atom/proton/charge/embedding checks.
No force tolerance is changed, and no chemical-model validation is claimed.
An explicit pilot-decision record and separately tested driver mode are being
prepared; do not silently substitute this for the old strict gate. All MM
atoms remain restrained at 10 kcal/mol/Å² to the original coordinates; QM
atoms are free and no reaction-coordinate restraint or MD is introduced.
Actual new job status must be recorded below before the next checkpoint.
The older 15:00 stopping time was superseded by this explicit user resume.

## PAUSED after enzyme force validation — 2026-10-24

User requested work until 15:00 America/New_York, then a break. A connection
disruption was followed by a clock check at 19:07:54 UTC (15:07:54 EDT); root
immediately paused scientific work. No new scientific jobs were launched after
that check. All delegates stopped. Only the two unrelated original tmux
sessions remain. Resume only when the user authorizes more work.

The approved 12 archived regression inputs/logs were committed and pushed as
3b1cd4f. New pilot code/config/tests/docs remain uncommitted pending complete
verification; unrelated untracked cif/ remains untouched.

Independent chemical audits passed for the primary 98 original QM atoms plus
8 link hydrogens (106 nuclei, 368 electrons) and a 104+8 boundary variant
(112 nuclei, 384 electrons). Both are charge −2 singlets. Config:
configs/pa6_es_dftb_pilot_20261024.json. The full original PA6 L4/tetramer/
solvent MM environment is retained; no forced water or TI1 classical parameters.

Periodic SANDER API NMA/water smoke passed after explicit setup corrections.
Its evidence and failed attempts remain in output_nyl12_es_dftb_api_20261024/.
The actual enzyme baseline produced finite energy/forces for all 103251 atoms,
correct 98+8 inventory and unchanged box; energy −332062.9867349793 kcal/mol.
This absolute energy is not a reaction barrier or chemistry validation.
Baseline: output_nyl12_pa6_dftb_pilot_20261024/api_primary_baseline_attempt1/.
Explicit controls: DFTB3/3ob-3-1, cut/qmcut 10 Å, electronic embedding,
adjust_q=2, QM Ewald/PME, SCF 1e-10, C–H links 1.09 Å, fixed box, no SHAKE.
The engine zeros MM link-partner charges and distributes the charge correction
uniformly over the remaining MM atoms, preserving the original −0.005001934 e
system charge residual. Original topology and protonation were not modified.

**Force gate FAILED overall.** All 26 isolated single-point evaluations
completed; 11/12 comparisons passed. C54 x at h=0.001 Å: analytic force
29.3573104, finite difference 29.0911276 kcal/mol/Å, error 0.2661828 above
tolerance 0.0493573. At h=0.0005 Å the error is 0.0005500 and passes.
Baseline energy repeats exactly; source hashes are unchanged. Do not loosen
the threshold or call the model validated. Next diagnose the discrepant endpoint
(cutoff/neighbor-list crossing versus SCC or another cause) with targeted
repeats and displacement convergence. These causes are hypotheses, not findings.
Evidence: output_nyl12_pa6_dftb_pilot_20261024/api_primary_force_check_attempt1/
force_check.json and scientific_acceptance.json; frozen source is retained.
No minimization, dynamics, TS search, PMF or hybrid-DFT chemical check ran.
The boundary variant has passed only the input audit, not an engine calculation.

Prepared restrained-relaxation driver remains review-only and fails closed on
the rejected force evidence. Its tests are unexecuted and additional gate
tightening was discussed but not implemented. Do not run it before review and
passing force validation. New QM-region tests initially passed 13 cases/87%
coverage; final whitelist changes have actual-input validation but no final
test rerun. Engine/controller mocked tests completed 44 passed, 97.06% combined
coverage, at 19:08:15 UTC before the stop message reached that delegate. No
full-project regression was run for these new files. Actual tests/job sources:
output_nyl12_pa6_es_audit_20261024/jobs/engine_mock_tests_attempt1/.

Native TI1 parameterization remains paused at scientific validation. See
docs/ti1_parameterization_status_20261024.md. No native RESP fit, validated
TI1 classical force field or TI1 protein refinement is accepted.

The local-agent mechanism project is deferred for later at the user's request.
Chinese-origin LLMs and their derivatives are prohibited. No model/runtime was
installed or downloaded. Notes: docs/local_agent_mechanism_project_20261024.md.
The two new status/agentic notes are saved locally but the laptop export has not
yet been refreshed for this pilot. Active journal: output_agent_actions/20261021.jsonl.

## WT PA6 ES recovery and hydration completed — 2026-10-24

User accepted the staged project and authorized progress today. Audited three
normally completed 100 ns WT runs with distinct seeds, ff19SB/GAFF1/TIP3P and
Thr OH/NH2. Nine full-solvent seeds extracted in tmux with source hashes;
native CIF residue identities checked. Native Asp272 OD1–Thr233 N = 2.83048 Å;
no catalytic role or literature homology assigned. See
docs/pa6_es_reuse_audit_20261024.md and
output_nyl12_pa6_es_audit_20261024/reuse_audit_attempt1/reuse_audit.json.
Installed DFTB3 water/NMA regression energies match references: engine QA only.
All nine seeds passed independent geometry/source checks. Full saved-frame
hydration analysis found no water O within 4 Å of occupied-site Thr N/OG1 in
any of 6000 frames (3 × 100 ns, 50 ps cadence); no simultaneous 3.5 Å contacts.
Water-centered positive control matches independent periodic counts. This is
conditional ES-model evidence, not a universal no-relay claim or deacylation
model. See docs/pa6_es_hydration_20261024.md and the portable PyMOL export in
output_nyl12_pa6_es_audit_20261024/seed_inspection_attempt4/.
559 tests passed, 13 skipped; new-code coverage 95.12%. Exact old-PDB phi and
raw-NetCDF values both retained after identifying precision tie rounding.
User clarified that n12c.cif comes from a coauthored manuscript in revision;
the formal citation will be provided after publication. Record it as the
supplied structural reference; do not invent a public accession/citation.
No new production MD/PMF, TI1 fit or TS accepted. Next define the locally
truncated QM core, boundaries and charge inventory, then bounded DFTB3/MM and
selected hybrid-DFT PES checks. See docs/pa6_dftb_validation_pilot_20261024.md.

## Proposed separate WT PA6 mechanistic project — 2026-10-24

User supports considering a full chemical-cycle DFTB QM/MM free-energy study
as a standalone paper. Focus PA6, optionally poly(MAXA), without expanding PA66.
Use hybrid DFT only for selected PES optimizations/checks, not free-energy
sampling. QM/MM validation and convergence are substantive research tasks.
Prior antechamber/AM1-BCC WT noncovalent PA6 L4 MD (3 × 100 ns) already exists;
corrected records show PA66 mainly loses local geometry/register, without
necessarily dissociating. Pin valid topologies/trajectories before reuse;
do not repeat simulations or change charges without an explicit rationale.
NylB is a comparator, not an authorized pivot; its published 2014 reactive
QM/MM study is DOI 10.1021/jz500323y. Scope and validation gates, including the
distinction between a reactive QM core and state-specific MM parameters, are
in docs/nyl12_mechanistic_project_scope_20261024.md. No new simulation started.
A full-cycle mechanism would support a paper; it is not required for every
exploratory sequence. A PMF maximum is not automatically a PES saddle point.
Updated scope and references are included in the existing laptop pull package.

## Literature/manuscript checkpoint — 2026-10-24

User requested preserving relayrefs and recent2025/2026 Baker/othergroup
lessons for Nyl12 redesign. R23--R33 and new BibTeX records now retained in
manuscript/. Introduction/Methods/sourceledger updated; proposals remain
labelled. New docs/recent_enzyme_design_lessons_20261024.md compares final
2025serine/2026metal, Oct14proteases, Oct15RFD3, nativefold/PET alternatives.
Primary-access/preprint limits explicit. Pseudo-ring excluded only in tested
fixedmodel, not a universal/transientrelay. No new QM/design/MD campaign.
Recommended boundedchemistry models and multistatescreening with fulloligomer
conformers/registers/termini; trueTS1 validation remains required for saddle/
barrier claims, not a prerequisite to every exploratory structural screen.
No new dependencies/software adopted. Laptop export includes updated drafts,
references, machine-readable BibTeX and both research documents.

## Mechanistic interpretation extended — 2026-10-24

User reconsiders frontwater and requests literature beyondNylC. Found positive
primaryproteasome QM/MM example Wei2013 DOI10.1021/jp405337v and competing
water-free protontransfer model in penicillinacylase Chilov2007. NylC2023
FEBS DOI10.1111/febs.16755 proposes Asp-Asp-Thr network for maturehydrolysis;
its precursorwater/autocleavage figures must not be conflated with catalysis.
No definitiveNyl12 mechanism assigned. See newmechanism evidence doc.
Oppositehemisphere diagnostic also0/225000 full andheavy-only. No newseed
orGaussian. Assess currentAF3-derivedconformation andmapliteraturecatalytic
residues before forcingwater orfittingassumedNH3+/oxyanion state.

## Resumed feasibility search completed — 2026-10-24

User resumed frontwater placement and parameterization/testing if feasible.
No acceptable front oxygen in fixed complete protein/substrate pocket:
0/225000 samples at2.6--3.2A N/OG1 contacts; heavy-only0/225000 as well;
expanded3.3A diagnostic0/302760. Best fullmargin−0.684A (LIGH23), heavy-only
−0.217A (substrateC70 plusAsp272 contacts). No new Gaussian or parameterfit.
RotatingNH3 or waterH cannot rescue an oxygen-heavyatom clash. This is a
finite fixed-conformer search, not an impossibility proof for the enzyme.
See docs/front_water_feasibility_20261024.md and
output_nyl12_front_water_20261024/. Full-pocket viewer uses two rejected
oxygen probes, not seeds. Next resolve protein/substrate geometry or choose
another validated conformation without bypassing TI1-FF validation requirement.
Eight SMDacceptance corruption tests nowpass, including lostbridge and reverse
protontransfer, with93.94% assessor linecoverage. Source hashes, journal and
all failed attempts retained. Earlier PAUSED/ACTIVE notes below are history.

## PAUSED at user request — evening 2026-10-23

No Gaussian process or tmux calculation from this effort remains running.
Other three existing tmux sessions are unrelated and untouched. Stop work.
Latest user XYZ is preserved exactly and checked; it clears the substrate but
clashes with Asp topology272 backbone in the original full pocket. Resolve
full-pocket geometry before treating that water placement as enzyme-compatible.
No seed/reference/FF has been accepted. No automatic restart or continuation.
SMD remains a diagnostic proposal with one user-rejected-placement pilot stopped;
that interrupted run does not establish whether continuum stabilizes TI1.
See detailed notes below and retained output roots. Active journal remains
output_agent_actions/20261021.jsonl. Unrelated untracked cif/ is untouched.

## User-supplied manual water seed checked — 2026-10-23

Saved exact supplied XYZ: output_nyl12_user_water_seed_20261023/native_seed_jmp.xyz.
Script/check_job/seed_checks.json preserve measurements and input hashes.
63atom inventory, graph and S/R/R stereo pass. WaterO-N2.67154A,
waterO-OG12.90269A; waterH-OG1 angle166.53deg, best NH-water117.73deg.
Nearest substrateheavy3.629A. However FULL original AF3-derived alpha/beta
pocket has severe clashes: ASP topology272 CA1.44314A, HA0.81751A from
waterO. This is model-context evidence, not experimental steric assignment.
Do not treat the isolated-cluster placement as enzyme-compatible. Native
pocket geometry/possible rearrangement must be resolved separately. User
coordinates remain unchanged; no new Gaussian calculation was launched.

## STOP requested: SMD seed water rejected by user — 2026-10-23

User inspected seed and requires water IN FRONT of amino terminus, away
from the substrate-side placement shown. Identity-verified Gaussian worker
4055951 sent SIGTERM; interruption record in SMD pilot. No replacement run.
Check job/status.json for terminal record. Preserve all files. This rejects
the water geometry, not SMD as a method. Earlier ACTIVE note is superseded.

The deterministic search was insufficient: only water moved, so original
NH3 proton directions biased donor selection toward H3; geometry screening
used covalent-radius contacts in truncated63atom cluster, not full-pocket
van der Waals clearance. Distance/angle satisfaction alone is not acceptance.
Next place/check against FULL substrate/protein context, define exposed face
in an atom-based frame, and allow NH3 H reorientation with water placement.
Do not restart Gaussian with another unreviewed side/clearance assumption.
User's five-minute departure constraint remains relevant; no extended work
or autonomous continuation requested. SMD/RESP protocol distinction preserved.

## ACTIVE SMD/water-bridge diagnostic — 2026-10-23

User requested a decision and tmux launch before leaving in5min.
New root: output_nyl12_smd_water_bridge_20261023/.
ACTIVE tmux: nyl12_bridge63_smd_20261023; pilot svp_pilot_01/.
63atoms, neutral singlet, wB97XD/def2SVP SCRF=(SMD,Solvent=Water),
16cores32GB,20cycles maximum, same frozen21heavy plus chi1.
No reactive-bond/proton/water restraints during optimization; no Hessian/RESP.

Prior water was not a validated bridge. New deterministic search moves ONLY
water O/H/H, retains every nonwater coordinate, and screens graph/clashes.
1800 orientations evaluated; selected N-Owater=OG1-Owater=2.70A,
N-H3-Owater=151.38deg, Owater-H-OG1=180deg, OH0.96A/HOH104.5deg.
Initial63atom geometry passes. Final acceptance adds both water heavy-atom
contacts2.5--3.3A and best NH-water / waterH-OG1 angles>=140deg, alongside
existing NH3, oxyanion, graph/stereo and convergence checks. A retained
state would be a constrained aqueous-continuum cluster diagnostic, not a
native enzyme minimum or accepted FF reference. Water continuum is not a
measured enzyme dielectric. Do not autoextend or fit parameters.

Ruff passes and real seed/source geometry checks pass. New bridge-gate unit
coverage has not been measured under the user's5min deadline. Earlier12tests
apply to the preceding gas-phase pipeline, not these added water criteria.
Source snapshots, input hashes and action journal retained. Job manifest pins
launched assessor; separate assessment_source_hash_at_checkpoint.json discloses
preformat digest in earlier staging source_manifest. No method silently changed.

Literature: Meadows etal2026 JCTC22,4012--4025 DOI10.1021/acs.jctc.6c00162
explicitly uses implicit water to stabilize glycine zwitterion for GAFF/RESP
parameter variants: https://pmc.ncbi.nlm.nih.gov/articles/PMC13130867/ .
Enzyme SMD dielectric4 precedent: https://doi.org/10.3389/fchem.2018.00249
(singlepoint enzyme energies; not proof SMD Hessians transfer to Amber).
Recommendation: solvated geometry is defensible; ESP charge protocol and
bonded fitting remain separate decisions. Do not adopt fully solvated ESP
charges or raw continuum Hessian as automatic ff19SB/GAFF2 parameters.
No replacement RESP or force-field fit launched. Diffuse/larger-basis and
actual protein-environment checks remain pending if TI1 survives.

## Pilot ended: TI1 state loss; trajectory available — 2026-10-23

The63atom donor/water pilot ended17:10:33EDT at20cycle limit, exit1.
No convergence and no TI1 acceptance: OG1-C54=1.90553A, N-H3=1.40361A,
OG1-H3=1.03419A, consistent with reverse proton transfer. WaterO-N=3.30394A
and waterO-OG1=3.53239A; an acceptable pseudo-ring is not established.
No continuation is running. Earlier ACTIVE entries below are historical.

Exported21complete Input orientation tables,63atoms and37explicitH each,
with both water H. Workstation snapshot:
output_nyl12_oxyanion_cluster_20261023/trajectory_snapshot_01/.
The usual rsync pull now includes
qm/ti1_oxyanion_cluster/optimization_progress/view_optimization.pml,
optimization_trajectory.pdb and .xyz, CSV metrics and acceptance report.
Run @view_optimization.pml then mplay from that directory. This is a
snapshot, not a live feed; fixed PDB bonds do not prove connectivity survived.
Source log snapshot/hash, exporter and tmux records preserved. Ruff and
PDB frame/atom/H roundtrip checks pass. No new QM or force-field run launched.

## ACTIVE: oxyanion-donor/water cluster pilot — 2026-10-23

User approved N-methylacetamide for Phe113 backbone NH, acetamide for Asn186
side-chain amide, and relay water, then selected omegaB97X-D/def2-SVP.
New cluster: output_nyl12_oxyanion_cluster_20261023/seed_01/; 63 atoms,
C15H37N5O6, charge 0, singlet, 37 explicit H. Components: 39-atom TI1,
3-atom water, 12-atom NMA, 9-atom acetamide. Original 42 seed coordinates
unchanged; protein frame verified without alignment; R at C54 retained.
Phe H--O49=1.7133 A; Asn HD22--O49=2.1559 A. Seed geometry screen passes.

ACTIVE tmux: nyl12_oxyanion63_svp_20261023, started 16:40:57 EDT.
Job: output_nyl12_oxyanion_cluster_20261023/svp_pilot_01/job/.
Gaussian input/policy/log in parent svp_pilot_01/. Uses 16 physical cores,
32 GB, ReadAllGIC, MaxStep=2, at most 20 optimization cycles; no Hessian,
extra dispersion or solvent. Short state-retention pilot, not a force-field
reference. All 37 H and 5 reactive/water heavy atoms move. Freeze original
12 scaffold heavy atoms, 9 added donor heavy atoms at protein positions,
and source AF3 chi1. No reaction-bond/proton restraint. Donor heavy freezes
represent held protein geometry, not a freely optimized isolated cluster.
Remaining protein electrostatics are absent; enzyme TI1 is not established.

Twelve focused tests pass in scripted tmux; Ruff passes. Source line coverage:
assessor 93.44%, donor helper 93.02%, cluster builder 94.52%, pilot preparer
91.94%. New gate checks all three ammonium N--H (0.85--1.25 A) and their
OG1--H distances (>=1.20 A), in addition to graph/stereo/oxyanion checks.
Final assessment runs even if Gaussian exits unsuccessfully. A cycle-limit
exit is not convergence; retention over a short pilot is not a proven minimum.
Inspect actual job status and log before continuing. No RESP, FF acceptance,
protein refinement, MD, AF3 or redesign campaign was started. Earlier stopped
water-free39 data remain preserved and rejected as a TI1 reference.

Seed construction/preparation/test scripts, source snapshots, hashes and all
failed preparation attempts remain in the new output root. Routine fragment
scripting delegated to Sol High and reviewed. Active journal unchanged:
output_agent_actions/20261021.jsonl. Unrelated cif/ and other tmux untouched.
The laptop package includes qm/ti1_oxyanion_cluster/view_seed.pml and seed,
plus pilot/input.gjf and pilot/pilot_policy.json. Next: assess state retention,
not automatic extension, charge fitting or force-field transfer.

## STOPPED: water-free TI1 reverses to reactant-like state — 2026-10-23

User requested stop after observing reverse proton transfer. The verified
worker3978424 receivedSIGTERM; jobended16:20:27EDT/20:20:27UTC,exit1.
All files preserved. No Gaussian job from this effort is now running and no
replacement model/method was launched. Earlier active-job notes are history.

Stopped19-frame analysis confirms ThrN–H3 lengthens1.010->1.488A while
OG1–H3 shortens1.993->0.916A; OG1–C54=2.139A, C54–O49=1.220A lastframe.
This is not a TI1 fitting reference. See
output_nyl12_waterfree_reference_20261023/state_loss_review_01/report.json.

Source-protein candidate oxyanion donors: Phe113 backboneNH (H–O49=1.713A)
and Asn186 ND2-HD22 (2.156A). Proposed next model restores relaywater plus
these actual donor fragments, retains shortsubstrate, and defines protein
boundary/embedding treatment. Exact caps still need mapping; no new run.
See docs/ti1_state_loss_and_qm_options_20261023.md for the measured evidence,
proposed~63atomlocalcluster, proton-transfer screening improvement and QM-cost
options. Basis/functional pivot is not yet adopted. KeepexistingTZVPdata;
considerlower-costbasis with samefunctional and targetedlarger/diffusebasis
checks. No RESP/forcefield/MD accepted. Livejournal unchanged.

## Trajectory inspection warning — 2026-10-23

A 14-frame snapshot is exported at
`output_nyl12_waterfree_reference_20261023/trajectory_snapshot_01/`.
At its last frame OG1-C54=2.02549 A and C54-O49=1.21644 A; these exceed
our TI1 acceptance windows, and frames13/14 fail the graph geometry screen.
This is an in-progress snapshot, not a converged chemical-state conclusion.
Do not fit parameters from this candidate without resolving state retention.
Input graph connectivity in the PDB is for display only. The movie, XYZ,
per-frame metrics and plain-command PyMOL viewer are in the laptop package.
The user requested proportionate work and fewer tokens/credits for simple
tasks: use existing tools first and avoid expanding small tasks unnecessarily.

## Water-free 39-atom candidate active — 2026-10-23

The user authorized removal of water for the local parameter reference.
The former 42-atom job FAILED at 13:49 EDT: five energy/gradient points,
23 min 52 s, l103 “Error imposing constraints”, exit 1. Preserve its entire
truncated_pilot_01 directory and failure_assessment_20261023.json. There is
no accepted geometry; this is a numerical constraint failure, not a conclusion
about TI1 stability. The old 103-atom job remains stopped by user request.

ACTIVE tmux: `nyl12_waterfree39_gic_16core_20261023`.
Directory: `output_nyl12_waterfree_reference_20261023/waterfree_pilot_01/`.
39 atoms, C10H23N3O3, 23 explicit H, net charge 0, singlet. Thr NH3+, O49-
and C54 R retained. Remove exactly water O/H1/H2 from original 42-atom seed;
all remaining coordinates unchanged. Do not restart from the failed endpoint.

Same method/resources: wB97XD/def2TZVP, 16 physical cores, 32 GB. Numerical
remedy: explicit Cartesian ReadAllGIC plus frozen chi1, MaxStep=2. Same 12
frozen heavy atoms, 27 mobile atoms, no reactive-bond restraint. AF3-seed
scaffold/torsions remain distinct from n12c reference values. GIC toy test
preserved frozen XYZ and torsion; eight focused tests pass (87% focused-run
coverage, preparer 94%, assessor 93%). Gaussian accepted the full GIC input.

Next inspect actual convergence and chemistry; a separate exact Hessian plus
gradient/constrained curvature follows only after geometry passes. No Hessian,
RESP fit, accepted force field, protein refinement, MD or AF3 job was started.
Water belongs to later enzyme/QM/MM validation with standard TIP3P in MM.
See docs/native_ti1_reference_20261023.md. Action journal remains
output_agent_actions/20261021.jsonl; unrelated cif/ preserved. Earlier active
42/103-atom notes below are historical and superseded.

## Truncated 42-atom pilot replaces 103-atom pilot — 2026-10-23

The user explicitly authorized truncation and more cores on violet. The old
worker PID 3922102 was identity-checked and sent SIGTERM; its job ended with
exit status 1. Preserve the old native_frozen_scaffold_pilot_01 directory,
including user_requested_interruption.json. It did not complete; this is
supersession rather than chemical rejection.

ACTIVE tmux: `nyl12_truncated_ti1_16core_20261023`.
Directory: `output_nyl12_truncated_reference_20261023/truncated_pilot_01/`.
The model has 42 atoms, formula C10H25N3O4, 25 explicit H, net charge 0 and
singlet multiplicity. It retains Thr NH3+, the oxyanion and explicit water.
C54 is R; Thr CA/CB are S/R. Retain two substrate carbons on each side; cut
C70-C78 and C75-C79 and add H at 1.09 A along the old bond directions.
All 40 retained atom coordinates are unchanged.

Gaussian uses 16 physical cores (CPU affinity 0-15), 32 GB RAM and 677 basis
functions. Method: wB97XD/def2TZVP optimization only, no CalcFC/Freq, extra
dispersion or solvent. Same AF3-seed chi1 and frozen-scaffold pilot policy:
12 frozen heavy atoms, 30 mobile atoms, reactive bonds free. This does not
replace the n12c torsion reference. Eight focused tests pass: 97% overall
coverage, 94% preparer and 95% assessor. See docs/native_ti1_reference_20261023.md.

Next assess actual convergence, chemistry and water geometry, then consider
a separate exact Hessian, gradient and constrained-curvature assessment.
No native RESP, force-field transfer or protein relaxation is accepted.
No MD: the user's proposal remains explicitly deferred. Active journal:
output_agent_actions/20261021.jsonl. Unrelated cif/ is preserved. The earlier
103-atom running-job note below is historical and superseded.

## Native TI1 pilot started — 2026-10-23

User inspected capped geometries and authorized the next step. Native seed
mapping establishes **R** at C54 (Thr CA S, CB R); do not select S to evade
R cap sensitivity. A source-mapped103atom/60explicitH charge0 native seed
passes initial graph/stereo/contact screening. Water comes from the earlier
physical water-only TIP3P restart, not the manual displacement or failed DFTB
endpoint. See docs/native_ti1_reference_20261023.md for exact construction.

Active tmux: `nyl12_native_ti1_pilot_20261023`; directory
`output_nyl12_native_reference_20261023/native_frozen_scaffold_pilot_01/`.
Gaussian wB97XD/def2TZVP,8CPU32GB;65mobile atoms,38frozen heavy atoms plus
seed chi1 constraint. This is a frozen-scaffold state-retention pilot using
explicitly labeled AF3-seed torsions, not a replacement for n12c-derived
parameter-reference torsions. No reactive-bond restraint, solvent or extra
dispersion. Engine started; no convergence or scientific acceptance yet.
Six checker corruption tests passed;95% assessor coverage. Tiny Gaussian
freeze test confirmed fixed coordinates unchanged. Preserve all source hashes.

The user suggested eventual explicit-solvent MD to find relay-water snapshots
for QM/MM, and explicitly said **do not do this yet**. Proposal documented;
no MD or new AF3 run. Parameter-transfer/protein-relaxation gate stays closed.
Next inspect pilot process and scientific checks; do not auto-start RESP,
production MD, or repeat the old DFTB repair. Active journal remains
output_agent_actions/20261021.jsonl. Unrelated untracked cif/ preserved.

# Standalone pipeline handoff — 2026-10-20

## Final follow-up checkpoint — 2026-10-23

All four capped-surrogate HF/MK ESP datasets completed and passed separate
scientific acceptance (8606/8843/8712/8897 grid points). No RESP fit or native
transfer occurred. The signed raw-Hessian R/−60 OG1–C* projection is
−5.88965426 kcal/mol/A² and was rejected as an MM parameter; it does not
overturn the separate constrained-stability result. See the final section
of docs/ti1_qm_results_20261023.md. No further calculations were launched
after the user requested a quicker wrap-up. The next work is larger/native
reference preparation and validated parameter transfer, not protein repair.

## Completed Gaussian references and next gate — 2026-10-23

User resumed and requested results, documentation and next steps. All four
37atom capped TI1 optimizations, separate frequencies and conditional-curvature
assessments finished successfully. Each has105 positive raw frequencies and
102 positive constrained internal directions; no proton transfer or graph/
stereochemistry failure was detected. These are constrained surrogate minima,
not native enzyme TI1 or force-field validation. See
`docs/ti1_qm_results_20261023.md` for all four bond/frequency values and scope.

R cap sensitivity FAILS: OG1–C* shifts0.05327A and C*–O− shifts0.03166A,
above the predeclared0.03A threshold. S length sensitivity passes. Do not select
the passing configuration or average references to hide failure. Short capped
Thr N–H/oxyanion contacts persist. Larger/native reference geometry and native
charge/transfer validation remain necessary before protein TI1 refinement.

Next work: independent frozen-output audit/final inspection exports, diagnostic
raw-Hessian bond-stiffness comparison and separate capped-surrogate HF/MK ESP.
No native charge transfer, force-field fit or protein relaxation is certified.
New work is under `output_nyl12_parameter_validation_20261023/`; fresh surrogate
ESP uses separate `surrogate_esp_01` calculation directories when launched.
Read their process/acceptance records to distinguish preparation from execution.
The final standalone ESP checker passed27 tests/100%line coverage; full regression
passed447 tests with13 optional workstation-fixture skips, plus Ruff/format.
Final source-pinned verification: `output_nyl12_parameter_validation_20261023/full_verification_02/`.

The prior AF3-water discussion is now consolidated in
`docs/af3_water_placement_decision_20261023.md`: both minimal AF3 water tests
failed geometry acceptance; physical preconditioning removed crowding but the
subsequent DFTB3/MM diagnostic failed TI1 preservation. Better AF3 restraints
are useful parallel work but not the present prerequisite. Leave the user's
AF3 fork untouched; no new AF3 run or unvalidated protein refinement is started.

The active append-only journal remains `output_agent_actions/20261021.jsonl`.
Retain permission-review timeouts/bwrap failures and all failed attempts.

## Ranking and LigandMPNN symmetry clarification — 2026-10-22

New decision record: `docs/design_ranking_and_iteration_20261022.md`. The user
reports single-mutant confidence/strain degradation relative to WT, and expects
coordinated LigandMPNN redesign to improve fold/context compatibility. Preserve
the observation without treating confidence as affinity or all WT-relative
deformation as bad. Proposed ranking: chemistry gates, independent ES geometry
robustness and multistate/access compatibility, then local confidence/physical
quality/strain and noncovalent interface_delta as support. Sequence convergence
is diagnostic. No ranker, matched no-MSA arm or design campaign is implemented.
At the user's wrap-up clarification, primary geometry reproducibility refers to
unguided noncovalent ES. Covalent TI1 geometry and TS saddle validation are
separate gates. Interface_delta is optional pending demonstrated incremental
value/accuracy; it may be omitted. Leave existing QM tmux jobs undisturbed.

Official LigandMPNN README/run.py confirm automatic homooligomer mode ties all
parsed chains. Mature Nyl12 needs explicit separate alpha and beta ties with
equal weights and homo_oligomer0, as already used in the recent October16
protected preparation. Sequence ties do not select a geometric point group;
D2-like scaffold evidence and experimental assembly verification remain distinct.
Methods and the laptop inspection export are updated at this checkpoint;
scientific job state must still be read from current status/acceptance records.

## Redesign strategy/accessibility checkpoint — 2026-10-22

The user supplied the more recent `/data/jpc/enzyme-substrate-modeling/` source.
October16 focused redesign uses explicit separate alpha/beta ties, equal
weights and homo_oligomer0. Use its protected_inputs, not superseded root/pilot
snapshots; these were preparations, not an executed design campaign. Exact
source/provenance and caveats are recorded in the new strategy document.

The user's latest requirements and careful options comparison are recorded in
`docs/redesign_backbone_symmetry_and_access_20261022.md`. Recommended first
campaign: focused LigandMPNN on WT backbones/WT mature lengths, then independent
AF3 and relaxation. RFD3 is optional; if needed, initially remodel only selected
short pocket/entrance loops at exact WT lengths. Exact loop ranges are not yet
selected. Preserve catalytic/maturation/core/interface context. Indels remain
disabled until a per-design correspondence map and all postprocessing tests
cover them; WT-specific exporters currently hard-code232/96 and Amber indices.

The enzyme is (alpha-beta)4, not eight identical chains. Tie alpha copies and
beta copies independently. User has working LigandMPNN symmetry examples and
will point to them when needed. Do not assume C4 from tetramer stoichiometry:
a source-pinned tmux audit of the retained scaffold is D2-like (three near180deg
rotations, orthogonal axes, CA RMSD0.051–0.116A). Verify experimental assembly
operators before production RFD3. Raw audit is
output_nyl12_parameter_validation_20261022/repeat_symmetry_audit_01/.

Access screening must include open-compatible apo/ES and catalytic states,
actual substrate/product sterics, matched WT controls and state-aware sampling.
Sequence symmetry is distinct from identical dynamics/occupancy. AF3 confidence
and minimization do not establish opening kinetics; no opening in a short run
is unresolved. LigandMPNN multiple-PDB mode is batching, not joint ensemble
design. RFD3 exposure hints are static and do not guarantee channel access.
No RFD3/LigandMPNN campaign or new AF3 run has launched here. User AF3 fork is
untouched. Independent standard Amber typing/coverage preparation can proceed
without fitting charges or relaxing a protein with an unvalidated force field.

## Amber workflow clarification and parallel checks — 2026-10-22 evening

The user asked whether this uses standard Amber parameterization with
antechamber/Gaussian and ff19SB plus GAFF/GAFF2. Source audit confirms the
legacy scaffold uses original GAFF, not GAFF2, and provisional AM1-BCC or
Gasteiger charges. The new coupled TI1 Gaussian/RESP fit has not run.
Do not describe the existing topology as Gaussian-RESP parameterized.

New candidate plan: retain ff19SB for unmodified protein, use GAFF2 for
ordinary substrate chemistry, and fit separate HF/6-31G(d) MK ESP/RESP charges
with AmberTools utilities and explicit coupled-state/cap constraints. Use
antechamber for typing/file preparation and parmchk2 for parameter inventory;
validate their assignments. First test the standard candidate; refit only
missing or demonstrably inadequate reactive terms. Modified-Seminario is an
optional targeted tool, not an antechamber replacement or an automatic refit
of all bonds/angles. No legacy topology was migrated and no new FF was fitted.
Full detail: docs/ti1_parameter_protocol_20261022.md, final section.

While Hessians run, installed MCPB angle routines passed 96 manufactured
pair-spring/compliance checks plus three shared-angle tests, covering units,
AMBER harmonic convention, scaling, endpoint reversal and rotation. This is
API/convention verification, not recovery of an arbitrary physical molecular
Hessian or validation of TI1 parameters. Exact script, process log and report:
output_nyl12_parameter_validation_20261022/mcpb_angle_conventions_01/.
Compact versioned evidence: provenance/ti1_angle_conventions_20261022/.

The prospective independent validation plan is now
`docs/ti1_heldout_validation_plan_20261022.md`: 22 initial bond/angle distortions
per accepted reference, fixed baseline torsions, exact-coordinate QM/full-MM
comparisons and explicit training/validation separation. No scan inputs or
Gaussian scan jobs have been generated. Torsions and native-state transfer
still need separate coverage. Existing scientific gates and live tmux jobs
remain unchanged; user AF3 work remains untouched.

## Resumed checkpoint — 2026-10-22, 20:36 EDT

The user explicitly resumed; the previous 21:00 UTC stopping point is historical.
Two optimizations have now converged normally and passed the geometry gate:
`attack_S_cap_phi_minus60` and `attack_R_cap_phi_minus60`. Their fresh frequency
calculations are active in the existing tmux controllers. R/-120 and S/-120
remain unconverged at the time of inspection. No completed TI1 Hessian or
conditional-stability assessment has yet been accepted. Do not fit yet.

An independent frozen-log S/-60 audit reproduced its accepted geometry in
`output_nyl12_parameter_validation_20261022/completed_reference_audit_02/`.
Its `result/` contains assessment.json and final_structure.xyz/.sdf with all
37 atoms and 20 explicit H. Reactive bonds: OG1-C* 1.508577 A, C*-O- 1.302946 A,
C*-N 1.472865 A, C*-C 1.547196 A. No proton transfer detected. The capped Thr
backbone N-H contacts the oxyanion at 1.581748 A; retain this possible cap bias
for comparison against the other cap orientation. This is a constrained
bonded surrogate, not the native charged TI1; relay water is absent by design.

The installed pymsmt MCPB bond extractor passed eight manufactured quadratic
checks for bohr/Hartree units, one-based indexing, rotation, block averaging,
frequency-scale squaring and the AMBER K*(r-r0)^2 factor of one-half. Actual
TI1 fitting and modified-angle validation remain pending. Do not substitute
an arbitrary installed Seminario wrapper. Attempt01 emitted a complex-cast
reporting warning; attempt02 explicitly recorded real/imaginary components,
all with zero imaginary residuals. Both are retained. Independent geometry
audit01 failed before execution due to missing PYTHONPATH; fresh audit02 used
explicit project src/root paths and passed without changing chemistry gates.

Compact outputs, exact helpers and both failures are versioned under
`provenance/ti1_completed_reference_20261022/`. Production code/dependencies
are unchanged; the previous full-suite 433-test result still applies to that
source revision. New checks ran as scripted tmux jobs with source hashes.
Active action journal remains `output_agent_actions/20261021.jsonl`.

Next: finish the four geometry/Hessian/conditional-curvature gates, compare cap
hydrogen bonding and sensitivity, then fit and independently validate bonded
terms. Native 103-atom geometry and separate HF ESP/RESP/held-out validation
follow. AF3-model MM relaxation remains gated on the full TI1 force field,
then DFTB3/MM optimization. The user's AF3 fork remains untouched.

## Parallel preparation checkpoint — 2026-10-22, 20:36 UTC

The user visually inspected the four interim optimization trajectories and
judged them good. This is visual approval only: convergence, chemistry,
cap-artifact and native-transfer gates remain unchanged. Work resumed for
approximately one hour from 20:06 UTC; arrange a stopping point by 21:00 UTC.
Existing authorized Gaussian jobs continue in tmux during the user's absence.

At 20:35 UTC all four optimizers remain active and unconverged, around steps
11/11/11/12 in R/-120, R/-60, S/-120, S/-60 order. No completed TI1 Hessian,
ESP fit, accepted force field or new protein refinement is available.

Completed while waiting:

- Implemented conditional stability from the actual Gaussian gradient/Hessian,
  including the nonlinear constraint-force (Lagrangian) correction. Three
  frozen torsions plus six removed rigid motions leave 102 internal directions.
  Raw energy-Hessian projection or raw frequencies alone are insufficient.
  Both finite-difference steps must pass the prospective numerical criteria in
  docs/ti1_parameter_protocol_20261022.md. Missing gradients fail explicitly.
- Four tested analysis controllers now wait in detached tmux sessions
  `nyl12_ti1_curvature_<case>_20261022`. Their retained commands/statuses live
  under `<case>/conditional_curvature_controller_01/job/`, beside the existing
  frequency controllers. After successful frequency-controller completion,
  the pinned analyzer independently rechecks source/input identity, process
  acceptance, geometry and torsions, then writes a fresh
  `<case>/conditional_curvature_01/assessment.json`. Failed prerequisites,
  negative/unresolved curvature or nonstationarity cannot certify parameters.
  The source snapshot is `pa6_local_bonded/conditional_curvature_source/`.
  Waiters use a 48-hour deadline and 60-second polling, one BLAS thread, and
  launch no new Gaussian jobs. Original pinned frequency sources are unchanged.
- Prepared R/S native graph crosswalks (14 heavy-atom correspondences) and
  explicit proposed RESP constraints for the native 103-atom charge-zero
  Thr-NH3+/PA6-L4/oxyanion/NME/water graph. No native 3D geometry or ESP exists
  yet. No surrogate charge transfer or arbitrary residue integer charges.
- Independently sourced ff19SB/TIP3P with tleap and verified all 11 proposed
  fixed atoms against clean exported NTHR/NME/TP3 units. Peptide C/O are the
  NTHR values +0.6163/-0.5722. Native free-atom charge target is -0.0441 e.
  This is a source-consistent proposal, not an executed/accepted charge fit.
  Whole-library CHCL3BOX parser warnings remain in the first retained log;
  clean target-unit extraction/verification provides the source acceptance.

Verification: **433 full-suite tests passed**, Ruff and changed-file formatting
passed. New/modified module coverage is 84–96% (combined 89%). Separate waiter
safety tests passed 14/14, 94% coverage, without invoking engines. All checks
ran as scripted tmux jobs. Compact evidence, exact helpers and source hashes:
`provenance/ti1_parallel_preflight_20261022/`; larger originals are under
`output_nyl12_parameter_validation_20261022/`. The active append-only journal
remains `output_agent_actions/20261021.jsonl`. No dependency was added.

Next: inspect optimizer/frequency/conditional-curvature acceptance, then assess
cap H bonding and cap sensitivity before extracting bonded terms. Validate
bonded candidates with independent geometries/scans; build accepted native
charge-fit conformers and execute separate HF/6-31G(d) ESP/RESP/held-out checks.
Only after the full TI1 force-field gate passes may AF3-model classical repair
resume, followed by DFTB3/MM. The AF3 fork remains the user's independent work.

## Interim TI1 validation — 2026-10-22, 19:45 UTC

User requested validation and inspectable Gaussian outputs. All four optimizers
remain active around steps 6–7, with force criteria unmet; no frequency or final
parameter acceptance is available. Frozen live-log prefixes were audited in
tmux under `output_nyl12_parameter_validation_20261022/progress_audit_01/`.
Its `result/` contains 7/7/7/8-frame XYZ trajectories, latest XYZ/SDF, input
and atom maps, all-frame geometry checks, CSV metrics and PyMOL viewers.
Every frame has 37 atoms and 20 explicit H. Latest frames pass the declared
geometry/stereo/torsion/proton screens, but remain unconverged.

Both S cases show capped Thr backbone N–H/oxyanion contacts (1.77–1.93 A),
so cap H bonding needs explicit review before accepting any fitted reference.
The earlier S/-120 frame 4 had C*–O 1.225 A and failed the predeclared bond
screen; latest C*–O is 1.291 A. Retain all trials and unchanged final gates.
Do not interpret the latest-screen pass as convergence, a free minimum,
validated parameters or proof of negligible cap artifacts. Inspect S/-60
against R/-60 first. Relay water is intentionally absent from these small
bonded surrogates, not missing from the enzyme refinement model.

Actual PyMOL checks verified every exported state, 20 visible H, 36 intended
display bonds, coordinate agreement below 2.4e-7 A, and repeatable viewer reloads.
Checks: 371 full-suite tests passed; 34 focused export/converter tests passed,
97% packager coverage. Raw check outputs are `progress_viewer_check_03/` and
`inspection_full_verification_01/` under the same parameter-validation root.
The laptop package adds `qm/ti1_bonded_reference_progress/`; the existing
zero-argument pull script is unchanged. Original package revision 07 is retained
when the refreshed bundle is promoted. The next real calculation-stage gate
remains completed-optimization acceptance, then the existing automatic Hessian
continuation. Fitting/native-state transfer validation cannot use these interim
snapshots as final references.

## Active TI1 parameter checkpoint — 2026-10-22, afternoon

The user explicitly resumed after the earlier immediate stop. Four actual
Gaussian constrained optimizations are now running in detached tmux sessions
on violet, each with 8 cores and 16 GB. Inspect their retained status/logs;
this statement describes launch, not successful completion or fitted parameters.

The user selected `cif/n12c.cif` as the geometry reference. Author chain A,
Thr233 gives psi 134.110276428 degrees and chi1 -58.964417878 degrees;
native phi is undefined at this N-terminal Thr. The prescribed small bonded
surrogate is ACE-THR(-O-C*(O−)(CH3)(NHCH3))-NME, 37 atoms, charge -1,
singlet. Compare explicit ACE-cap phi choices -120/-60 degrees and both
attacked-carbon R/S configurations. Only phi, psi and chi1 are frozen;
reactive bonds and hydrogens are free. This capped anion does **not** supply
transferable charges for the native NH3+/oxyanion TI1.

All four deterministic starting geometries passed the predeclared geometry
screen. Prepared inputs, exact source/CIF snapshots, hashes and per-case jobs:
`output_nyl12_ti1_parameters_20261022/pa6_local_bonded/`.
Cases are `attack_R_cap_phi_minus120`, `attack_R_cap_phi_minus60`,
`attack_S_cap_phi_minus120`, `attack_S_cap_phi_minus60`; each actual run is
`optimization_01/`, containing `input.gjf`, `optimization.log`, retained
launcher and `job/status.json`. Tmux session names are
`nyl12_ti1_<case>_20261022`. Use the existing passwordless `sg g16 -c ...`
wrapper, explicit Gaussian profile and descriptive scratch directories.
No sudo or permission changes were made.

Optimization method is wB97XD/def2TZVP with UltraFine grid, Tight SCF/Opt,
CalcFC and XQC, no extra dispersion. Do not launch frequencies/ESP until
normal completion **and** TI1 protonation/connectivity/stereochemistry pass.
Four source-pinned frequency continuations also run in tmux as
`nyl12_ti1_freq_<case>_20261022`; their status/logs are under
`<case>/frequency_controller_01/job/`. They wait for the optimizer's final
process record, verify unchanged input hashes and apply the tested scientific
gate. Only an accepted reference proceeds to fresh `<case>/frequency_01/`,
with a separate Hessian checkpoint, formchk output and geometry cross-check.
Rejected references stop with retained acceptance records. There is no
automatic ESP fit, parameter certification or protein refinement.
Frequency and HF/6-31G(d) MK inputs use separate checkpoint names. Raw
frequencies of this constrained reference do not certify a free minimum;
conditional curvature, bonded fitting, held-out validation and native-state
RESP/transfer checks remain outstanding. The protein-refinement gate is closed.
See `docs/ti1_parameter_protocol_20261022.md` for the prospective acceptance
criteria and scope.

Latest full verification passed **368 tests**, Ruff and 87% combined coverage;
each new scientific module/entry point measured 82–99%. The separate
continuation safety tests passed **9 tests**, Ruff and 93% coverage without
invoking engines. Retained outputs are `full_verification_01/` and
`controller_verification_01/` under `output_nyl12_parameter_validation_20261022/`.
A real neutral-reference checkpoint conversion succeeded
in tmux under `output_nyl12_parameter_validation_20261022/real_fchk_probe_01/`;
the first parser check found a fixed-width character-array bug and was retained
as a failure. The corrected reader handles those records, then correctly
rejects this final-HF checkpoint because it has no Cartesian Hessian. No real TI1
Hessian has yet been accepted. Active append-only journal remains
`output_agent_actions/20261021.jsonl`; all original artifacts are preserved.

## Resumed 2026-10-22: viewer correction and TI1 validation

User is independently improving the AF3 contact-restraint fork and will report
when ready. Leave that fork untouched; no new AF3 tests until their update.
Water hydrogens exist in all 100 saved frames. Viewer scripts hid them, and a
semicolon in a comment caused the reported PyMOL SyntaxError. Corrected viewers
were tested twice in actual PyMOL 3.1.0: 100 trajectory states/one final state,
O/H1/H2 inventory, both O–H bonds, two visible water H, and no parser errors or
state duplication. Coordinate files are unchanged. Raw before/after reports and
PyMOL console logs are under `output_nyl12_viewer_fix_20261022/`.
Next scientific work remains coupled TI1 force-field validation. Source audit
confirms guessed junction terms, unfitted reactive charges and arbitrary O49
charge balancing; neutral free-NH2 Thr-NME QM artifacts are reusable provenance
but do not validate the covalent NH3+/oxyanion state.

The user clarified substrate termination on 2026-10-22: retain a real charged
PA6 terminus to support L2 release; a capped XA/MA fragment may be appropriate
for local copolyamide chemistry, recalling the Dishner JACS model. This concerns
substrate caps, separately from the peptide cap on Thr. See
`docs/substrate_termini_decision_20261022.md`. Preserve the existing full PA6 L4
reference and its terminal chemistry; source inspection confirms NH3+/COO−
termini. Do not silently cap it. The user subsequently supplied the original
capped copolyamide SMILES, preserved verbatim in the decision note. Scripted
RDKit/tmux audit confirms neutral acetyl-X-A|M-A-N-methylamide, C29H47N5O5;
slash cut C17/O36/N18 in zero-based original SMILES atom order. Original input,
script, report and job records are in
`output_nyl12_parameter_validation_20261022/substrate_termini_01/`.

Latest user clarification: PA6 L4 and L2 were both detected by MS, supporting
the hypothesis that soluble L4 can be a substrate yielding two L2 molecules.
Keep this first PA6 reference, but detection alone does not establish turnover.
Both charged termini represent a small oligomer; a long-polymer end model
would instead retain one real N- or C-terminus with a neutral distal boundary.
That is a separate future comparison. Only certain soluble low-molecular-weight
oligomers are detectable in the user's MS activity assay; do not infer the full
cleavage mechanism or total polymer conversion from that readout alone.
No coupled TI1 reference, QM calculation or AF3 refinement has been generated
in this audit. Next define its stabilization before generating QM fitting data.

## New scientific prerequisite — TI1 force-field validation

The user now requires validated TI1 parameters before classical AF3-model
relaxation, followed by DFTB3/MM. Read
`docs/ti1_forcefield_validation_gate_20261021.md` before any further refinement.
The current reactive topology is provisional and does not satisfy this gate.
The user suspects the poor initial geometry led the QM/MM structure toward a
reactant-like state; this is plausible but is not an accepted ES assignment.
The original failed run has been reproduced in a separate trajectory-only
replay:100 full-system frames, endpoint coordinate difference0.0 A.
No new model refinement was run after this prerequisite was established.

## Latest evening checkpoint — 2026-10-21

**Terminal update before disconnect:** attempt05 stopped after its first 100-cycle
QM/MM stage at 22:37 Eastern. No calculation from this work remains active.
Water contacts pass (maximum shortfall 0.029899 A), but OG1–C54 stretched to
2.195122 A, the directional water bridge failed, and RMS 0.63057 exceeded 0.5.
Later stages were withheld. This is a retained failed TI1 repair.
Diagnostic canonical-chain PDB:
`output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation/01_hydrogens_water/optimized.pdb`.
Next review state-preservation during early water relaxation and the DFTB3
Hamiltonian/starting geometry before defining another experiment. Do not rerun
unchanged or relax acceptance. The former running-state notes below are historical.

Read `docs/local_repair_checkpoint_20261021.md` first. The user authorized work
for about45 minutes after21:56 Eastern, with wrap-up around22:41. Active tmux
job: `nyl12_local_preconditioned_20261021`, attempt05 under
`output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/`.
Original AF3 oxygen is preserved during hydrogen construction. User-approved
TIP3P water-only pre-minimization completed with all nonwater atoms fixed exactly;
water O moved1.408215 A by minimization. It exhausted200 cycles above its gradient
target and is labeled coordinate preconditioning. DFTB3/MM subsequently started
with converged SCC evaluations. Inspect actual job/stage status on resume;
no accepted local minimum is yet claimed. Earlier attempts and failures retained.

## Active work — 2026-10-21, after evening resume

User authorized AF3 rough-water placement -> explicit hydrogens -> restrained
local DFTB3/MM repair. Begin with the two-restraint alpha/beta PA6 L4 model.
Both AF3 predictions completed; both independent water screens failed.
Attempt_01 has a 1.673 A water–C62 clash despite satisfying N/OG1 bounds;
attempt_02 with 28 additional exclusions breaks the N bridge (4.847 A).
See `docs/relay_water_repair_20261021.md` for the newly authorized separate
repair stage and water naming/Lennard–Jones requirements. Preparation is in
progress; no local minimization has been launched at this writing.
The earlier 18:00 pause was followed by explicit user authorization to resume.
Active journal remains `output_agent_actions/20261021.jsonl`.


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

## Model-switch handoff — 2026-10-21, evening

User wants GPT-6.1 Sol High for continuation. The rough-water restrained local
minimization remains authorized; no repair minimization has launched. Preserve
all current uncommitted work. New draft modules are `local_repair_model.py`
(explicit octamer topology subset -> 330-residue alpha/beta model; AF3 heavy
coordinate mapping; retained TIP3P parameters), `local_water_seed.py` (bounded
1.5 A oxygen displacement and water H orientation), and `local_hydrogens.py`
(local bonded-frame transfer from original topology). All under
`src/nyl12_covalent/`. These are development code, not accepted models.

Luna hydrogen helper agent added `tests/test_local_hydrogens.py`: initial run
reported 93% helper coverage but two incorrect expected-coordinate assertions;
agent corrected them, then could not rerun because tool approval timed out.
Water-seed test agent could not execute any command; no tests were added.
Root has not yet reviewed or tested these modules. Source/import integration,
terminal OXT handling, exact atom mapping, covalent inventory, QM boundary and
water parameters must be checked before use. No new dependencies were added.

Next: finish and test preparation, implement a separate local DFTB3/MM repair
runner and geometry audit, then run in tmux. Proposed stages: initial H/water
relaxation with protein/substrate heavy atoms restrained; local pocket relaxation
with outer protein restraints; then a separately recorded bridge-restraint
release. Use only two broad N/OG1–water distance restraints. Record flexible
selection and constants explicitly. These stages are not implemented yet.
Dimer QM mask is `:233,329,330|(:234&!@C,O)`, charge0, expected109 atoms,
sole link CA234–C234; verify before launch. Existing runner/audits use octamer
indices and must not be applied unchanged to this dimer. Do not silently relax
legacy scientific acceptance to force a launch.

Sources: topology/coordinates in
`/data/jpc/enzyme-substrate-modeling/output_nyl12_covalent_qmmm_20261019/pa6_l2_ti1_af3_neutral_beta/ti1.prmtop`
and sibling `ti1.inpcrd`; AF3 source is
`output_nyl12_af3_water_test_20261021/pa6_l4_alpha_beta/attempt_01/prediction/nyl12_pa6_l4_ti1_alpha_beta_water_contact/nyl12_pa6_l4_ti1_alpha_beta_water_contact_model.cif`.
Extract original residues1–232,929–1024,1313–1314 without changing parameters.
Sander: `/data/jpc/code/ambertools26/venv/bin/sander`.

AF3 audits are complete under
`output_nyl12_af3_water_audit_20261021/pa6_l4_alpha_beta/attempt_{01,02}/result/`.
Both failed water acceptance. Audit/export and optional exclusion-input code
and tests remain uncommitted from prior work. Laptop export has not yet been
refreshed with these diagnostic CIFs. Commit/push remains authorized at a
verified stopping point; do not claim it was completed.

Frequent automatic approval-review timeouts precede process creation; they are
not safety rejections. Default sandbox still fails bwrap RTM_NEWADDR. Record
failures honestly in `output_agent_actions/20261021.jsonl`; repair logging before
mutations. Neither unrelated tmux session was changed. Current tools may work
on a short single retry; do not repeatedly escalate unchanged rejected actions.

## Model 23 E-S checkpoint — 2026-10-26

Use the accepted static Boltz **Model 23** as the current PA6 L4-to-L2 E-S
starting pose. A fresh periodic system was built under
`qmmm/pa6_model23/00_system/` with ff19SB, TIP3P, and the exact historical PA6
ligand library/frcmod. The resulting system has 106,497 atoms, 106 NaCl pairs
at approximately 150 mM, and the four expected neutral beta-chain Thr termini.
The reviewed E-S QM region remaps exactly: 98 real QM atoms, eight link atoms,
106 QM nuclei, and charge -2.

There is **no completed PA6 L4 RESP/GAFF2 parameter set** in the audited
artifacts. The current ligand files are the historical GAFF1 set associated
with the earlier AM1-BCC workflow. They are acceptable for this coordinate and
DFTB3/MM geometry pilot, but must not be described as RESP/GAFF2. Generate and
validate a separate multi-conformer RESP/GAFF2 set before production classical
MD or TI1 classical MD.

The fully unrestrained classical negative control converged after 39,503 cycles
but degraded the catalytic pose (Tier B; Thr233 OG1 displacement 1.555 A and
psi change -64.707 degrees relative to the aligned reference). It is rejected
as a DFTB starting point and retained at
`qmmm/pa6_model23/clutter/01_full_unrestrained_classical_negative_control/`.

The accepted preparation instead relaxed LEaP-added solvent/ions and all
hydrogens while restraining only solute heavy atoms. It converged after 30,391
cycles and retained the Model 23 Tier A geometry: d1 2.884 A, d2 2.885 A, d3
2.871 A, aBD 94.490 degrees, and |phi| 99.041 degrees. The accepted restart is
`qmmm/pa6_model23/01_solvent_relaxation/solvent_relaxed_final.ncrst`; reports
are `geometry.json` and `alignment/alignment_report.json`. This scoped
solvation restraint is preparation only. The next QM/MM optimization has no
positional or reaction-coordinate restraints.

The fully unrestrained periodic DFTB3/MM E-S calculation is prepared, but was
**not launched before the planned reboot**, at
`qmmm/pa6_model23/02_dftb3_es_optimization/`. Its `minimize.mdin` explicitly
uses `ntr=0`, `nmropt=0`, and DL-FIND `active='@*'`; the exact command,
environment, hashes, and audited region are in `preparation.json`. After a
controlled reboot, first verify `nvidia-smi`, confirm no unrelated jobs, and
launch with `scripts/run_nyl12_job.py` in tmux. Although PMEMD.CUDA succeeded
when the correct Amber/CUDA environment was sourced, NVML currently reports a
driver/library mismatch, so the host should still be rebooted at this stopping
point.

The hASNaseIII GitHub repository was audited at commit
`8d1c076bc72393d4d7d794a8cf1d72bc92b3ff42`; see
`docs/hasnaseiii_repository_audit_20261026.md`. It provides no runnable scripts
or license, but its Amber/Gaussian six-CV string inputs are directly useful as
a mechanistic design reference. They argue for a multidimensional Nyl12 path
that couples proton transfer, OG1 attack, carbonyl rehybridization, and C-N
cleavage rather than further one-dimensional coordinate driving. Do not copy
the repository inputs; implement and cite an Nyl12-specific version.

## Active Model 23 DFTB3/MM E-S optimization — 2026-10-27

Violet was rebooted successfully. NVIDIA driver/NVML 580.196.07 is healthy and
both RTX 6000 Ada GPUs were idle at preflight. All recorded preparation hashes
matched, and the QM-region audit remained 98 real atoms plus eight link atoms
(106 QM nuclei), charge -2. The calculation is fully unrestrained (`ntr=0`,
`nmropt=0`, DL-FIND `active='@*'`).

The optimization is running in tmux session
`nyl12_model23_dftb_es_20261027`. Calculation files are under
`qmmm/pa6_model23/02_dftb3_es_optimization/`; launcher provenance and status
are under `provenance/attempt_01/`. Initial DFTB3/MM evaluation and at least 30
DL-FIND cycles completed without initialization errors, and `RESTRAINT=0.0`.
Trajectory frames are written every 25 optimizer cycles to `optimization.nc`.
Do not infer convergence or chemical acceptance from the live restart. After
process completion, inspect `status.json` and `minimized.mdout`, verify DL-FIND
convergence, regenerate catalytic-geometry/alignment reports, and visually
inspect the complete optimization trajectory before accepting the endpoint.

Live inspection after more than 3,000 cycles found active-site degradation;
see `qmmm/pa6_model23/02_dftb3_es_optimization/diagnostics/`. The frozen
11:27:49 EDT snapshot is Tier B rather than Tier A: d1/d2/d3 =
3.345/3.267/3.065 A, aBD = 99.788 degrees, and projected phi = 123.800
degrees. Thr233 heavy-atom RMSD to the fitted `n12c.cif` scaffold increased
from 0.623 to 1.355 A, OG1 displacement from 0.189 to 1.618 A, and psi
difference from -30.01 to -60.14 degrees. This closely resembles the rejected
fully unrestrained classical basin. Current strict nonconvergence is dominated
by individual bulk-water hydrogens while the RMS gradient already passes.
Do not accept this live endpoint merely if DL-FIND later reports convergence.

The active-site measurement helper was corrected to accept explicitly solvated
topologies with the validated 1,314-residue solute prefix; it previously
required exactly 1,314 total residues and failed on the 30,303-residue system.
The fix is covered by a solvated-topology regression test.

At the 2026-10-27 12:17 EDT snapshot after more than 5,100 cycles, degradation
continued: d1/d2/d3 = 3.407/3.092/3.019 A, aBD = 99.030 degrees, and projected
phi = 123.857 degrees. Thr233 heavy-atom RMSD is 1.424 A, OG1 displacement
1.811 A, N displacement 2.255 A, chi1 difference +25.52 degrees, and psi
difference -57.92 degrees. The structure remains only narrowly Tier B and must
not be accepted as the preferred E-S reference.

The newly implemented CB–OG1–C54 metric is 119.085 degrees in the accepted
Model 23 start, 111.642 degrees near cycle 3,200, and 109.659 degrees near cycle
5,100. Andjelkovic et al.'s deposited hASNaseIII TS1 and TS2 structures give
135.211 degrees (OG1–C 3.365 A) and 116.949 degrees (2.470 A), respectively.
Thus the angle is useful but not sufficient: degraded Nyl12 geometries can
retain a plausible value. The official SI provides no Michaelis-complex
coordinates. See `docs/hasnaseiii_repository_audit_20261026.md`.

The next recommended Nyl12 step is the bounded pilot in
`docs/template_conditioned_es_plan_20261027.md`: use `n12c.cif` as a tetramer
template, generate PA6 L4-to-L2 models with minimal register/site contacts,
release those contacts before acceptance, and use staged MM plus inner-shell
DFTB3/MM refinement. Do not repeat whole-system strict minimization dominated
by bulk-water maximum components. Nyl12 remains justified by experimental
activity, observed PA6 products, its crystal structure, and stable prior PA6
MD; NylB should be a lower-complexity control if the bounded pilot fails, not an
immediate replacement.

## Template-conditioned Models 7 and 8 Amber rescue — 2026-10-27

The user authorized a bounded Amber rescue of the two raw Tier B seeds, Models
7 and 8. Both completed successfully in independent provenance-tracked tmux
jobs. The diagnostic workflow was fully unrestrained and unsolvated OpenMM
minimization with ff19SB plus the retained historical PA6 L4 AM1-BCC/GAFF1
parameters; it must not be described as RESP/GAFF2 or production MD.

Model 7 changed from d1/d2/d3 = 2.603/2.887/2.818 A, aBD = 89.837 degrees,
CB-OG1-C54 = 121.676 degrees and phi = -115.312 degrees (Tier B) to
3.347/2.862/2.871 A, 85.022 degrees, 107.073 degrees and -123.167 degrees
(still Tier B). Model 8 changed from 2.532/2.800/3.002 A, 87.409 degrees,
118.400 degrees and -106.455 degrees (Tier B) to 3.096/3.055/2.971 A,
60.135 degrees, 121.414 degrees and -118.938 degrees (reject). Model 8's
initial potential energy was 3.55 million kJ/mol, indicating severe clashes
after force-field system building; it should not be promoted. Model 7 is also
not promoted because it did not reach Tier A. Full outputs and provenance are
under `models/nyl12_template_conditioned_es/02_amber_rescue_relaxation/`.

## Active template-conditioned Model 7 DFTB3/MM pilot — 2026-10-27

The Model 7 periodic system was rebuilt from its Amber/OpenMM endpoint with
ff19SB, historical PA6 L4 AM1-BCC/GAFF1, TIP3P and 150 mM salt. Direct reuse of
the earlier unsolvated topology was rejected because only one of four beta
Thr N-termini had been neutralized there, changing the validated atom mapping.
The periodic rebuild has all four neutral catalytic Thr termini.

New solvent, ions and hydrogens were relaxed with solute heavy atoms restrained.
The first two 10,000-cycle GPU segments reached their cycle limits; the third
continuation converged after 45 additional cycles (RMS = 0.00984). The
pre-QM/MM active-site geometry remained Tier B: d1/d2/d3 =
3.347/2.861/2.869 A, aBD = 85.229 degrees, CB-OG1-C54 = 107.174 degrees and
phi = -123.143 degrees.

The established QM core was retained rather than enlarged: 98 real atoms,
eight link atoms, charge -2, comprising Thr233/Thr234, the Phe113 backbone
donor, Asn186 side chain, Asp272, Asp274 and the local reacting PA6 amide. The
prior Model 23 result does not establish a missing QM atom or QM water. The
Model 7 optimizer instead corrects the previous excessive-coordinate problem:
all protein and substrate atoms plus complete waters within 8 A of the
catalytic QM-component residues are movable (19,976 atoms; three waters in the
starting structure), while only distant bulk solvent is fixed. There are no
energetic or reaction-coordinate restraints (`ntr=0`, `nmropt=0`, observed
`RESTRAINT=0.0`).

The calculation is active in tmux session
`nyl12_template_model7_dftb_es_20261027`. Files are under
`qmmm/pa6_template_model7/02_dftb3_es_optimization/`; exact hashes, command and
environment are under `provenance/attempt_01/`. Initialization succeeded and
the optimization reached at least cycle 57. Do not accept the endpoint until
process completion, convergence checks, catalytic geometry/alignment analysis,
and trajectory inspection.
## CURRENT — PA6 L4 zwitterion RESP/GAFF2 generation active, 2026-10-27

The user authorized production-quality PA6 L4 parameters before enzyme MD.
The locked noncovalent substrate graph is the full 79-atom zwitterion with
N56 NH3+, C53/O48/O51 COO-, net charge zero and singlet multiplicity. This is
the primary soluble, solvent-exposed oligomer state near neutral pH; it is not
a claim that termini buried in crystalline nylon are always charged. A neutral
NH2/COOH sensitivity state, if needed, requires a separate parameter set.

The reproducible workflow is under
`parameters/pa6_l4_zwitterion_resp_gaff2/`. Four bound-like conformers comprise
template-conditioned Model 7 plus one snapshot from each of three independent
historical 100-ns WT trajectories. Inputs preserve canonical atom order and
freeze 25 declared heavy-atom backbone/amide torsions during
omegaB97X-D/def2-TZVP/SMD(water) optimization; all other coordinates move.
Focused preparation/graph tests pass (13/13), and Ruff/diff checks pass.

Four eight-core Gaussian 16 jobs are active in tmux sessions named
`pa6_l4_resp_dft_<conformer>_a2_20261027`. Exact inputs, hashes, scripts and
status are retained under each conformer's `job_control_attempt2/`. An initial
wrapper attempt failed before Gaussian because the vendor profile referenced
unset `PERLLIB` under Bash nounset; failed `job_control/` records remain and no
scientific steps were produced. The narrow wrapper fix was preflighted before
attempt 2. After geometry acceptance, run fixed-geometry gas-phase
HF/6-31G(d) MK ESP calculations, a shared multi-conformer two-stage RESP fit,
GAFF2 typing/parmchk2 review, nylon torsion validation, and isolated/solvated
validation. Historical GAFF1/AM1-BCC parameters remain a comparison control.

## CURRENT — PA6 RESP geometry jobs restarted without CalcFC, 2026-10-27

The four omegaB97X-D/def2-TZVP/SMD jobs in
`parameters/pa6_l4_zwitterion_resp_gaff2/01_dft_geometry/` were intentionally
interrupted after more than three hours. All initial SCFs had converged, but
none had completed the full initial Hessian requested by `CalcFC`, and no
geometry step had begun. The exact tmux sessions exited with status 130; logs,
checkpoints, hashes and both earlier job-control attempts remain. No process
from those attempts is active.

The replacement inputs are in `02_dft_geometry_def2svp/` and use
omegaB97X-D/def2-SVP, SMD water, 25 frozen heavy-atom backbone/amide torsions,
`Opt=(Tight,ModRedundant,MaxCycles=300)`, eight cores and 16 GB, explicitly
without `CalcFC`. Four new provenance-tracked tmux sessions named
`pa6_l4_resp_def2svp_<conformer>_20261027` are active. Def2-SVP is the geometry
basis only; fixed accepted geometries will still receive gas-phase
HF/6-31G(d) Merz-Kollman ESP calculations. Before certification, compare one
accepted conformer against def2-TZVP to quantify geometry sensitivity. The
restart generator and metadata-driven acceptance checks have focused test and
Ruff coverage; the production inputs independently pass atom-count, route,
constraint-count and no-`CalcFC` checks.

## CURRENT — PA6 L4 RESP fit complete; GAFF2 validation active, 2026-10-28

All four omegaB97X-D/def2-SVP/SMD constrained geometries completed normally
and passed graph, electronic-state, frozen-torsion and terminal-proton
acceptance. Their fixed-geometry gas-phase HF/6-31G(d) MK calculations also
completed normally. AmberTools 26 `espgen` outputs were validated and combined
in the declared order: template Model 7, then historical MD replicas 1--3.

The shared four-conformer two-stage RESP fit completed under
`parameters/pa6_l4_zwitterion_resp_gaff2/03_multiconformer_resp/results_attempt3/`.
The stage-2 ESP relative RMS is 0.03133; corresponding atom charges are
identical across conformers and the rounded 79-atom charge sum is -0.000004 e.
GAFF2 typing and `parmchk2` completed. There are no missing mass, bond, angle or
nonbonded entries. Seven analogous dihedrals and three impropers were emitted;
the improper with penalty 6.0 and the default improper require explicit review,
and nylon torsion scans remain mandatory. The candidate is not MD-certified.

The first minimization validation used the pre-QM MOL2 coordinates and exposed
a 1.196 A H33--H36 clash; its nonconverged artifacts are retained as failed
scientific evidence. Attempt 2 transfers the accepted template-Model-7 DFT
coordinates without changing RESP charges or GAFF2 types. LEaP reports zero
errors and zero warnings. The implicit-solvent XMIN continuation converged at
RMS 0.00441 kcal mol-1 A-1. Explicit TIP3P minimization and a staged 170-ps
sanity test completed. During the final unrestrained 100 ps, mean ligand
heavy-atom RMSD was 2.58 A, the former H33--H36 contact remained at least
2.90 A, and no amide isomerization occurred. Mean amide deviations from
planarity were 7.2--8.6 degrees, with rare excursions near 30 degrees. Exact
metrics and hashes are in `04_gaff2_validation/attempt2/short_md_summary.json`.
The candidate passes short-MD sanity but remains uncertified until the key
nylon torsions and flagged `parmchk2` impropers are benchmarked against QM.

## CURRENT — integrated PA6 L4 torsion recovery active, 2026-10-29

Six of eight bidirectional wB97X-D/def2-SVP/SMD(water) bonded benchmarks have
completed. The carbonyl-side forward scan completed all 25 points; its reverse
scan remains active in tmux session
`pa6_scan_internal_neutral_carbonyl_side_20261028`. The original amide-omega
forward scan stopped after ten optimized points with an L9999 error
termination, before the reverse direction began. Those failed artifacts were
preserved unchanged.

A tested, non-overwriting recovery was generated under
`parameters/pa6_l4_zwitterion_resp_gaff2/05_qm_torsion_benchmarks/02_relaxed_scans/internal_neutral/amide_omega/recovery_20261029/`.
It continues forward from the last accepted geometry using 45 increments of
5 degrees (46 optimized points including the starting geometry) and starts a
fresh reverse scan from the original seed using 72 increments of -5 degrees
(73 points). Both use `CalcFC`, `MaxCycles=300`, separate checkpoints/scratch,
and provenance-recorded tmux sessions
`pa6_amide_omega_recovery_forward_20261029` and
`pa6_amide_omega_recovery_reverse_20261029`. Do not analyze or refit until all
active scans terminate normally and the original/recovery forward segments are
merged without double-counting their shared endpoint.

## HOST SAFETY — azalea Gaussian installation, 2026-11-06

Gaussian 16 is already installed on azalea at `/usr/local/g16`; `jpc` belongs
to the `g16` group, the vendor profile resolves the executable, and the g16
SHA-256 is `9dd1b1a7495c313954b5243ea427389adac55e55829ab4934f89b8e46fb0c8d5`.
The ISO is not required for normal calculations. `/etc/fstab` currently has a
required loop-mount entry from missing `/mnt/CMBStorage/Gaussian16.iso` to
`/mnt/gaussian`; this caused `mnt-gaussian.mount` to fail and systemd to enter
maintenance/emergency behavior at boot. `/mnt/CMBStorage` is a directory on
the root filesystem, not a mounted storage volume.

Do not mount the missing ISO, edit `/etc/fstab`, reboot azalea, or make the ISO
a workflow prerequisite. Any boot-configuration repair requires a separate,
explicit user authorization. Gaussian jobs should use the installed
`/usr/local/g16` tree and established `g16` group wrapper only.

## CURRENT — RESP/GAFF2 WT QM/MM restart, 2026-11-04

The active, source-pinned checkpoint is
`docs/qmmm_resp_gaff2_restart_20261104.md`. WT modular PA6 MD replicas contain
many joint Tier A/B plus native-Thr-torsion frames. Replica 1 frame 1553 is the
first near-pass seed. Fully unrestrained MM minimization is queued behind the
current violet GPU job, with conditional 8 and 12 A movable-shell DFTB3/MM
branches. They fail closed if MM convergence, Tier A/B geometry or native
Thr233 gates fail. The modular QM region has been rebuilt and audited as 98
real atoms + 8 links, charge -2; old whole-ligand indices are not reused.

The prior unbiased TI1 validation did revert to E--S under the old setup.
DFTB3 proton-transfer accuracy remains unvalidated for this chemistry; the
small gas-phase omegaB97X-D/def2-SVP comparison remains required. Do not begin
new coordinate driving until an accepted E--S endpoint and Hamiltonian gate are
available.

The frame-1553 MM continuation completed on 2026-11-05. A 20,000-cycle budget
converged after 9,317 additional cycles (RMS 0.009969), but the endpoint stayed
Tier B and failed the native Thr233 OG1-displacement gate (0.871 A versus 0.75
A). CB--OG1--C stabilized near 120.26 degrees after falling from 141.28 degrees
in the MD snapshot. The conditional QM/MM branches were not relaunched. See
`docs/qmmm_resp_gaff2_restart_20261104.md` and the immutable continuation
acceptance record for the exact decision.

The full three-replica native-Thr screen completed on 2026-11-06. Replica 3
contributed 157 joint Tier A/B plus native-torsion frames but no displacement-
gate pass. Across all replicas, only four frames passed, all in replica 1 from
77.65--80.10 ns. Frames 1554, 1562, and 1603 are the declared next MM seeds;
frame 1572 is retained as an alternate. See
`qmmm/pa6_resp_gaff2_wt/07_selected_es_candidates_20261106/selection_decision.json`
for exact metrics, limitations, and input hashes. These correlated seeds must
still pass fully unrestrained MM endpoint acceptance before any DFTB3/MM launch.

Frames 1554, 1562, and 1603 were subsequently extracted and subjected to
fully unrestrained segmented MM minimization. All converged in the second
10,000-cycle segment, but none passed the complete endpoint gate. Frames 1554
and 1603 remained Tier B and missed only the native Thr233 OG1-displacement
limit (0.830 and 0.782 A, respectively, versus 0.750 A); frame 1562 also lost
Tier-B geometry (`d1 = 3.576 A`). No DFTB3/MM job launched. Frame 1603 is the
declared near-pass; do not loosen the gate post hoc without an explicit
scientific decision. See `docs/qmmm_resp_gaff2_restart_20261104.md` and
`qmmm/pa6_resp_gaff2_wt/09_mm_qmmm_selected_candidates_20261106/`.

The user subsequently authorized a documented frame-1603 sensitivity test.
Only the starting OG1-displacement ceiling is explicitly widened to 0.800 A;
all other gates remain unchanged, and standard acceptance remains recorded as
false. Independent 8-A and 12-A DFTB3/MM optimizations are active in tmux as
`nyl12_f1603_nearpass_qmmm08_20261106` and
`nyl12_f1603_nearpass_qmmm12_20261106`. Artifacts and authorization are under
`qmmm/pa6_resp_gaff2_wt/10_frame1603_nearpass_qmmm_20261106/`. Final endpoint
promotion still requires the standard, unmodified QM/MM acceptance criteria.

## CURRENT — window-24 hybrid-DFT/MM basin diagnostic, 2026-11-09

The chemically correct endpoint of asynchronous window 24 was preserved and
hashed before launching a hybrid-DFT/MM diagnostic. It has `Thr233:HG1` on the
N-terminal nitrogen, but `Thr233:OG1`--substrate-carbonyl-C is still 2.104 Å;
therefore it is a proton-transferred pre-TI1 structure, not a covalent TI1.

The calculation is under
`qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths/02_asynchronous_es_to_ti1/03_window24_wb97xd_def2svp_mm/`.
It retains the 98-real-atom QM region (+8 links), charge -2, and the established
12 Å movable shell, with no reaction-coordinate or Cartesian restraints.
Gaussian 16 uses ωB97X-D/def2-SVP, an ultrafine grid, tight/XQC SCF, 16
threads, and 32 GB. No additional dispersion correction is applied.

Amber's external-Gaussian interface rejects QM Ewald/PME. The active diagnostic
therefore uses periodically imaged MM point charges within 20 Å of any QM atom
(16,000+ charge sites in the first generated input), while retaining periodic
Amber MM--MM treatment. This is adequate to test local basin behavior but is
not yet a quantitatively validated DFT/MM energy protocol. Require a 15/20/25 Å
embedding-cutoff sensitivity test before interpreting relative energies.

The persistent session is
`nyl12_w24_wb97xd_def2svp_mm_20261109`. The wrapper first runs one Gaussian
energy/gradient and checks both Gaussian and Sander completion; only then does
it start the 500-cycle DL-FIND/L-BFGS optimization. At this checkpoint the
first ωB97X-D SCF was active. Two earlier launches terminated before any
calculation: first because the vendor Gaussian profile conflicted with strict
shell nounset handling, then because the tmux shell lacked Amber's runtime
library path. Both launcher defects were corrected and retained in the action
journal; neither produced a scientific result.

Interpretation must remain limited. Reversion to E--S would reject a stable
proton-transferred minimum in the window-24 basin at this DFT/MM setup, but it
would not rule out a TI1 minimum starting from an OG1--C geometry near
1.5--1.6 Å.

## CURRENT — MBAR pipeline smoke tests, 2026-11-09

A typed PyMBAR 4 analysis path now performs per-window statistical-
inefficiency estimation and decorrelation, MBAR overlap/effective-sample
diagnostics, and a histogram FES with analytical uncertainties. It uses the
Amber restraint convention `k*(xi-target)^2`; the restraint bias is therefore
included consistently in the reduced potentials. Outputs are explicitly
watermarked and labeled as pipeline smoke tests, not converged or promotable
PMFs.

Focused tests pass (5/5; 87% coverage of the new analysis module). The
completed direct-proton-transfer 2-ps data set produced 8,000 raw and 4,268
decorrelated samples. Its MBAR overlap scalar is 0.0137 and minimum
bidirectional adjacent-window overlap is 0.164. The immutable asynchronous
block-1 snapshot taken at 20:03 EST contains 205 samples in each of 25 windows
(5,125 raw; 2,591 decorrelated); its corresponding diagnostics are 0.00367
and 0.0791. These short profiles are useful only to validate the pipeline and
identify inadequate sampling/overlap.

The asynchronous campaign remains active on azalea in tmux session
`nyl12_async_umbrella_30ps_20261109`. At 20:03 EST all 25 windows were at
2,900--3,300 of 10,000 block-1 steps (mean 3,224; 1.612 of 5 ps). The observed
rate predicts block-1 completion near 22:00 EST, subject to workstation load.

An updated synchronized snapshot was taken at 21:09 EST with 615 reaction-
coordinate samples per window (15,375 raw; 7,236 decorrelated). The provisional
MBAR overlap scalar is 0.00469 and the minimum bidirectional adjacent-window
overlap is 0.100. The profile reaches about 19.7 kcal/mol before declining by
roughly 1 kcal/mol near the product-side boundary around -0.7 A. This endpoint
downturn must not be interpreted as a TI1 basin. A window-24 trajectory audit
showed that all four saved heating frames retained HG1, H1 and H2 on the
N-terminal nitrogen with OG1--C = 2.096--2.197 A. H1 then transferred from N
to OG1 before the first saved equilibration frame. All 20 equilibration and 80
production frames have HG1 and H2 on N, H1 on OG1, and OG1--C = 2.971--3.690
A; none retains the intended three-N--H protonation or OG1--C <= 2.5 A. The
composite coordinate remains near its target because the increases in
OG1--HG1 and OG1--C compensate. Product-side sampling must therefore be rebuilt
from the last valid heating restart using explicit N--H proton-state guards and
OG1--C as the attack coordinate, or with a proper two-dimensional PMF.

## Stepwise attack publication recovery, 2026-11-11

02_nucleophilic_attack stopped after accepted attack_01 solely because the runner
exported update_0002 twice with plots enabled. No second scientific window had
been attempted. Acceptance gates native convergence, all proton identities and
scissile connectivity passed. Runner now publishes once and can continue from an
accepted checkpoint without promoting or recomputing it. Ruff and 26 focused
tests passed; runner coverage 84%. Original snapshot/logs remain intact; stopped
status and input metadata copied to publication_recovery/. New immutable snapshot
source_snapshot_publication_recovery/ and job_publication_recovery/ retain hashes.
Tmux model30_neutral_attack_publication_recovery_20261111 resumed the original
authorized OG1-C protocol, original start reference, step size 0.10 A maximum,
endpoint 1.50 A, 2000 cycles and final reaction-unrestrained release. Actual
attack_02 process verified starts from accepted attack_01 with exactly one
OG1-C restraint and no N-H restraint. No duplicated movie frames or publication
from resume itself. Scientific outcome remains pending; inspect status and gates.

## Flipped-model direct-acyl TS search, 2026-11-13

The original DL-FIND dimer search was seeded from the final-cycle-499 NEB
highest-energy image 8, with an initial tangent from images 7 and 9. The
SCC-guarded `dftb_telec=300` K diagnostic continued that same initial NEB
seed for 60 dimer cycles and passed the prior cycle-50 SCC failure. It did not
converge a saddle: maximum gradient was 0.040982 against 0.0001, and the saved
geometry was product-like (Thr233 OG1--substrate C 1.401 A; HG1--leaving N
1.024 A; scissile C--N 2.107 A). Its energy is not a validated barrier.

The 200-cycle SCC-guarded continuation in tmux session
`model30_dimer_scc300_continue_20261113` starts from that cycle-60 endpoint,
not directly from the NEB image. The dimer direction was reconstructed from
the saved mode trajectory; the 12 A movable region (4,165 atoms), 300 K
electronic temperature, 1e-10 SCC tolerance, and dimer algorithm are retained.
This is a new DL-FIND search from saved coordinates and mode; prior L-BFGS
history is not carried over. Inputs, hashes, status, and native output are in
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/06_dimer_ts/scc300_pilot_20261113/continuation_inputs_cycle60_20261113/`.
The guard stops at the first SCC failure. This continuation instead stopped
near dimer cycle 112 with Sander exit -11 (SIGSEGV) inside PME reciprocal-space
`fill_charge_grid`, without an SCC failure. At cycle 112 the maximum gradient
was 0.034233 versus a 0.0001 target, the step was essentially zero, and the
reported dimer curvature was positive (0.00121). It was not a converged TS.
The saved restart contains finite coordinates; its largest atomic displacement
from the continuation start is 0.870 A. Both start and saved restart have the
same 5,486 atoms outside the primary periodic box, so the crash cannot yet be
attributed to a newly displaced atom crossing the box. The precise PME crash
trigger is unknown. No automatic restart was launched; native output is kept.

Keep the dimer search as the immediate diagnostic. Product-like movement may
indicate a saddle toward acylenzyme or escape from the intended saddle; more
cycles alone do not establish a valid TS. The installed Sander/DL-FIND exposes
P-RFO with Bofill Hessian updates, but a full Hessian over the unchanged 12 A
region would have 12,495 Cartesian coordinates and be costly. A GPR Hessian
update is not exposed by this interface, so it is not a drop-in keyword change.
If dimer stalls or finds an unwanted mode, revisit a bracketed NEB image and
reactive mode before a bounded alternative search. A candidate needs one
imaginary mode, endpoint connectivity, and a 0 K DFTB3 recheck; 300 K smearing
is an SCC diagnostic, not final TS validation.

### Bracketed image-7/image-8 TS reseed, 2026-11-13

The 300 K image-8-seeded dimer continuation stalled by cycle 100 and then
crashed in PME near cycle 112; it must not be resumed as a TS candidate. Native
NEB cycle 499 was not converged (energy residual 5.2219e-4 vs 1e-5) and had no
climbing image. Image 7 was reactant-like (OG1--C 2.045 A, HG1--leaving N
1.836 A); the highest-energy image 8 was already product-like (1.474 and
1.047 A). A fresh seed at the image-7/image-8 midpoint gives OG1--C 1.7595 A,
OG1--HG1 1.4571 A, HG1--leaving N 1.4370 A, and scissile C--N 1.6940 A.
The local audit found no nonreactive heavy-atom contacts below 2.0 A. Its
initial dimer direction is the local 7-to-8 tangent, displaced 0.01 A on the
rotation mask. The seed preparation is scripted, hashed and tested (5 focused
tests passed, 84% helper coverage, Ruff passed).

The new DFTB3/MM dimer search completed normally at its 200-cycle limit in
tmux session `model30_dimer_bracket78_20261113`; artifacts are under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/`.
It retains the 12 A movable region and original QM/MM settings, with 300 K
electronic smearing and a 200-cycle SCC-guarded budget. Sander exited 0 with
no SCC failures or crash, but DL-FIND did not converge: final max gradient
3.4118e-4 vs 1e-4, max step 6.0723e-3 vs 4e-4. Dimer curvature remained
negative (-0.12285). The final reactive geometry shifted toward a proton-
transfer saddle rather than concerted acylation: OG1--C 2.0678 A, OG1--HG1
1.4204 A, HG1--leaving N 1.2110 A, scissile C--N 1.5651 A. It is not a
validated TS or barrier.

The 200-cycle continuation in tmux session
`model30_dimer_bracket78_continue_20261113` finished normally inside the same
07_dimer_ts directory at `continuation_cycle200_20261113/`. The saved cycle-200
restart and dimer mode were reconstructed with source hashes; initial native
energy matched the prior endpoint (-718.1334165 Eh). DL-FIND **converged at
cycle 159** with all five criteria passing: max gradient 7.1496e-5 (target
1e-4), max step 3.0105e-4 (target 4e-4), and final energy -718.134029878265
Eh. Dimer curvature was negative (-0.12048); no SCC failure occurred. Final
geometry: OG1--C 2.0851 A, OG1--HG1 1.4414 A, HG1--leaving N 1.2183 A,
scissile C--N 1.5376 A. Displacement along the saved dimer mode couples HG1
transfer to leaving N, OG1--C shortening and C--N lengthening, but a frequency
calculation and downhill endpoint connectivity have not yet been performed.
This is a converged dimer stationary-point candidate at 300 K electronic
smearing, **not yet a validated first-order TS or activation barrier**. Do not
restart for more optimizer cycles. Next: zero-smearing evaluation, mode and
endpoint checks; also retain the possibility that this saddle represents an
early proton-transfer-dominated step rather than the full concerted reaction.

### TS validation progress, 2026-11-13

The combined, stripped VMD optimization movie has 71 frames and is at
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/vmd_progress/latest/view_ts_optimization.vmd`.
The laptop pull script includes this bundle. A native DFTB3/MM one-cycle
stationarity check at 0 K electronic temperature converged immediately at
-718.134029878173 Eh, matching the 300 K dimer endpoint to the reported
precision. It did not calculate a fresh Hessian or imaginary frequency.
The original optimized E-S energy was -718.2079132 Eh, so the candidate is
+46.3625 kcal/mol on this potential-energy comparison; this is not a
validated kinetic or free-energy barrier.

The Python Sander API failed on this ff19SB CMAP topology. A native Sander
one-cycle analytic-gradient adapter passed a source-audited force conversion
and a 0 K finite-difference directional-curvature check (2.03% difference).
Its acceptance record is under
`continuation_cycle200_20261113/native_gradient_0K_onecycle_20261113_attempt2/scientific_acceptance.json`.
Two numerical mass-weighted downhill branches from the converged dimer mode
are being run in tmux to test endpoint connectivity. This is a provisional
descent pilot, **not a formal Hessian-validated IRC**. The first launch under
`irc_0K/job` stopped before an energy evaluation because `AMBERHOME` was
missing; its evidence is retained. The corrected launch uses
`AMBERHOME=/data/jpc/code/ambertools26/venv`, tmux session
`model30_ts_irc0k_attempt2_20261113`, and the same scientific inputs under
`irc_0K/job_attempt2` and `irc_0K/paths_attempt2`. Endpoint chemistry and
first-order saddle character remain to be checked.

The corrected numerical downhill pilot subsequently finished normally with
20 accepted steps in each direction. The mode-plus endpoint is -718.1348814
Eh (OG1--substrate C 2.076 A, HG1--leaving N 1.163 A); mode-minus is
-718.1374740 Eh (OG1--C 2.094 A, HG1--leaving N 1.324 A). Both are still
near the starting OG1--C separation. Thus these short paths demonstrate
downhill proton-transfer movement, not connectivity to E-S and acylenzyme.
Their source-native evaluation files and endpoint records are in
`irc_0K/paths_attempt2/`.

Two unrestrained 0 K DFTB3/MM endpoint minimizations, each allowing the same
4,165-atom 12 A active region to move, started from those branch checkpoints
in tmux sessions `model30_ts_mode_plus_min_20261113` and
`model30_ts_mode_minus_min_20261113`. Inputs, trajectories and launcher status
are in `irc_0K/endpoint_minimizations/mode_{plus,minus}/`. A bounded native
finite-difference Hessian low-mode probe was also started in tmux session
`model30_ts_mode_probe_attempt2_20261113`, with at most 240 analytic-gradient
evaluations, under `irc_0K/mode_probe_attempt2/`. Its first launch failed
before an energy call due to a relative topology path; that evidence is
retained in `irc_0K/mode_probe_attempt1/`. The corrected attempt uses the same
TS, Hamiltonian and 12 A shell. Until the two lowest modes converge and both
unrestrained endpoints are classified, do not call the candidate a validated
first-order TS or the 46.36 kcal/mol potential-energy difference an activation
barrier.

For context, inserting 46.36 kcal/mol *as if it were* an activation free
energy into the Eyring equation at the 75 C assay temperature (348.15 K)
gives approximately 5.7e-17 s^-1 for a single elementary step. A 20--22
kcal/mol free-energy barrier at that temperature instead gives about
2.0--0.11 s^-1. These rates are conditional illustrations, not rates inferred
from the crude-lysate polymer assay: the computed quantity is a QM/MM
potential-energy difference, the stationary candidate has not passed
Hessian-index validation, and the assay does not specify an elementary
single-turnover rate. Alternative protonation and conformational pathways can
change the barrier itself, not just add a small thermal correction.

**Endpoint result:** both 0 K unrestrained minimizations converged. The
mode-minus branch returned to E-S-like chemistry at -718.207958911553 Eh
(OG1--HG1 0.981 A, OG1--C 3.101 A, scissile C--N 1.360 A). The mode-plus
branch reached acylenzyme-like chemistry at -718.201496183098 Eh (OG1--C
1.329 A, HG1--leaving N 1.018 A, scissile C--N 3.053 A). The 0 K candidate
TS is 46.3912 kcal/mol above the newly minimized E-S endpoint and 42.3357
kcal/mol above the acylenzyme endpoint; the acylenzyme endpoint is 4.0554
kcal/mol above E-S. This establishes *local downhill endpoint connectivity*
for the direct route, not a minimum-energy path or activation free energy.
The direct four-membered proton-transfer/acylation route remains highly
unfavorable on this DFTB3/MM potential surface. Consider alternative
protonation and proton-transfer routes before interpreting it as the enzyme's
productive mechanism.

The first bounded 0 K Hessian spectrum probe completed 120 native-gradient
evaluations but exhausted its 60-HVP cap before ARPACK converged the two
lowest active-shell modes. Its `mode_probe_attempt2/` and launcher records are
preserved. Thus **first-order saddle character is still unverified**. An
improved mode-seeded, bounded retry is being prepared; its result must be
reported separately from the endpoint connectivity test.

The later mode-seeded 0.01 A sqrt(amu) attempt also exhausted its 200-HVP
limit after 400 native gradients without two converged eigenpairs; its
`mode_probe_attempt3/spectrum.json` explicitly records an incomplete result.
DL-FIND exports the analytic force in `path_force.xyz` with only seven decimal
places. A bounded 0.1 A sqrt(amu) displacement-sensitivity retry is active in
tmux `model30_ts_mode_probe_attempt4_20261113`. The existing dimer-mode energy
curvatures at 0.01 and 0.1 A agree to about 1%, supporting this larger
displacement as a numerical diagnostic. It remains unproven that the full
4,165-atom active-shell Hessian has exactly one negative eigenvalue.

A compact 99-frame endpoint-minimization VMD trace is at
`irc_0K/vmd_progress/latest/view_downhill_minimizations.vmd`, with 19,976
display atoms and three retained local waters. The workstation-to-laptop
pull script includes it and the verified relative-energy/distance plot
`irc_0K/plots/endpoint_minimization_progress.png`. The trace is stitched
E-S-to-TS-to-acylenzyme from endpoint optimizations; it is not a formal IRC.

A **fresh regular 0 K native NEB** is running in tmux session
`model30_direct_0k_neb_20261113` under `irc_0K/new_neb/`. It uses the
converged E-S/acyl endpoints, exact TS center seed, identical 104,922-atom
topology and periodic box, and *exactly identical* coordinates outside the
4,165-atom active shell. Its 13 trajectory-derived source structures have
maximum adjacent changes of 0.3064 A (OG1--C), 0.3523 A (C--N), and 0.3052
A (HG1--leaving N), with no new nonbonded heavy contacts below 2 A in seeds
or segment midpoints. DL-FIND interpolates to 16 regular images plus its
reserved climbing slot (`neb_nimage=17`); climbing is disabled initially
(`neb_climb_test=0.0`). The previous 499-cycle band remains archived as
unconverged and must not be merged with this path. Native `nebpath.xyz` is
the coordinate authority. Apply force/energy and geometric acceptance before
calling any profile a converged minimum-energy path; later CI activation
requires a stable regular band and a separate documented continuation.

**2026-11-13 Hessian continuation queue:** The 0.1 A sqrt(amu) low-mode
probe (attempt4) was still advancing through native gradients without an
electronic-convergence failure. Its ARPACK Krylov state cannot be resumed
mid-process, so it remains running. A conditional tmux queue in
`irc_0K/mode_probe_retry_queue/` checks its terminal `status.json` and
`spectrum.json`. Only if attempt4 exits normally but is incomplete will it
launch attempt5 with the *same* stationary coordinates, topology, active
shell, electronic settings, 0.1 A sqrt(amu) difference step, and tolerance;
the cap rises from 200 to 600 Hessian-vector products (400 to 1,200 native
gradient evaluations), and ARPACK maximum iterations from 80 to 240.
The queue skips attempt5 if attempt4 converges and stops on an attempt4
process failure. The exact script is `irc_0K/queue_mode_probe_attempt5.sh`;
all attempt outputs remain separate. Neither attempt establishes first-order
saddle character until two lowest modes and their residuals are validated.

### Update — direct NEB outcome and trace-seeded higher-resolution attempt, 2026-11-14

The first regular 0 K direct NEB in `irc_0K/new_neb/` ended unconverged after
499 accepted L-BFGS steps and 7,000 energy evaluations (the native energy
evaluation cap). Image 7 was an unstable high-energy outlier; the final report
still had residual tangential/perpendicular forces above `tol=1e-4` around
images 5–8. Preserve that band as a failed attempt; do not use its energy
maximum as a barrier or merge its coordinates into the new seed.

A fresh trace-seeded run is active in tmux
`model30_trace_seed25_step005_20261113`, under
`irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/`. It reuses the
13 saved trajectory-derived restart structures from the validated 295-frame
mode-minus and 517-frame mode-plus endpoint traces, the same 104,922-atom
topology, 98-atom QM region at charge -2, 4,165-atom active shell, fixed
endpoints, and 0 K DFTB3/MM Hamiltonian. Only `neb_nimage` changed 17→25,
`maxcycle` 1000→4000 (28,000 native energy evaluations), `maxstep` 0.1→0.05
bohr, and the copied restart paths. Climbing remains disabled.

Preflight verified all mapped atom identities, the same Thr HG1 proton
transferring from OG1 to substrate leaving N 19929–HN 19937, retained Thr
N-terminal H1/H2 and substrate N-HN distances near 1.02 A, exact source
restart/trace hashes, and no new heavy-atom contacts below 2.0 A in the 13
source images or their 12 linear midpoints. The first native path frame also
preserves the retained proton bonds. At that point, the run was active and
post-start checks were pending; see the completion update below.

### Update — trace-seeded direct NEB completed and audited, 2026-11-14

The native 25-slot run exited 0 and DL-FIND reported all regular-band
convergence criteria met at cycle 477. The final native report contains 24
regular images; its maximum is **35.80895 kcal/mol at image 8 relative to
image 1**. The climbing image remained disabled. Native coordinates are
`calculation/nebpath.xyz`; do not use the unusable numbered native NetCDF
restart outputs as path coordinates.

The full-system reconstruction was checked against the 4,165-atom active mask
and fixed 104,922-atom topology. Atom-name/residue mapping and all 24 images
passed endpoint identity and adjacent reactive-coordinate continuity gates.
Image 1 has Thr233 OG1–HG1 0.981 A, OG1–scissile C 3.101 A, scissile C–N
1.360 A, and HG1–leaving N 2.446 A. Image 24 has Thr OG1–HG1 2.430 A,
OG1–scissile C 1.329 A, scissile C–N 3.053 A, and HG1–leaving N 1.018 A.
The substrate carbonyl C–O, leaving N–HN, and Thr N-terminal H1/H2 bonds stay
intact throughout the band. Maximum adjacent changes are 0.388 A (OG1–C),
0.290 A (C–N), and 0.515 A (HG1–N). This validates the intended proton-
transfer/acyl-enzyme endpoint identity and continuity; it does not complete
substrate hydrolysis to measured PA6 L2 product. No separate full-system
nonbonded clash sweep was performed on every optimized image.

The portable 24-frame VMD bundle, energy profile, distance table/plot, and
scientific acceptance record are at
`irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/vmd_final_24img/final_24img/`.
The pull script now includes it. The profile maximum remains a regular-band
potential-energy maximum, not a validated PES barrier or PMF; a climbing-image
TS refinement and first-order saddle validation are still separate work.

#### Qualified comparison with earlier candidates

- The dimer stationary-point candidate is **46.3625 kcal/mol above its
  E–S reference** in the documented one-cycle 0 K evaluation. It is a single
  proton-transfer-dominated stationary-point candidate (OG1–C 2.085 A,
  C–N 1.538 A, HG1–leaving N 1.218 A), obtained with a 300 K electronic-
  smearing dimer search. It is not the maximum of a path profile, and its
  first-order saddle character was not established. Comparing its point energy
  directly with the 35.81 regular-band maximum is not an apples-to-apples
  barrier comparison.
- The earlier 16-image direct NEB used the same 0 K DFTB3/MM Hamiltonian,
  topology, active shell, QM charge, and endpoint pair, but stopped at its
  7,000-evaluation cap after 499 L-BFGS steps. Its report had an unstable
  image-7 outlier at 70.56 kcal/mol and residual forces above tolerance around
  images 5–8. The completed run uses 24 regular images, a denser trace-derived
  seed, a smaller maximum step, and a larger optimization budget; the optimizer
  converged and the reactive distances vary continuously. The lower 35.81
  value reflects a different, better-relaxed discretized band and structure,
  not evidence by itself for a lower physical activation barrier.

The durable run exporter is `export_final_24img.py`. Its identity, continuity,
endpoint, and profile checks are summarized in the bundle's
`scientific_acceptance.json`.

### Update — corrected-shell ASH274 E–S accepted, 2026-11-17

The matched ASP272/ASH274 seed has a force-converged, accepted Tier A
DFTB3/MM endpoint at
`qmmm/pa6_flipped_model30/07_ash274_comparison/06_qmmm_optimization/replica_1_tier_A_seed/corrected_shell/`.
The endpoint retains ASH274 HD2 and neutral Thr233 N. Its d1/d2/d3 are
3.07665/3.06706/2.94462 Å, attack angle 116.807°, and signed phi +97.925°.
The corrected movable shell has 4,650 atoms, including water 1576; its oxygen
shifted 0.1465 Å. This replaces the inherited shell that froze 111 residues
within 12 Å. The earlier accepted run was an incomplete relaxation and is not
the accepted endpoint.

The final, audited water geometry is recorded in
`corrected_shell/water_1576_full_hbond_audit_final.json` (the earlier
`water_1576_full_hbond_audit.json` is retained as interim provenance).
Thr233 OG1–HG1···water O is 167.9° with H···O 1.811 Å; water
O–H2···Asp272 OD1 is 173.4° with H···OD1 1.533 Å. There is no corresponding
bridge to Thr233 N or substrate leaving N. This local endpoint geometry does
not demonstrate water occupancy or reactivity.

A single-water-in-QM reoptimization (102 atoms, charge −1) in sibling
`water_qm_reoptimization/` converged and passed the recorded Tier A and
proton-inventory gates. Its accepted `optimized.ncrst` is the source of the
bounded path below. The water geometry alone does not establish reactivity.

### Current — bounded water-to-Asp272 TI1-like path, 2026-11-17

Tmux session `water_to_asp272_ti1_like_20261117` launched a 25-window
DFTB3/MM coordinate-drive diagnostic in
`qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/`.
Revision 4 of `protocol.json` passed no-launch hash, atom-map, restraint and
target-step preflight. It drives Thr929 HG1 through WAT1576 to Asp272 OD1
while shortening OG1–substrate C; the 12 Å movable shell remains unchanged.
The runner tests force convergence, individual proton identities, bond
geometry and continuity after every window. Each accepted window updates
`vmd_progress/view_path.vmd`, `energy.svg`, and `distances.svg`, which the
laptop pull script includes. An accepted endpoint triggers a separate
unrestrained release. Stop status and scientific outcome are **pending**.
This route targets Asp272-COOH with neutral Thr N, so any product is an
alternative TI1-like state, not canonical TI1 with Thr N as NH3+. Treat the
biased profile as a PES diagnostic, not a PMF or validated barrier.

**Window-5 checkpoint and continuation (2026-11-17):** Windows 1–5
force-converged. Window 5 retained all declared proton identities, with
OG1–carbonyl C = 2.843 Å and physical energy 3.755 kcal/mol above the
unbiased E–S endpoint. The first runner stopped because an unrelated,
approximately 6 Å ASH274 HD2–water O distance changed by 0.445 Å and was
incorrectly included in a 0.35 Å reaction-continuity gate. Protocol revision
5 limits that gate to seven reactive distances; all measured distances remain
reported. The original window-5 acceptance is retained as
`window_005/acceptance_revision_1.json`. Revision-5 preflight validated all
five existing windows, and the same tmux session name was relaunched with
`--resume` from window 6. The completed outcome is recorded below.

### Current — water-relay scan complete; next decision, 2026-11-17

All 25 windows finished with force convergence. No job remains active in this
path. Final status is `endpoint_negative_or_incomplete`; the gated unbiased
release did not run. Both relay transfers occurred, yielding water-bound HG1
and Asp272-bound H2 while ASH274 retained HD2 and Thr N stayed neutral.
OG1–C ended at 1.852 Å (target 1.650 Å). Physical relative energy increased
monotonically to +54.0235 kcal/mol; bias-inclusive energy is +58.1547.
This scan neither found TI1 nor demonstrated its absence. It tests a
both-Asp-protonated branch, not the Thr-NH3+/oxyanion intermediate.

Read the **Controlling decision — after completed water-relay scan** at the
top of `docs/NYL12_MECHANISM_ROADMAP.md` before launching further calculations.
Next proposed work is a bounded 5 × 5 relaxed map of direct Thr OG1→Thr N
proton transfer versus OG1–carbonyl C attack on the same accepted 102-QM-atom
E–S seed, with the two targets varied independently. Follow informative
points with sparse matched-embedding method checks and candidate releases.
Do not claim a stable TI1 is mechanistically mandatory; a coupled path to
acylenzyme remains an alternative. This plan is documented but not launched.
Updated manuscript interpretation is in `manuscript/RESULTS_DRAFT.md`.
