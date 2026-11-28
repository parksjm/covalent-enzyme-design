# Nyl12: active plan through 30 November 2026

**Confidence scoring qualified, 28 November:** ipSAE2 0.1.11 is installed in
an isolated, locked uv environment at the reviewed Git revision. WT and
design 06 passed native-token mapping, all-four-ligand reporting and 896
independent directional numerical checks. Supplementary scoring uses a
separate postprocessing job; active Boltz/Amber and QM/MM/US inputs remain
unchanged. Scores describe prediction confidence, with no PA6-specific cutoff
or catalytic-activity claim. See [qualification and usage](IPSAE2_QUALIFICATION.md).

**Relay sampling update, 28 November:** the common analyzed main horizon is
70 ps across 57 windows. The new pooled fit includes all 43 eligible add-on
histories: 33 corrective activation windows at 10 ps, four targeted checks at
25 ps and six additional initialization checks at 30 ps. With the fixed 5 ps
discard once/history, the pool contains 82,000 frames and has a projected
maximum of 26.28 kcal mol⁻¹ at s=0.6875. The original 57-only estimate is 26.12;
its early/late maxima are 25.53/26.78. Neither estimate is accepted as a
converged activation barrier. Pooled uncertainty remains unresolved because
the new grid histories retain only one 5 ps block each; gray intervals in the
new comparison plot apply only to the original 57. All 33 grid and four target
workers continue on Violet; main continuation to 100 ps runs on Azalea. The
single 75 ps-complete main window waits normally for its 56 peers at the scheduler
barrier, without an error. Geranium's six 30 ps checks are finished.

The refreshed main VMD bundle covers70 ps/window with 1,596 frames and the
required camera/display policy. It explicitly includes the first 5 ps for
visual inspection, whereas PMF statistics discard them. Supplemental movies
remain separately identified. The unchanged movie API is
`vmd_progress/current/view_path.vmd`; the updated laptop pull now follows the
extended-current PMF monitor, not the stale60ps monitor. The readable dated
plot is `plots/relay_pooled_70ps_20261128.png` with matching PDF/CSV/JSON.

**Matched design comparison completed, 28 November:** all six Boltz predictions
finished. WT seed 20261126 relaxed to four positive-phi Tier A sites. Design
03 seed 20261126 relaxed to C/B/B/B (three Tier B), and its second seed has
no Tier A/B sites. WT seed 20261127 and both design-06 predictions have severe
interchain overlaps and were rejected before Amber. All three launched
minimizations force-converged. Neither design has demonstrated an improvement
over the accepted WT model. Native ipSAE2 is available for all six, with no
catalytic threshold; do not assign four unique occupied sites to an invalid
assembly. [All-site geometry and confidence report](../design/nyl12_relay_ts_ensemble_20261126/broader_pilot/original_backbone_matched_20261127/comparison/update_20261128/README.md).

**Latest PI clarification, 27 November:** The current design workflow excludes RFD3. Its candidate pool comprises the
12 fixed-backbone LigandMPNN designs from two chemically coherent MD-derived
TS-like parents, with the 101-position protection mask and identical sequences
across all four copies. Thermal conformational asymmetry is allowed; exact D2
coordinates are not imposed on MD snapshots. Parent acceptance independently
requires reaction mapping, proton inventory, ligand graph and catalytic geometry;
multiple chain conformations alone do not establish TS compatibility. Historical
RFD3 remodels and their 12 sequences remain diagnostic artifacts outside this
current pool; no replacement generation or RFD3 qualification is required.

Fixed-backbone LigandMPNN design followed by independent template-free,
paired-MSA Boltz prediction and matched Amber relaxation can assess
sequence-driven repacking and assembly; it does not promise improved activity.
Only clear evidence of a backbone limitation would motivate a future controlled
RFD3 comparison evaluating coupled perturbations of A/E, B/F,
C/G and D/H, with independent physical-movement, clash, assembly, fixed-context
and chemical checks. It is not implemented or qualified for production, and
exact coordinate symmetry is not a universal acceptance requirement.

QM/MM PMFs are **not prerequisites for experimental testing**. They are an
optional matched comparison before testing or a mechanistic explanation
afterward. PI nomination can proceed from the available structural/dynamic
evidence, with appropriate experimental controls and bounded claims.

**Sign correction, confirmed by the PI on 27 November:** the desired Boltz
orientation is **positive phi** under the unchanged project definition. The
earlier assistant claim that the target was negative was incorrect and is
withdrawn. Preserve all measured signs and historical cohort results. New
E–S selection must explicitly use `target_phi_sign=+1`, calibrated to
`approved_positive_pose_reference_20261127.json`; no phi restraint is added to
the ongoing reaction simulations.

Updated 27 November. This is the short operational plan. The mechanism roadmap
and manuscript retain detailed evidence and dated history; live JSON receipts
establish current execution. Do not reload the full project history for a
routine update. The objective is a defensible acylation pathway and a bounded
TS-ensemble redesign pilot, not a barrier adjusted to a desired value.

**Superseded next-round proposal, PI discussion on 27 November:**
The later explicit RFD3 exclusion above governs the current workflow; the
following proposal is retained as dated history, not an execution plan.
 keep every current
QM/MM, US and Boltz job unchanged. Prefer a small comparison of original
TS-like backbones versus gentle whole-assembly partial diffusion, allowing
active-site loops to move while retaining catalytic, maturation and assembly
anchors. Keep WT chain lengths and broad, four-copy-tied LigandMPNN sequence
freedom; no loop indels or de novo fold campaign. Qualify the global mobile
mask and resulting geometry before expansion. Retain TS-derived structural
input for RFD3; external crystal templates for prediction are a separate input.
The preferred next prediction screen is MSA-based and free of structural
templates, provided a matched WT control recovers the native pocket. Earlier
MSA-only WT failures remain a calibration limit, not design-inactivity evidence.
Test precursor-derived alpha/beta MSA pairing as one specific improvement;
retain eight mature protein chains in the structure prediction. Existing
independent alpha/beta MSAs cannot be paired by synthetic row-label matches.
The working local ColabFold installation and databases were recovered from
enzyme-substrate-modeling. A local 328-residue precursor search succeeded and
provided 1,387 valid homolog pairs; no installation, download or public-server
request was needed. Template-free WT controls passed native paired-MSA and
dummy-template feature checks and completed inference with exit 0 on the
authorized DGX GPUs 2/3 after the earlier queues finished. Both canonical
exports pass, and matched Amber follow-through is active in persistent tmux.
Before relaxation the seeds had 4/4 (one A, three B) and 0/4 positive-phi A/B
sites. The first is now independently accepted with four Tier A sites; the
second failed assembly/minimization, as detailed below. See [the verified local MSA
recipe](LOCAL_MSA_WORKFLOW.md). Use the current corrected six-model run as a retained diagnostic,
not an implicit decision to keep crystal templates in every future screen.


**Paired-WT Amber outcome, 27 November:** seed 20261126 converged at RMS
0.047399 kcal mol⁻¹ Å⁻¹ and retained four positive-phi sites, improving from
one Tier A/three Tier B before relaxation to four Tier A afterward. This is
an independently qualified WT baseline for the corrected sequence-only
24-sample pilot. Seed20261127 is a failed assembly control: the native model
has 576 interchain heavy-atom pairs closer than 1 Å (minimum 0.1225 Å). Amber
stopped with repeated line-minimization failure and RMS 53.676; its distorted
endpoint must not be continued or treated as relaxed. Preserve the failure
and do not claim two-seed reproducibility or use a mixed accepted bundle.
Future prediction screens must check severe interchain overlaps before
spending minimization time. Existing QM/MM and US jobs remain unchanged.


**Relay checkpoint, evening 27 November:** the current pooled fit is
`checkpoint_011_main60ps_supplement10`: 57 main histories at 60 ps, four earlier
checks and six Geranium checks, totaling 67 histories and 64,100 retained
frames. Its provisional projected maximum is 26.054 kcal/mol at s = 0.6875.
This is a point estimate; the main 57-only uncertainty band is not pooled
uncertainty, and convergence is not established. The six Geranium checks
are now included; none of the 33 corrective-grid histories is eligible yet.
First grid contributions are estimated for Saturday 28 November morning
(approximately 08:00–10:00 EST), with the planned full 100-history assessment
targeted Sunday night/Monday 30 November, subject to throughput and diagnostics.

The main movie is current at 60 ps/window (57 histories/6,840 frames). Its
automation does not include the new grid or Geranium histories; the older
supplemental export originally contained only four checks. Updated separate
post-discard old4/Geranium6 exports now contain 800/600 frames. Persistent
verification03 jobs and `relay_post5ps_index_collect03_20261127` own all-frame
source/time/coordinate/box validation, atomic publication, index and pull
updates. Original exports remain immutable; derived DCDs correct only cell
record ordering for VMD compatibility. These are two supplemental viewers,
not automatic movie coverage of all 100 histories. Grid export follows its
accepted-data threshold.

Azalea's original 57 histories completed 60 ps and have now resumed toward 100 ps.
The active Violet grid worker
has no safe block-drain/resume hook, so an ad hoc migration is not performed
while preserving all current sampling blocks. Its parent also owns reserved
paused jobs whose lifecycle must not be bypassed. Host-specific engine pins
differ. The prepare-only assessment is retained under
`analysis/activation_support_preparation_20261127/` in the relay campaign.
The resumed main 57 now utilize Azalea through the existing accepted-prefix
worker. Future grid workers should expose a tested drain-after-block and
restart interface.

## Priority 1 — finish the current relay sampling assessment

- Keep the existing central-amide relay Hamiltonian, topology, proton inventory,
  300 K protocol and movable region unchanged. No new NEB or mechanism branch.
- Main 57 histories are running 60→100 ps on Azalea (57 single-thread workers;
  56 logical CPUs, slight oversubscription avoids starving one window). Violet
  runs the 33 corrective
  grid plus four earlier checks; Geranium runs six independent-start checks.
  Planned pool: 100 histories, main 100 ps and supplements 30 ps each. The
  original main 60-ps analysis remains an immutable earlier checkpoint.
- The automatic eligible-data updater incorporates each accepted supplement
  after 10 ps production, discarding its first 5 ps once. Latest checked fit:
  main 60 ps plus four earlier checks and six Geranium checks, 67 histories,
  64,100 retained frames, projected maximum 26.054 kcal/mol at s = 0.6875. It remains unconverged and point-only;
  the displayed original-main confidence band is not pooled uncertainty.
- Geranium contributions are now included. Expected first grid contributions:
  28 November morning. Next common 65-ps main milestone is about four hours
  after the 02:52 UTC launch (roughly 02:00 EST on 28 November) ; 100 ps is projected
  in ~32 hours (roughly 06:00 EST on 29 November). The grid remains the slower cohort.
  Comprehensive corrected analysis is
  expected Sunday night/30 November, subject to measured throughput.
- Final assessment must include local/global overlap, physical-coordinate and
  projection sensitivity, weight support, proton inventories, temporal halves,
  independent starts and 5/10 ps block sensitivity. An operational precision
  goal is <1 kcal/mol late-checkpoint/initialization disagreement, accompanied
  by all other diagnostics; no single number proves convergence.
- If gaps/drift remain on 30 November, release the supported profile and an
  explicit unresolved-region/continuation list. Do not smooth away the problem,
  change discards to lower the barrier, or label a provisional profile converged.

Live source: `qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/`
`full_path_300K/analysis/current_eligible_progress_20261127/current/`.
Final source: the same campaign's `analysis/final_candidate_20261127/`
`monitor_state.json` and `watcher_job_attempt03/`. The tmux session is named
`relay_final100_geranium_20261127`; it is not a directory name.


## Priority 2 — qualify the audited 101-position design cohort

- **Current protections:** canonical 113, 156, 186, 221–235 and 272–276 plus M1 and
  the recurring inter-repeat contact mask. The union is 101 protected / 227
  variable precursor positions: 24 core policy plus 77 additional interface
  positions. K156 and L235 are fixed. The complete list/reasons/figure are in
  [the SI inventory](figures/nyl12_sequence_protection_20261127/SI_caption_and_inventory.md).
- **Negoro audit:** NylC K189A has cleavage-normalized hydrolysis activity
  ~0.6% of parent, supporting precautionary K156 retention. Exact Lys chemistry
  and transfer to Nyl12 remain unproven. Other catalytic/maturation counterparts
  in the 2012/2018/2023 papers are covered by the current mask; stability/
  processability modifiers are not automatically called essential. See
  [the gap-aware audit](nylc_critical_residue_audit_20261127.md).
- **Calibration:** paired-MSA template-free WT seed 20261126 converged at
  RMS 0.047399 and has four positive-phi Tier A sites. Seed20261127 failed
  assembly/minimization and remains excluded; no two-seed robustness claim.
- **Current execution:** `broader_pilot/k156_protected_pilot_20261127/`, session
  `nyl12_prot101_pilot_20261127`, owns a replacement 24-raw-sample run: six on
  each of two original and two historically accepted RFD3 parents, unchanged temperature/
  seeds/engines. Fresh 101-position flags and ONE frozen package serve both
  generation and scoring. All four TS and E–S native parser checks and GPU
  prerequisite passed. Generation and all cross-scores are complete. Review
  the existing 12 original-parent designs using their completed scores; the
  12 RFD3-parent designs remain outside the current baseline comparison.
  No regeneration or new scoring calculation is needed.
- The preceding 24 `relay_corr01_*` sequences all changed K156 and are
  ineligible. The score-only diagnostic was stopped after 118 evaluations;
  data/failures are preserved. Do not silently revert old sequences or
  reuse their identities for replacement designs. Earlier broad designs
  also remain diagnostics under restored 275–276 protections.
- **Next steps and budgets:** review the existing 12 original-backbone designs;
  no replacement generation or RFD3 qualification. Then small matched template-free
  Boltz cohort → canonical/clash checks → matched Amber relaxation → all 48
  positive-phi/Tier A/B measurements. Native-Thr similarity is diagnostic.
  At most two initial designs enter the matched two-seed WT/design prediction
  comparison; expensive MD is limited to two or three finalists and QM/MM/US
  further limited. No prediction, MD or QM/MM design campaign is automatically
  released by sequence likelihood. See the [full scripted protocol](NYL12_TS_ENSEMBLE_REDESIGN_PROTOCOL.md).
- Amber is the current matched relaxation backend. OpenMM is an available
  alternative that needs matched control qualification, not an extra mandatory
  step. No new loop indels, capped-force-field development or whole-assembly
  diffusion launch is included here.

## Register and assay rules

- PA6 L4 is a manageable local proxy for a longer, flexible polymer. L4 itself
  is observed as a product; its subsequent L4→L2+L2 turnover remains a hypothesis.
- Original contact-guided models and the current relay use C54–N59, between
  units 2 and 3: **L2+L2 for L4**. The newly relaxed WT controls use C55–N57,
  between units 3 and 4: **L3+L1 for L4**. Distinguish bond-unit indices from
  fragment lengths in every table. Record all candidate registers explicitly.
- Preserve the agreed general-cleavage objective: any internal amide can pass
  geometry screening. For direct comparison to an existing TS/path, retain
  its coherent register/atom map or explicitly rebuild the comparison. Do not
  reuse central-amide QM masks on a terminal-amide model or pool their PMFs.
- The manuscript reports PA6 L2/L3/L4, with L3 prominent. Its MRM list does not
  include a PA6 L1 target. Therefore unreported L1 is not proof that L1 was
  absent, and a finite-L4 split does not establish long-polymer selectivity.
- Lysate data have uncertain active-enzyme concentration; K59E duplicates are
  discordant. These measurements support provisional activity comparisons,
  not precise kcat/barrier calibration. Purified WT plus a few finalists and
  complementary amine-release assays remain the experimental follow-up.

## Loop and capped-model work: preserve context, limit new scope

- The qualified RFD3 partial-diffusion CIFs contain both loop regions, but
  their current fixed/mobile masks fix alpha 95–103 and beta 289–300 in all
  four copies. The 24 mobile residues are elsewhere in the pocket. Current
  remodeled parents therefore do not test remodeling of these loops. A
  subsequent fixed-length 95–103 experiment needs a revised native mobile
  mask, with its starting loop coordinates retained; no CIF deletion is
  required for this purpose. Keep predictor-template coordinate omissions
  and true sequence-indel/inpainting specifications distinct. No new RFD3
  run is launched while the corrected six-model screen is in progress.
- First true loop-length target: **WT residues95–103**; next **289–300**.
  In the observed A/E pocket these are neighboring D-alpha and F-beta loops,
  respectively. Keep tetramer context and exact four-copy sequence edits.
- Historical loop/tail-omitted CIFs are inventoried in
  `docs/evidence/commercial_polyamides_20261127/loop_mask_template_inventory.json`.
  Their filenames are not exact mask definitions; some omit additional loops
  or contain PEG/water. Reuse the strategy through explicit audited masks.
- Defer new length-changing RFD3 runs until the present design screen is read.
  If essential, qualify one 95–103 loop insertion/deletion arm; do not start
  289–300 simultaneously. Regenerate WT-reference/design/Amber mappings,
  protections, cleavage/terminal-charge maps and topology indices before use.
- Preserve the existing capped/register model collections. Legacy capped GAFF
  payloads exist; their exact historical charge method is not yet certified.
  They are not equivalent to the high-quality modular RESP/GAFF2 release.
  No new capped RESP parameterization or capped MD/PMF comparison before
  30 November. Do not mix charge/FF releases in quantitative comparisons.
- Defer deacylation expansion, broad new protonation/pKa campaigns, LLM
  benchmarks and agentic-system expansion. Preserve all existing results.

## Deliverables and usage after the deadline

By 30 November: supported acylation PMF plus convergence audit; chemically
mapped, clustered TS-like design ensemble with validation status explicit;
bounded design-model results and a justified shortlist or documented no-pass
outcome; updated manuscript/SI and presentation. A verified first-order saddle
or equilibrium TS ensemble must not be claimed without its own validation.

From 1 December, target approximately tenfold less AI usage through scripted
execution and event-driven review. This is an operating target, not a billing
guarantee: keep deterministic tmux controllers/analysis running; use concise
JSON/status receipts; consult this brief and changed evidence first. Reserve
LLM review for milestones, failures and PI decisions. Do not run an autonomous
LLM reasoning loop or repeat completed software/resource audits. Delegate only
bounded tasks when total cost/quality improves, without full-history forks.
Keep the append-only action journal and immutable scientific provenance.

Evidence: [manuscript extraction](evidence/commercial_polyamides_20261127/),
[design policy](design_ranking_and_iteration_20261022.md),
[mechanism evidence/roadmap](NYL12_MECHANISM_ROADMAP.md),
[SI methods and diagnostics](../manuscript/SUPPORTING_INFORMATION_DRAFT.md).
