## PRIMARY OBJECTIVE — qualify relay sampling for TS-based redesign, 27 November 2026

**Latest design/experimental gate, 27 November:** RFD3 is outside the current critical path. Review the existing 12 original
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

PMFs are optional before experiments, not a nomination prerequisite. With PI
review, test a bounded WT/design cohort using available structural and dynamic
evidence; use matched QM/MM PMFs before testing when feasible or afterward
to investigate measured improvements. Current QM/MM/US jobs continue.

**Current PI decisions, 27 November:** K156 and L235 are protected. The current
mask has 101 fixed / 227 variable precursor positions. All 24 `relay_corr01_*`
designs change K156 and are ineligible; diagnostic outputs are preserved and
only that score-only job was stopped. No QM/MM or US job was interrupted.
The Negoro mutation/mapping audit precedes replacement design launch.
Azalea's relay campaign is extended from 60 to 100 ps/window with 57 single-
thread workers; Violet's corrective grid continues. Use the authoritative
[TS-ensemble redesign protocol](NYL12_TS_ENSEMBLE_REDESIGN_PROTOCOL.md) and
[SI mask inventory](figures/nyl12_sequence_protection_20261127/SI_caption_and_inventory.md).
Older K156-variable/90 ps proposals below are superseded, not active policy.

**Protection revision after PI review:** K156 is variable for the corrected
24-sample pilot. Its former protection was based on NylC Lys189 analogy and
local electrostatic context, not established Nyl12 essentiality. L235 is mature
beta residue 3 and remains a conservative N-terminal-context protection for
this pilot; its individual catalytic importance has not been demonstrated.
The active minimum is 113, 186, 221–235 and 272–276 plus existing structural anchors.
The preceding 101-position mask remains an immutable historical artifact;
new preparation must generate and hash its own corrected mask.

**Precursor-paired template-free WT controls, 27 November:** the existing
Violet ColabFold/MMseqs environment and indexed databases are operational.
A local search of the exact 328-residue WT precursor supplied 1,387 retained
same-precursor homolog pairs, split at query-aligned 232|233. Both matched WT
predictions completed natively with exit 0 on DGX GPUs 2/3; actual dataset
checks accepted coherent eight-chain pairing and stock dummy-template features.
No installation or public-server request was needed. Both canonical exports
passed full-precision CIF and true mapping-roundtrip checks; matched Amber
follow-through is active in persistent tmux. Before relaxation, seed 20261126
has one Tier A and three Tier B positive-phi sites; seed 20261127 has none.
These are geometric observations, not predicted activity. The source, jobs and acceptance receipts are
under `broader_pilot/wt_paired_precursor_controls_20261127/`. See
[the recovered local MSA procedure](LOCAL_MSA_WORKFLOW.md).

**Previous restored hard sequence protections (superseded for K156):**
preserve canonical 113, 156, 186,
221–235 and 272–276 in all four copies, plus existing structural/interface
protections. This includes the entire PI-specified alpha tail 221–232.
Recovered historical 275–276 protections are conservative policy, not proof
of direct catalytic participation. The hardened prospective mask contains 101
fixed / 227 mutable positions. All 24 existing broad designs mutate 275;
relay_02 has V275M/L276T and relay_23 has V275I. Retain their current jobs as
diagnostics, but do not advance them or silently repair their sequences.
New designs require new identities and revalidation under the restored mask.
Sequence protection does not by itself freeze RFD3 backbone coordinates.


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
supplemental export contains only four checks. An updated supplemental
export must be qualified and published separately; do not claim it already
updates automatically for all 100 histories.

Azalea's original 57 histories have completed. The active Violet grid worker
has no safe block-drain/resume hook, so an ad hoc migration is not performed
while preserving all current sampling blocks. Its parent also owns reserved
paused jobs whose lifecycle must not be bypassed. Host-specific engine pins
differ. The prepare-only assessment is retained under
`analysis/activation_support_preparation_20261127/` in the relay campaign.
Use idle Azalea for nonintrusive upcoming analysis when feasible; future
workers should expose a tested drain-after-block and restart interface.


**Corrected 24-sample LigandMPNN pilot launched, 27 November:**
`broader_pilot/protection_corrected_pilot_20261127/` uses the same two original
and two accepted 5 Å-remodeled TS backbones, six samples each, four-copy
sequence ties and the new 100-fixed/228-mutable mask. K156 is variable; L235
is retained conservatively. Accepted paired-MSA WT seed 20261126 is the
single-seed calibration control; seed 20261127 remains a failed assembly.
New design IDs use `relay_corr01_*` and do not replace the previous 24 designs.
The first preflight stopped on copied Python bytecode hashes; a second
preflight encountered retained effective-manifest paths. Neither generated
sequences. Both failures and their manifests/releases remain preserved.
`manifest_attempt02.json` excludes only generated caches, verifies 397 actual
source/input hashes and disables bytecode writes.

The project launcher owns session `nyl12_protection_corrected_20261127_attempt02`
and `job_attempt02/status.json`. All four TS native parsers and the E–S parser
passed; the GPU prerequisite exited 0 and first native generation began on
Violet GPU UUID `GPU-ec7e1af2-0974-2893-d34b-411eee82bbd6`, capped at two CPU
threads. The bounded run includes deduplication and matched TS/E–S scores;
no new RFD3, Boltz/AF3, MD or QM/MM variant campaign is released by this step.
Existing QM/MM and US continue unchanged. Source/release proofs are retained
with the campaign, including the independently accepted WT artifact.


**Generation complete; score-only recovery, 27 November:** all four native
generation calls exited 0 and produced 24 raw / 24 unique `relay_corr01_*`
sequences. The combined job then exited 1 because its nested historical
scoring adapter still protected K156. This was a sequence-permission mismatch
before candidate inference, not a force-field or generated-sequence failure.
The generated pool and all five successful WT scores remain immutable.
`scoring_recovery_attempt01/` supplies the coherent current policy package
and a score-only continuation using the same engines, weights, templates,
mask, seeds and exact WT decoding orders. It does not regenerate sequences.
The new continuation passed 9 regressions (81% statement coverage) and Ruff;
source/input hashes and separate process/acceptance records own recovery.


**Score-only continuation active:** session
`nyl12_corrected_scores_resume01_20261127` and
`scoring_recovery_attempt01/job_attempt01/status.json` own the remaining
120 matched candidate evaluations. The first two native candidate scores
exited 0 under the coherent K156-variable policy and matched WT-order checks.
Five accepted WT scores are reused. All 24 generated sequences remain
unchanged; score-only validation has 9 passing tests and 81% statement coverage.

**Pose-sign correction, 27 November:** the PI-approved Boltz orientation is
positive phi under the existing mathematical definition. The assistant's
earlier negative-only target was wrong and is withdrawn. Numerical phi values,
legacy tiers and historical negative-cohort analyses remain unchanged. Active
E–S model selection requires explicit `target_phi_sign=+1`, not an inferred
meaning of “flipped”. See the approved reference manifest and figure guide.

**Active execution priorities:** read the concise
[plan through 30 November](NYL12_THROUGH_NOV30.md) first. It supersedes older
prospective task lists below; native JSON receipts establish current execution.
The lower dated sections retain methods, evidence and historical decisions.
The latest checked relay fit uses main 55 ps plus four 15 ps checks, 61
histories/57,800 retained frames, projected maximum 26.028 kcal/mol; it is
still provisional. The new six-model Boltz screen uses WT and two existing
designs with matched seeds and audited loop-coordinate-masked templates.
True loop-length changes and new capped-substrate parameterization are
deferred while acylation sampling and this bounded screen take priority.

**Template alignment repair, 27 November:** all six original relaxed
models completed, including both relay_23 seeds. All four design endpoints
have zero positive-phi Tier A/B sites; original WT endpoints have zero and
two. The initial WT-only discrepancy was in alignment records. Actual native
template tensors are identical across all six predictions: Boltz does not
clip to the record boundaries, and alpha query17–222 receives features,
including all 36 intended masked loop positions95–103 across four copies.
Trimming the WT record leaves tensors unchanged. Thus original conditioning
was shared but failed the intended loop omission in every case.
An isolated native-feature correction is being prepared, with source/query
mapping and actual mask verification before inference. A rigorous comparison
requires six replacements (two WT and the same four designs), not corrected
WT against original designs. Preserve sequences, seeds, MSAs, weights,
unguided PA6 and relaxation protocol; keep installed Boltz unchanged and
retain original outputs as exploratory evidence. No new sequences or
mechanistic branches are added. Existing QM/MM and relay US continue. Native
receipts establish launch/completion. The original diagnostic is under
`wt_alignment_repair01/`; the corrected workflow is
`template_alignment_repair01/`. All six native dataset feature paths now pass:
1,196 observed protein residues per case, bit-identical native features versus
an independently author-indexed reference, intended masked positions absent
and original processed source hashes unchanged. Independent adapter tests:
14 passed, 100% statement coverage. Prediction/relaxation acceptance remains
separate from this input qualification.

**Current design steering, 27 November:** prioritize Boltz for the next
matched design-model screen; AF3 need not be run for every sequence. New
designs advance only with a positive-phi (PI-approved mapped pose)
Tier A/B site after unrestrained relaxation. Evaluate all 48 candidates with
this combined eligibility before choosing a best site, preserving every
site's data. Positive-phi WT Tier A controls remain valid observations and
satisfy this pose and geometry requirement; a flipped WT pass is
not a prerequisite for testing the new sequences. Keep the fixed-length
24-design pilot intact. A small one-loop L−1/L/L+1 reconstruction arm is
proposed separately, with identical indels in all four copies and fully
regenerated numbering/protection/topology maps. This arm and its added budget
have not been launched. See [the design policy](design_ranking_and_iteration_20261022.md).

**Matched template-control relaxation, completed 27 November:**
`wt_template_amber_unrestrained_20261127` finished both controls and their
canonical four-site audits with exit zero at 15:25:39 EST on Violet GPU 1,
without changing relay-US scheduling. It reuses the validated four-ligand
`relax_model` protocol: ff19SB, modular PA6 RESP/GAFF2 release, TIP3P 12 Å
octahedral solvent, 150 mM salt plus neutralization, four neutral Thr termini
and anionic Asp272/Asp274 at all four sites. One minimization stage permits
10,000 cycles (1,000 initial SD), RMS target 0.05 kcal mol⁻¹ Å⁻¹,
`ntr=0,nmropt=0,ntc=1,ntf=1`; every atom is movable. There is no restrained
solvent stage. Source/code/parameter and engine hashes are frozen under
`wt_protein_template_controls_20261127/amber_relaxation_20261127/`.
Both native minimizations met the matched RMS target: AF3 0.047702 and Boltz
0.045807 kcal mol⁻¹ Å⁻¹. Postprocessing used explicit matching dry topologies,
restarts and conversion maps, recanonicalized both snapshots, measured all 48
combinations and published paired PyMOL bundles with fixed-identity deltas
plus best/nearest assignments. Boltz moved from four rejects to four Tier A
sites (attack angles 100.46–108.26°, both oxyanion contacts 2.84–2.99 Å).
AF3 moved from four Tier C to four rejects: attack angles 73.59–79.68° and
Asn contacts 4.16–4.32 Å. Selected site/ligand/amide identities did not change;
all three amides at all four sites were searched. Boltz retains positive phi
(+84.20 to +88.23°); AF3 retains negative phi. Thus the Boltz Tier A outcome
is a positive-phi WT geometric control and matches the PI-approved
positive-phi pose selection. Its Nα–leaving-N distances nevertheless shorten
to 3.77–3.97 Å versus 5.76–6.39 Å for relaxed AF3; retain this evidence when
assessing substrate directionality. No MD or QM/MM validation is implied.
See SI S10 and per-case `postprocess/before_after_site_summary.json`. Process
completion, force convergence and geometry acceptance remain separate. The
pull script includes each bundle only after its atomic manifest appears.
LigandMPNN cross-scoring had completed naturally with exit zero at 14:34:58 EST;
no existing job was stopped to launch these controls.

**Canonical inspection and current pooled fit:** the updated
`scripts/pull_nyl12_inspection.sh` includes the four original WT controls,
PyMOL loaders, biological and mature-local numbering, and sequential Amber
heavy-atom inputs with explicit conversion maps. Thr233 is E1 in mature-local
exports and residue929 in the corresponding Amber ordering. Four-L4 controls
have 1,328 dry residues; one-L4 TS parents have 1,316. Re-export from matching
topology/restart/conversion metadata after Amber, then restore canonical
identity before RFD3; never reuse a residue offset without its map.
**Site selection:** Amber residue 929 identifies one catalytic Thr in this
ordering; it does not identify the best substrate pose. With four unguided
substrates, evaluate all four physical catalytic sites and all internal amides,
then select by the agreed catalytic geometry criteria rather than proximity
to residue 929 alone. Retain the selected alpha/beta chains, ligand identity,
scissile amide and topology-specific atom/residue map for downstream work.
Build restraints, QM masks and viewer labels from that selected map; never
silently reuse the historical residue929/C54/N59 assignments for another site.
The two completed MSA controls are also published in their own optional
`wt_msa_controls_20261127/inspection_export_20261127` bundle. Both completed
protein-template controls are also available in their separate
`wt_protein_template_controls_20261127/inspection_export_20261127` bundle:
all eight WT predictions can now be inspected without replacing earlier data.

The accepted-data updater `relay_current_eligible_20261127` checks every 180 s
and fits only changed accepted sources. Its first pooled fit includes the
57 main histories at 50 ps and four earlier checks at 10 ps: 61 independent
histories, 59 distinct biases, 51,700 retained frames and a maximum of
25.913 kcal/mol at s = 0.6875. The pooled result is point-only: those supplemental
histories retain only one 5 ps block after the one-time 5 ps discard. The displayed
95% band belongs to the original 57-history fit, not the pooled estimate.
The new 33-grid and six Geranium checks are still ineligible. They will enter
as each history reaches the declared production/data-acceptance threshold.
The 100-history final analysis and longer-block/initialization checks remain
required; neither the current maximum nor narrow conditional intervals prove
convergence.

**Supplemental eligibility ETA, measured at 14:28 EST on 27 November:** all
six Geranium checks completed 5 ps equilibration and are 26–27% through the
first 5 ps production block. Their recent rate is about 0.96 s/step, giving
first/all six 10 ps production eligibility around 19:00–21:00 EST Friday
27 November. All 33 Violet grid windows completed 1 ps equilibration and are
42.5–58.8% through the following 4 ps equilibration; no grid production is yet
eligible. Their 2.47–2.89 s/step rate gives first grid contributions around
07:00–10:00 EST Saturday 28 November, and all 33 around 08:00–12:00. Allow
several hours of rate uncertainty. Each independent history enters after two
consecutive accepted 5 ps production blocks, with the first discarded once;
the 180 s updater then detects changed accepted inputs. The existing four
checks already contribute. Eligibility is not convergence. Exact native stage
progress, rate arithmetic, source identity and 15 ps checkpoint estimates are
in [the ETA audit](../qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/supplemental_eta_audit_20261127/status_eta.json).

**Current decision:** the 33-window activation supplement is necessary but
insufficient by itself. Separate late-activation initialization checks and
TS-region checks are now queued. Use Violet and Azalea for the main relay work; the PI additionally authorized
six CPU checks on Geranium on 27 November, leaving its other user’s GPU work
untouched. Preserve every original trajectory. The detailed construction, diagnostics,
restraints, ESS definitions and corrective protocol are in
[Supporting Information S1–S8](../manuscript/SUPPORTING_INFORMATION_DRAFT.md).

| Task | Verified state | Next checkpoint |
|---|---|---|
| Original 57 relay windows, Azalea | All reached 50 ps; continuing to 60 ps | Automatic common-horizon PMF/movie; main 60 ps expected Friday evening |
| 33-state activation supplement, Violet | Running; 5 ps equilibration + 30 ps production | Common 10 ps production: point/support diagnostic; 15 ps: first eligible conditional block intervals |
| Four earlier reciprocal/product checks, Violet | Native 10 ps complete; continuing in place to 30 ps | Independent lineage statistics; no repeated discard on continuation |
| Six new mixing checks, Geranium | Running in `relay_mixing6_geranium_20261127`; native compatibility passed | 5 ps equilibration + 30 ps production; compare independent initializations |
| Final most-extensive relay analysis | Monitor `relay_final100_geranium_20261127` active; no fit yet | All main 60 ps + all supplements 30 ps; global MBAR, initialization/longer-block checks |
| Broader redesign | Both 5 Å remodeled parents accepted; all 24 raw samples generated and unique; native cross-scoring completed with exit zero | Audit ranking and matched WT relaxation; variant prediction still requires a usable WT baseline |
| WT predictor/Amber controls | Raw controls and matched endpoints canonicalized/published; relaxed Boltz has four Tier A sites, relaxed AF3 all reject | Inspect positive/negative phi classes separately; static geometry still requires dynamics/mechanistic validation |
| Current supplemental MBAR | Four completed 10 ps checks integrated: 61 histories / 59 distinct bias potentials with main 50 ps | Automatic accepted-data updater; new 33-grid and six mixing checks enter when individually eligible |

**Exact redesign gates, updated 27 November:** the exported parent01 5 Å
adapter now passes the unchanged scientific audit. It proves the native
source-centered BF16 signature before restoring only fixed source coordinates,
the typed source graph and eight omitted fixed OXT atoms; mobile coordinates
are not repaired. Its 2 Å alternative retains genuine Arg118 NH1–Ile122 O
(1.648 Å) and NH2–Ile122 CB (1.791 Å) clashes and remains rejected. The second
TS parent03 also passed its independent post-run audit (mobile Cα RMSD
0.314 Å; source sequence, atom inventory and typed ligand graph retained).
The old failed controller remains stopped. The new source-hashed queue has
generated all 24 raw samples: six on each of the two original TS parents and
their two accepted 5 Å remodels. All 24 are unique and preserve exact four-copy
sequence ties and all 99 protected positions. Cross-scoring against the four
TS templates and matched E–S reference completed with exit zero; variant predictions remain
held for predictor calibration. No additional RFD3 expansion was launched.
The 229-position mutable mask, 99 protected positions and four-copy sequence
ties are prepared. Full PMF convergence is not an extra prerequisite for this
authorized exploratory pilot. The failed Boltz WT geometry is a separate
predictor-validation issue, not evidence that WT is catalytically inactive.
Detailed evidence: [gate audit](../design/nyl12_relay_ts_ensemble_20261126/broader_redesign_gate_audit_20261127.md).

**Completed WT control comparison:** both corrected AF3 jobs finished at
12:25:55 EST. All four models passed software/atom inventory, and independent
protein contacts confirm A/E, B/F, C/G and D/H physical partners. All 16
assigned sites and all 192 amide/ligand/site combinations fail even Tier C.
The protein-only Phe113 N–Asn186 ND2 separation is 9.92–20.34 Å in predictions
versus 3.83–3.88 Å in n12c, establishing a pocket-validation failure independent
of ligand register. These are unrelaxed, query-only/template-free controls;
they do not invalidate the separate MD/QM/MM TS parents or measure activity.
Investigate protein fold/assembly conditioning with homologous protein MSAs
and a protein-only crystal-template control while retaining unguided PA6
placement. No MSA-related cause has yet been established, no threshold is
relaxed. Two matched-seed protein-MSA controls subsequently completed native
inference and canonical coordinate/metric roundtrips; neither yielded an
assigned Tier A/B site. A conditional workflow collects and audits these
outputs before one protein-only crystal-template follow-up per predictor.
Both template controls completed native inference and canonical/Amber
roundtrips. AF3 finished at 13:51:16 EST; its postprocessing was repaired after
the original watcher mistook the normal top-level plus seed/sample CIF copies
for an ambiguous prediction. The retry explicitly selects the native top-level
model and verifies identical atom identities/coordinates against the sole
sample; original failures are retained. Publication finished at 14:24:03 EST,
without rerunning inference. No template control was relaxed.

AF3 now yields four Tier C sites: OG1–C = 2.666–2.695 Å, carbonyl-O–Phe-N
= 2.527–2.538 Å, carbonyl-O–Asn-ND2 = 3.443–3.466 Å. Attack angle alone
(80.48–81.68°) fails Tier A/B (90–120°/85–125°); signed phi is negative.
The Phe-N–Asn-ND2 separation is 3.715–3.729 Å versus crystal 3.831–3.883 Å,
showing recovery of the native pocket. Boltz restores pocket spacing but all
48 amide/site comparisons reject: its attack angles are 56.02–66.27°, Asn
contacts 4.015–4.156 Å and signed phi positive. All internal amides and sites
were searched; the best and nearest assignments coincide in these controls.
The remaining issue is reactive ligand orientation, not the earlier gross
pocket deformation. Thresholds remain unchanged; no activity improvement or
relay-ready geometry is inferred. Full metrics are in SI S10 and
[the template audit](../design/nyl12_relay_ts_ensemble_20261126/broader_pilot/wt_protein_template_controls_20261127/audit_summary_20261127.json).

**Terminology and secondary-feature audit, 27 November:** `activation` means
Thr233 hydroxyl proton uptake by its own N-terminal Nα coupled to nucleophilic
addition; `delivery` means Nα-to-substrate-amide-N proton transfer coupled to
C–N cleavage/acylation completion. These are sampling-segment labels, not
verified elementary steps or enzyme maturation. The 50 ps temporal-half
barriers, 24.877/26.964 kcal/mol, have the same peak location s = 0.6875:
the height shifts by 2.087 kcal/mol. The s ≈ 0.4 feature is a shoulder, whereas
the .5375-to-.5875 hump/trough drop is 1.086 kcal/mol and persists conditionally
in paired bootstraps and both halves. It is nevertheless strongly coordinate
sensitive; 37.4% of trough-region weight in the frozen-45 ps audit projects
exactly onto vertices. The shared window033 TI1-like reference is at
s = 0.537776, not the shoulder. See SI S1/S9 for explicit definitions and
the horizon-specific evidence. No separate TI1 free-energy basin is established.

**Design decision:** retain the main-maximum ensemble and identify a small
separate earlier addition/proton-transfer set as a possible secondary-bottleneck
control within the existing budget. Select by chemistry and structural
clustering, not s alone; do not label every s = 0.4–0.6 frame a TS. Evaluate
whole sequences across the candidate bottlenecks plus E–S/acylenzyme
compatibility. Stabilizing TI1 alone can raise its escape barrier. This
clarification launches no additional campaign and does not bypass backbone,
sampling or TS qualification gates. Current scientific jobs continue.

**Why extra sampling:** hard nearest-polyline projection collects finite
probability at corners, and some early bins depend on single highly weighted
frames. The new 3×11 grid samples q_acyl = −1.25/−1.50/−1.75 Å and
q_act = −1.0…1.0 Å by 0.2 Å beside the original center line. Both coefficients
remain 25 kcal mol⁻¹ Å⁻² in Amber's U = kΔ² convention. This does not directly
repair the independent slow rearrangement in activation windows 26–28:
HG1–leaving-N mean changes +1.53–1.62 Å between the first/last retained 5 ps,
with zero histogram intersection, while HG1 remains Nα-associated.

**Six additional checks:** A026 and A028 each receive independent original
PES/near and late-MD/far starts at the same target bias; D008 starts from a
late D010 leaving-N-associated HG1 structure, and D009 from late D007 with
HG1 associated with Nα. All use 5 ps equilibration plus 30 ps production,
unchanged Hamiltonian, masks and proton inventory. Original/far HG1–leaving-N
seed distances are 2.077/3.792 Å (A026) and 2.001/3.345 Å (A028).
Actual D010 donor distances are 1.200 Å to leaving N and 1.452 Å to Nα.
Host-specific Sander/library identities are frozen separately. Geranium’s
Sander and all 16 used 3ob files match Violet. The native 20-step D008 test
passed all SCC, finite-temperature/checkpoint and native-bias parity checks.
The exact idle Azalea queue was retired before six Geranium CPU workers
launched; no GPU allocation or other user process changed. Prior unlaunched
Violet/Azalea preparations and the queue-cancellation receipt are retained.

**TS-region assessment:** the latest completed 50 ps reconstruction has a
25.910 kcal/mol sampled maximum at s = 0.6875: 738 frames, Kish weight ESS
362.15 and conditional pointwise 95% interval 25.286–26.495 kcal/mol. Its
trajectory-half maxima are 24.877/26.964 kcal/mol, a 2.087 kcal/mol difference.
The detailed frozen-45 ps geometry audit places peak support mainly in
D007–D010: neighbor overlap 0.22–0.25 and five-distance early/late shifts
≤0.063 Å (histogram intersections 0.77–0.96). These are useful exploratory
TS-like geometries, not a verified equilibrium transition-state ensemble.
Local stability does not establish the barrier or reactant normalization.

**Analysis/acceptance:** retain a separate 5 ps production discard per
independent lineage; continuations retain the original lineage. Cross-evaluate
all actual native biases and retain all valid proton configurations. Examine
differentiable path-coordinate/physical-2D sensitivities, local/global overlap,
frame-weight concentration, chemical retention, independent-start agreement,
and temporal stability. New final-analysis uncertainty will include
count-preserving 5/10 ps circular moving blocks; the former 5/2.5 ps
nonoverlapping-block intervals remain historical conditional diagnostics.
Apply the same final estimator to the direct route before quantitative
mechanism comparison. Do not change discard/coordinate to obtain a desired
barrier. An operational precision goal is changes below about 1 kcal/mol
across successive late checkpoints and independent-start comparisons, with
uncertainty covering slow modes; this is not a universal convergence proof.
No automatic primary promotion or converged-barrier label is authorized.

**Schedule, measured at 10:38 EST on 27 November:** Azalea runs at roughly
1.39–1.49 s/MD step; Violet's new grid spans 2.18–2.88 s/step. Main 60 ps is
expected about 19:00 EST Friday; new-grid 10 ps around Saturday late morning;
new-grid 30 ps Sunday afternoon/evening. The six Geranium checks started at 11:03 EST Friday without the Azalea wait;
the provisional 24–30 hour native estimate places completion Saturday midday
to late afternoon, subject to steady-throughput refinement. These are throughput estimates, not
convergence guarantees. The earliest comprehensive convergence assessment is
Sunday night/30 November. A publication-ready PMF cannot yet be promised by
30 November; exploratory TS-ensemble design can proceed separately once the
structural qualification gates pass.

**DGX qualification, 27 November:** `bsd-a100` is reachable from Violet after
PI enrollment of the existing public key. Verified hardware is eight A100
40 GiB GPUs, 128 physical CPU cores/256 threads and about 1 TiB RAM. Other-user
work occupies GPUs 0–1; select idle GPUs 2–3, **at most two**, by pinned UUID
and refreshed occupancy checks. Existing Boltz 2.2.1/Torch 2.9 CUDA12.8 passes
native small-kernel tests on both; installed AF3 imports correctly. No software
or model downloads were needed. Home/root had only 16.7 GiB free and the PI is
cleaning it; use existing weights and a bounded staging footprint.

Prepare the two original WT seeds in both predictors as host-qualified,
unguided baseline controls: four WT predictions within the existing eight-model
comparison budget, serialized by GPU. The redesigned-variant gate remains
held until native backbone acceptance; no new variants or long MD campaign
are authorized by this control-first scheduling choice. Existing DFTB3/MM
and US stay on their hosts and continue unchanged. Sander/pmemd are not
qualified on the DGX; GPU availability does not accelerate the current CPU
Sander implementation. Paths, exact runtime identities, storage and launch
constraints are in [the DGX runtime note](BSD_A100_RUNTIME.md).

**Live evidence:** all paths below are within
`qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/`.
`analysis/activation_support_preparation_20261127/priority_group/` records the
37 running priority tasks and guarded pause/restoration of the 40 direct-US
and two hydration processes. `analysis/pooled_support_20261127/` watches the
controlled frozen-45 ps + 37-supplement diagnostic (94 independent lineages),
which does not replace the final most-extensive analysis. `mixing_checks_20261127/`
and `analysis/mixing_checks_preparation_20261127/` hold the six queued inputs,
actual seed geometry and host-specific execution. `analysis/final_candidate_20261127/` holds the active final100 watcher, with
an isolated Geranium-derived mirror and explicit host/hash binding so stale
Azalea metadata cannot enter the fit. `analysis/ts_sampling_audit_20261127/`
and `plots/reaction_distance_statistics_45ps_20261127/` retain the quantitative
TS/distance evidence and source hashes. Adding the first four short checks
improves the last two bins' weight ESS to 54.81/15.01, without proving equilibrium.

**Design qualification:** both RFD3 native runs finished. Fixed-coordinate
discrepancies are BF16 quantization; a residue-map error and false Val
methyl-label stereochemical flag were diagnosed. Only an in-memory precision
normalization makes the 5 Å arm pass; the 2 Å arm still has two Arg118 clashes.
No raw output or acceptance gate changed and no broad sequence/prediction
campaign has launched. Evidence is in
`design/nyl12_relay_ts_ensemble_20261126/backbone_qualification/diagnostics_20261127_v2/diagnosis.json`.

The dated entries below are historical checkpoints and retain their original
run states; this entry supersedes their active-status descriptions.

---

**Execution and 35 ps update, 26 November (late evening):** all 57 relay windows
completed 35 ps and continue. The common fixed-discard/200-bootstrap analysis
gives a sampled maximum of **25.445 kcal/mol at s=0.6875**, with conditional
pointwise 95% interval **24.758–26.070**. The temporal-half maxima are
23.976/26.575 kcal/mol; sampling is **not converged** despite connected global
overlap. Low support persists near s=0.21–0.34 and the product tail. Four
targeted diagnostic runs remain separate from the main 57-window data set.
The receipt-discovery bug that picked up nested diagnostic receipts is fixed.
`relay_auto_40to60_20261126` now watches common 40/45/50/55/60 ps, using four
analysis workers and automatic PMF/movie publication; it does not repeat 35 ps.
Movie `full_path_300K/vmd_progress/current/view_path.vmd` now contains
**3,990 frames (70/window)**, with checked hashes and display-policy probes.

The two RFD3 qualification targets are actually running/queued sequentially on
**Violet GPU 0** in `relay_rfd3_qualification_retry01_20261126`. The first native
attempt failed before diffusion because P6N/P6I/P6C collided with unrelated
installed CCD entries. A reversible input-ID alias repair passed the complete
native CPU feature pipeline; no ligand atoms, typed bonds or coordinates were
changed. The model omits eight protein-terminal OXT atoms by construction;
raw outputs are retained and the declared adapter restores only these fixed
source atoms after checking their retained terminal anchors/alignment.
Structural qualification is **pending**, not already accepted.

`relay_broad_gated_20261126` is actually queued and waiting on both native
qualification acceptances. It verifies source/output hashes, invokes the safe
release-only expansion gate, then runs at most eight total backbone outputs,
24 broad sequence samples and the eight-model unguided AF3/Boltz comparison.
The frozen broad mask permits **229/328 precursor positions** (172 alpha / 57 beta),
protecting 99 including initiator M1; exact four-copy sequence ties remain.
The fixed local scoring subset has 18 reactive-pocket positions; global 229-site
scores are kept separate. Current inputs, source snapshots, tests and live
status are under `design/nyl12_relay_ts_ensemble_20261126/`:
`backbone_qualification/` and `broader_pilot/controller_job/`.
No expanded generation/prediction is claimed completed before the gate passes.
The current MVP is sequential; no pinned source was altered to add parallelism.
Assay references R61/R62 and explicitly prospective TNBS methods are now in
`manuscript/REFERENCES.md` and `manuscript/METHODS_DRAFT.md`.

**Current design/assay revision, 26 November:** the PI approved general PA6
hydrolysis rather than dimer specificity, requests much broader mutable
sequence space and proposes TS-conditioned RFD3 partial diffusion before
LigandMPNN. The five-site pilot below is retained as a completed baseline.
The revised plan in `docs/design_ranking_and_iteration_20261022.md` supersedes
prospective contact-guided/single-register selection: unguided AF3/Boltz
validation, every chemically mapped PA6 amide eligible, exact copy-sequence
ties, protected catalytic/maturation/interface core, no initial indels.
Proposed ceiling: two RFD3 qualification outputs, expanding to eight outputs
total only if accepted; up to four remodeled parents and 24 new broad-mask
sequence samples split between original/remodeled backbones. At most 30 new
predictions total includes AF3-versus-Boltz qualification and confirmation;
downstream remains two designs+WT with two MD replicas each, six QM/MM
optimizations, three RCD paths and US only for the strongest design+WT.
No new production campaign or installation was launched by this plan update.
The current relay US and sampling checks remain unchanged.

RFD3 source inspection supports nonsymmetric partial diffusion in full mature
tetramer context followed by exact sequence tying; geometric symmetry and
sequence symmetry are distinct. Custom four-residue PA6 input/graph retention
and full GPU denoising remain qualification gates. The inspected
`redesigner-agent` checkout demonstrates LigandMPNN symmetry on merged
328-residue chains, not a retrieved RFD3 execution artifact. Preserve the
current eight-chain mature map and authentic chain breaks.

**Experimental interpretation:** the PI reports current calibrated detection
of L2/L3/L4 only; larger oligomers can be missed and assays are likely to remain
in crude lysates for initial screening. **The PI now confirms that purification
of WT plus a few promising designs should be feasible**, and favors adding a
primary-amine-release assay even if reaction buffer must change. Plan WT+2–3
finalists for purified validation after qualifying the assay on WT.
Optimize/report short-product yield as a practical outcome
without equating it to total cleavage or intrinsic catalytic efficiency.
Record expression, processing and stability confounding. Retain lysates for
throughput and the small purified panel for cleaner comparison; cell-free
synthesis remains an optional later platform. Prioritize TNBS qualification
alongside calibrated MS: an amine-free reaction buffer (phosphate candidate),
WT activity/stability and paired MS across old/new buffers, recorded pH versus
temperature, enzyme-only/substrate-only/time-zero/reagent controls, maturation
checks, calibration/quantification limits and pre/post-cleanup spike recovery.
Protein removal and non-nylon consumables need validation; report any size cutoff.
These readouts do not measure all cleavage or intrinsic kcat/Km. Consider
broader LC–HRMS, validated amine-release
analysis (Tris/lysate interference must be removed) and SEC/GPC as complementary
readouts; analytical post-hydrolysis of the isolated soluble fraction is an
unvalidated option requiring recovery/completeness controls. See the detailed
plan's “Assay constraints change what success means” and “Primary-amine assay
development” sections and references. No new assay has been validated/launched
by this documentation update.

**Bounded TS-ensemble redesign, 26 November (PI authorized):** tmux
`relay_ensemble_generation_20261126` completed on initially idle Violet GPU 0,
sampling eight sequences from each of three actual relay-US parents at
s=0.672/0.725/0.745 (24 raw total). Mutable canonical positions are
V111/I112/D114/A187/V188; beta chains remain fixed. Preserve exact four-copy
sequence symmetry, catalytic/maturation identities and inter-repeat interfaces;
use the occupied repeat for sampling weights, retaining the entire assembly.
No duplicated artificial TS ligands. Full parent/sequence/source provenance
lives in `design/nyl12_relay_ts_ensemble_20261126/`.
Native exit 0 and sequence symmetry/protection acceptance passed: **four unique
sequences** (A187S/V188I; A187S; I112L/A187P; I112L/A187S), with complete raw
lineages retained. Cross-scoring all three parents plus matched E–S also
completed in `relay_ensemble_crossscore_20261126` with exit 0 and acceptance;
A187S leads the descriptive sequence-model ranking, not a validated activity
ranking. The actual four-design pool reduces initial Boltz work to
at most 10 predictions including WT; confirmation would bring it to at most
22, within the larger ceiling below. No arbitrary quota-filling generation.
The detailed resource-capped funnel is in
`docs/design_ranking_and_iteration_20261022.md`, dated section
“Bounded relay-ensemble redesign campaign”: up to 8 designs+WT initial Boltz
screen, top3+WT confirmation, **30 predictions total**; at most 2 designs+WT
with two 100 ns MD replicas each; at most 6 QM/MM optimizations and 3 matched
RCD paths; US/MBAR only for the strongest design and matched WT. These are
ceilings, not guaranteed advancement or an activity ranking. Only the small
sequence pilot is launched now; follow-ups require actual selection results.
Contact guidance aids model construction and is removed for relaxation;
guidance-free sensitivity seeds and matched WT controls remain mandatory.
Keep SolubleMPNN surface redesign and dimer/exo specificity as subsequent
separate objectives. Current objective is improved PA6 cleavage with product
profiles measured, without assuming the modeled L4→L2 event is established
experimentally. The 57 main relay windows and four targeted tests continue
under their existing scheduling and acceptance rules.

**Relay 30 ps and retrieval update, 26 November:** all 57 Azalea windows are
already in stage09 (30→35 ps; last measured production 30.90–30.95 ps), continuing
toward 60 ps. The fixed-cutoff 30 ps PMF sampled maximum is 25.042 kcal/mol
(pointwise conditional 95%: 24.285–25.697); temporal-half maxima 23.578/26.149
remain inconsistent. The updated stripped movie has 3,420 frames (60/window)
at relay `full_path_300K/vmd_progress/current/view_path.vmd`. A broad compact-JSON
sync was overwriting its new metadata with stale remote 20 ps metadata; only
the derived synchronizer was repaired. Bundle checksums remained valid over
subsequent sync cycles. Native simulations were untouched.

Endpoint diagnosis: nearest-polyline projection maps about 17% of sampled
frames exactly to vertices, including 32% near s≈0.1625. This contributes to
persistent early structure, but does not explain every wiggle. The last two
bins extend beyond s=1: their weight ESS values are 1.25 and 1; the final bin
contains one frame and has unresolved CI. These are not established acyl-basin
features. The separately retained, low-priority one-refit diagnostic reproduced the
primary PMF exactly. Vertex projections carry 32.7% of the s≈0.1625 bin weight
but only 0.203% near s≈0.1875; the raw 17% global vertex fraction is therefore
not a weighted explanation of the whole profile. One frame contributes 89.0%
of the penultimate bin weight and 98.4% at s≈0.2375. The original PMF/grid/biases
remain unchanged; new binning sensitivities are diagnostic only.
See `analysis/fixed_burnin_20261126/endpoint_diagnostics_20261126/`.

A separately specified differentiable path-coordinate diagnostic removes exact
corner coincidences but retains early oscillations; it does not replace the
primary coordinate or revise the barrier. A subsequent endpoint follow-up finds
that three delivery-026 frames carry 90.46% of the globally normalized unbiased
MBAR weight; the
penultimate and terminal bin weight ESS values are 1.248 and 1.000, so a product
minimum is not established. Activation windows 021/022 have biased 2D histogram
intersection 0.556 but mean unrestrained q_del values −2.663/−1.913 Å; reciprocal
starts are needed to assess relaxation. The activation-030/delivery-001 junction
remains connected at minimum overlap 0.110; this is not a failed junction.

Four reciprocal/targeted windows were launched with 1 ps initial equilibration
using fresh 300 K velocities, 4 ps further equilibration, and 10 ps diagnostic
production:
`activation_021_from_022`, `activation_022_from_021`,
`delivery_027_qdel_1p50`, and `delivery_028_qdel_1p75`, under
`full_path_300K/additional_sampling_20261126/`. They were still in first-1-ps
equilibration at
the latest status. Any new production windows require matching at least 30 ps
before final PMF comparison. The full direct-US set (40 native engines) is
paused together, as are two hydration engines; the 57 Azalea relay workers were
left untouched. A tested independent monitor is waiting for all four targets to
finish before resuming the 42 direct/hydration engine identities together;
group restoration is queued, not yet completed.

The Marigold retrieval benchmark completed: current lexical evidence delivery
20/20, dense19/20, reranking20/20; ordering improved but no coverage gain. Current
retrieval stays unchanged, model servers remain off, and BioReason-Pro 2 is
explicitly excluded. The PI selected LigandMPNN; no RFD4 pilot is queued.
Mechanism-agent benchmark and expanded NVIDIA collection assessment are saved
in their existing project directories. No scientific dispatch was enabled.

---


## Current checkpoint — 26 November 2026

- **Mechanism-agent readiness:** Marigold and Lavender are reachable and idle;
  no model servers were started during this check. Lavender requires the
  benchmark's retained `deployment/nano_known_hosts` with strict host checking;
  its missing default-host record is not a reason to bypass verification.
  The 20-case comparison is complete: Nemotron is qualified for supervised
  read-only evidence/citation synthesis; Mistral failed chemistry qualification.
  The paired RAG pilot completed 20 responses and the lead factual review:
  9 correct, 5 minor and 6 major findings; no critical findings. Valid citation
  checks did not catch route/release identity, historical-status and water-metric
  conflation. These include retrieval/curation failures, not only model errors.
  The original answers are retained; the stale top-level status is reconciled.
  The refreshed index has 17 curated primary sources and a separate 29-source
  reference view; six problematic natural-query source groups are now retrieved
  in the top six. One supervised hydration-report test completed with controller
  exit 0, a TEST-only awaiting-human packet, null approval and one trusted mail
  submission. Lead review accepted the core observations/proton ledger but
  requires corrections for absence/status overclaims. Both Nano model servers
  and the local tunnel are verified off; 43 Violet and 57 Azalea Sander processes
  continued. Validation: 157 existing tests and 13 new guards passed, with 88%
  runner coverage. Scientific dispatch and autonomous planning remain disabled.
  No response to the test email is needed. See mechanism-agent's
  `runs/supervised_handoff_prep_20261126/completion.json` and its lead review.
  The fresh 10-scenario paired curated/retrieved scientific-decision benchmark
  is preparation-only pending lead confirmation; it does not change model
  authority or the November 30 PMF/TS priority. See mechanism-agent's
  `docs/ROADMAP.md` and `docs/QUALIFIED_SCIENTIFIC_CONTROL.md`.

- **Live launch verification:** Violet has 40 native direct-US workers plus
  three native hydration runs; Azalea retains 57 native relay-US workers.
  Violet has 24 physical cores / 48 logical CPUs. These 43 single-thread jobs
  share SMT resources; hydration has lower scheduling priority than US.
  Direct windows continue in their existing folders as production blocks
  07–12, under `continuation_60ps_20261126/plan.json`; no new NEB was started.
  The derived-only `direct_us_matched_pmf_20261126` watcher analyzes complete
  35/40/45/50/55/60 ps checkpoints with the common fixed-discard/MBAR/bootstrap
  method, then publishes comparisons when the matching relay horizon exists.
  It uses four CPU workers when fitting; current phase is waiting for 35 ps.
  Native original inputs/30 ps data are retained and new SKF/engine/topology
  identities are pinned explicitly. Source manifests, reviewed plans and
  tmux/process/scientific receipts are retained in the same campaign.
  Hydration now comprises the main 51→101 ps continuation and two additional
  25 ps stochastic branches from accepted 26/41 ps checkpoints. Their plans,
  movies and separate new-segment statistics are under hydration
  `replicas/replica_02_from_26ps_extend04/` and
  `replicas/replica_03_from_41ps_extend07/`. Shared prehistory is not independent
  sampling. All three stop for a candidate water or a failed acceptance gate.
  Earlier hydration blocks took about 4.1 h per 5 ps; the initial timings under
  the new concurrent load imply roughly 5–8 h per block, with substantial
  startup uncertainty. No resources are assigned to verbena,
  tulip or geranium. The pull script includes curated new outputs.

- **PI-approved execution through 30 November:** relay convergence and a
  defensible transition-state ensemble for redesign remain the primary goals.
  Keep all 57 relay windows on Azalea and the queued 60 ps continuation running.
  Extend the 63 existing direct windows from 30 to 60 ps on Violet using at most
  40 one-CPU workers, prioritizing completion of each common 5 ps checkpoint.
  Continue accepted coordinates/velocities with the same native biases,
  proton inventory and Hamiltonian; preserve all earlier inputs and samples.
  Retain Violet CPU capacity for relay analysis/initialization checks and the
  three bounded one-CPU hydration streams. No new NEB. Tulip and geranium are unavailable;
  no work is assigned to verbena while the two existing hosts suffice.
  Native launches have now been verified above; completion of 60 ps is not
  the convergence acceptance gate.

- **Initialization checks explained:** direct US used NEB seeds; relay US used
  two coordinate-driving branches joined at the identical window033 restart.
  That seed-generation difference is acceptable; we use one global MBAR fit,
  not independently zeroed PMF pieces. After extended-sampling diagnostics,
  select a few problematic windows and run each at the SAME bias center from
  two structures approached from opposite sides. Compare whether both settle
  into consistent chemical and environmental distributions. This tests memory
  of initial conditions, not reverse reaction rates. Prioritize the slow
  activation025/026 region, chart junction and weak acyl basin if those issues
  persist. Do not launch another full independent campaign automatically.

- **Deacylation hydration:** 51 ps accepted total (1 ps heating + 50 ps at
  300 K), with acyl chemistry retained. No near-attack water was found; the
  minimum water-O–acyl-C separation was 4.939 Å at 30 ps, angle 72.1°.
  The final nearest separation was 6.314 Å. PI approved a one-CPU continuation
  on Violet to 101 ps total in ten 5 ps blocks, stopping early for a candidate
  or failed numerical/chemical/boundary gate. Waters remain MM during this
  reconnaissance; a candidate must be promoted to QM and reoptimized before
  deacylation chemistry is modeled. At the current endpoint both the Thr
  N-terminus and product amine are neutral NH2; HG1 resides on the product
  amine. The earlier derived annotation calling the N-terminus protonated is
  stale and must be corrected explicitly in the new analysis provenance.
  If 101 ps produces no candidate, inspect product occlusion/water access
  before commissioning another blind extension. Completing the second
  half-reaction by the deadline depends on obtaining a defensible hydrated
  starting state; do not compromise the primary relay/PMF work to force it.

- **Deacylation and product protonation — hypothesis update, 26 November:**
  the retained neutral product NH2 is a valid acylation endpoint; its eventual
  solution protonation is a separate question. NH3+ should predominate when
  solution pH is below the conjugate-acid pKa, but the bound-state pKa, donor
  access and timing are unresolved. Do not assign an active-site protonation
  state from a bulk primary-amine pKa alone, particularly across temperatures.
  Keep the current hydration runs and proton inventories unchanged.

  A closed, atom/charge-balanced cycle can end as:
  `E–O–C(=O)–R + H2O + R′NH2 → E–OH + RCOO− + R′NH3+`.
  Here R/R′ retain their original polymer ends and spectator charges. HG1 has
  already supplied one product-N hydrogen during acylation. The incoming
  water's two hydrogens can ultimately supply Thr–OH regeneration and the
  additional product-amine proton; the water oxygen supplies the new carboxyl
  oxygen. This ledger does not establish a transfer order or unique proton
  identity in an exchanging hydrogen-bond network.

  Compare two hypotheses only after a suitably hydrated starting state exists:
  **(i)** a base activates water during attack, then its conjugate acid returns
  a proton to Thr OG1 as the ester cleaves; the resulting carboxylic acid can
  subsequently protonate the product amine, directly or through water/buffer;
  **(ii)** Asp272 or Asp274 accepts a water proton during attack and transfers
  it to the retained product NH2. In (ii), explicitly identify how the other
  water proton, potentially on the nascent carboxylic acid, reaches Thr OG1.
  Require regeneration of the selected catalytic microstate, or account for
  any residual protonation change; do not leave an unexplained Thr alkoxide
  or consume a catalytic acid/base stoichiometrically. Asp involvement remains
  a hypothesis, not a result inferred from proximity or pKa alone.

  The baseline is Thr233 Nalpha-mediated water activation and proton return,
  following the substrate-hydrolysis precedent in
  [Andjelkovic et al. (2023; R41)](https://doi.org/10.1021/acs.jcim.3c00900).
  Screen accepted hydration checkpoints for water–Nalpha/Asp272/Asp274 and
  product-amine relay contacts. Compare the baseline with the best supported
  Asp route from the same hydrated snapshot before expanding alternatives.
  Allow water deprotonation and attack to proceed together rather than
  requiring an isolated hydroxide minimum. Promote the attacking water and
  any transferring relay water to QM before reactive calculations; present
  MM waters cannot report these proton transfers.

  Keep a per-state O/H/charge ledger through water attack, ester cleavage,
  Thr–OH/base regeneration and product protonation. Internal redistribution
  conserves the closed-system total charge. Adding a proton from buffer or
  comparing states with different proton counts requires an explicit proton
  reservoir/pH thermodynamic leg and consistent countercharge treatment.
  Final bulk protonation may occur after product release; it is not required
  to occur inside the pocket to complete catalytic bond cleavage. No barrier
  or catalyst-regeneration claim follows from this proposed scheme.

- **Direct-route reanalysis completed:** both mechanisms now use the shared
  fixed-5-ps-discard/all-frame MBAR estimator, 3-CV projection convention and
  200 contiguous-block bootstrap refits. At the matched 25 ps horizon, sampled
  maxima are 50.920 kcal/mol for direct-NEB-seeded sampling versus 24.807 for
  relay-seeded sampling. The full 30 ps direct estimate is 50.947 kcal/mol.
  The direct maximum has only 15 frames; bootstrap support is 194/200 at 25 ps
  and 195/200 at 30 ps, so no finite 95% band is assigned at that bin. Neither
  curve is a converged activation barrier. Direct HG1–Nα nearest-owner visits
  occur in about 11% of its retained biased frames; this is not a pure direct
  mechanism ensemble. The old 47.05 kcal/mol direct result used a q_acyl
  marginal and is historical, not a like-for-like path-density estimate.
  Compare `plots/comparison_25ps_vs_relay_25ps_20261126.png`, available in both
  campaign plot directories; the direct `plots/pmf_latest.png` shows all 30 ps.
  Source-bound reports, model audit and summary are under the direct campaign's
  `analysis/fixed_burnin_20261126/`. Native US and its 60 ps continuation were
  untouched. Next: assess the extended relay data and independent/basin support;
  do not restart the expensive historical direct campaign solely for a smooth
  figure or interpret these maxima as a mechanism-specific rate ratio.

- **PI-required matched comparison:** both mechanisms must use the same
  versioned preparation, physical model, analysis rules and acceptance gates.
  Burn-in means an initial analysis discard; native data are preserved. At the
  current checkpoint, 25 ps production minus 5 ps additional discard gives
  20 ps analyzed per window, after separately excluded heating/equilibration.
  Changing the postprocessing rule requires reanalyzing every compared path.
  The direct campaign remains a historical control with matched postprocessing:
  topology, selections and frozen-coordinate identity have now been checked,
  but historical direct-era 3ob hashes are unavailable and SCC tolerances differ.
  Each route retains its physical reference path; sampled channels are mixed.
  Different numerical paths/window counts can be necessary; equal runtime
  alone is not matched convergence. The controlling specification is
  `qmmm_umbrella_sampling_and_pmf_protocol_20261109.md`, section "Matched pathway
  comparisons and publication acceptance". Next: complete 60 ps, diagnose
  windows025/026 and basin bracketing, then choose targeted bridge/independent-
  start tests from those diagnostics. Do not start another broad campaign.
  The new all-frame E–S/acyl minima are numerically interior to the declared
  search intervals and have occupied neighbors. The product minimum remains
  statistically weak (27 frames, local weight ESS 1.24; next bin one frame).
  Do not reuse pass/fail labels from the older automatically thinned estimator
  as though they were recalculated for this estimate. The next endpoint task
  is to establish adequate product-basin support/stability, not to fill a
  purely geometric gap that the all-frame histogram no longer has.

- **Primary PMF protocol revision:** use an explicitly declared additional
  5 ps production burn-in once per continuous window, all subsequent frames,
  and 200 within-window contiguous 5 ps block-bootstrap MBAR refits for
  pointwise 95% intervals. Native 1 ps heating and 5 ps equilibration are
  already excluded. Compare 0/2.5/10 ps discards, early/late halves and 2.5 ps
  bootstrap blocks. Neither the fixed discard nor the intervals prove
  equilibration. Do not choose the analysis giving the lowest/smoothest PMF.
  Archive the automatic-cutoff estimates below as sensitivity results.
  At 14:15 EST all 57 windows had completed 25 ps and moved to stage08;
  the corrected Azalea guardian had completed all stage07 derived records.
  The source-frozen `relay_fixed_pmf_parallel_20261126` session on Violet
  analyzes each new complete checkpoint using eight idle CPU cores. It has
  no simulation-control interface. Source/results:
  `full_path_300K/analysis/fixed_burnin_20261126/`; figures remain in
  `full_path_300K/plots/`. Both previous derived publishers were deliberately
  retired; their logs and frozen source remain. Keep native US, receipt
  validation, synchronization and the queued 60 ps extension running.
  Acceptance still requires temporal/orthogonal stability, overlap, adequate
  effective sampling and independent initialization checks. Use the same
  analysis protocol for the eventual matched Mechanism A/B comparison.
  **Completed 25 ps result:** 22,800 post-cut frames; highest sampled PMF
  24.807 kcal/mol at s=0.6875, conditional pointwise 95% interval
  23.938–25.553 kcal/mol. There are 43 fully bootstrap-supported bins out of 45
  populated bins and no disconnected resamples. However, the 5–15 and 15–25 ps
  profiles peak at 23.008 and 26.318 kcal/mol; statistical bands do not capture
  this temporal instability. Windows 025/026 have leaving-N–HG1 correlation
  scales g×dt slightly above 5 ps. The global target weight ESS is 1.72
  (not temporally corrected). No converged barrier is accepted. Inspect
  `plots/relay_25ps_fixed_burnin_20261126_pmf.png` and the corresponding
  `_sensitivity.png` / `_block_sensitivity.png`; `plots/pmf_latest.png`
  now points to this result. Keep 60 ps sampling and prioritize slow-coordinate,
  endpoint-bracketing and independent-initialization checks before rate claims.

- **Late-cutoff audit:** the 15–18.2 ps discards are not certified physical
  equilibration times. Taking the latest optimum across 13 observables and
  then applying the largest common-suffix g is overly aggressive in several
  examined windows. Window002 has stable five-ps means but a 17 ps fast-search
  cutoff, moving to 6.85 ps with full lag summation. Window023 shows an initial
  5–6 ps rearrangement and a later excursion; 18.2 ps isolates a final quieter
  suffix. A 5 ps sensitivity discard retains 18 versus 6 configurations there
  under the same common-g thinning. Window024 has a genuine late shift in
  the biased acylation coordinate, so a universal early cut is also unqualified.
  HG1 ownership and both original N-terminal N–H bonds remain intact throughout
  these eight series. Do not silently alter the running analysis or choose
  cutoffs to obtain a desired barrier. Before quantitative acceptance, compare
  prespecified burn-in sensitivity and temporal/orthogonal stability using the
  extended trajectories. Full source-bound plots and interpretation:
  `RELAY_EQUILIBRATION_AUDIT_20261126.md`. At 13:41 EST all 57 Sander processes
  remained active in the fifth segment; full 25 ps coverage is forecast around
  14:10 EST, followed by derived validation and automatic PMF analysis.

- **Archived automatic continuous-window PMF, superseded later 26 November:** the PI authorized reanalysis
  after concatenating each window's uninterrupted production segments, with
  one equilibration discard and a correlation factor estimated over the
  complete retained series. On the same 20 ps/window data, this retains 2,840
  frames instead of 3,150. At s=0.6875 the largest sampled value changes from
  22.154 to **25.320 kcal/mol**, relative to the unchanged sampled E–S plotting
  reference. Matched analytical pointwise 95% intervals are 20.924–23.384
  versus 23.919–26.721 kcal/mol. The maximum absolute curve difference over
  common occupied bins is 4.523 kcal/mol. This is sensitivity to analysis
  policy, not an independent-replicate significance test or a validated ΔG‡.
  Global target weight ESS changes from 6.73 to 29.29; 46/57 windows still
  retain fewer than 100 selected frames. Eight windows retain less than 5 ps
  after their single equilibration cut (minimum 1.80 ps).
  Supplementary 200 within-window bootstrap refits of the globally decorrelated
  frames give 24.119–26.762 kcal/mol at the maximum's coordinate; 37 bins have
  support in every replicate. These bootstrap units differ from the archived
  5 ps block resamples, so the main comparison uses matched analytical errors.
  Both estimators remain nonconverged and the basin-reference gate unresolved.
  Source-bound results are under `analysis/continuous_window_20261126/`;
  `plots/comparison_20ps_continuous_vs_blockwise_20261126.png` shows both curves.
  The native US workers were unchanged. The source-frozen derived session
  `relay_continuous_pmf_dynamic_20261126` published this selection method and
  is now deliberately retired. Its source, results and termination record
  remain evidence. Keep the native-data sync, accepted-receipt guardian and
  cumulative MBAR watcher running; use the fixed-burn-in publisher above.
- **Sampling extension prepared and queued:** bounded continuation from 30 to
  **60 ps/window**, in the same logical campaign/window directories and 5 ps
  production segments. Continue checkpoint velocities, thermostat, Hamiltonian,
  biases, frozen region and timestep without reheating. New workers must wait
  for all 57 baseline 30 ps endpoints and their accepted lineage; do not edit
  the live plan or merely increase its pinned block count. Record original-plan
  identities for stages through 08 and extension-plan identities thereafter.
  All **342 new production inputs** and **342 accepted historical stage
  lineages** passed read-only preflight on Azalea. Inputs match each window's
  stage08 input byte-for-byte after normalizing only the thermostat seed;
  seeds are distinct and recorded. Restraints, topology, QM charge, method,
  parameters and movable region are unchanged. New native stages are 09–14
  in the existing window directories. The waiting session is
  **`relay_us_cont60_wait_20261126` on Azalea**. It starts at most 57 CPU
  workers after all 57 stage08 endpoints and their strict data receipts are
  accepted and no campaign Sander process overlaps; it runs to 60 ps/window
  and stops. Its wait is bounded to 24 hours. Preparation launched no extra
  Sander process; all 57 baseline Sander processes remained active.
  The authoritative extension controls are
  `continuation_60ps_20261126/gate_status.json`, `plan.json`, `root_review.json`
  and `preflight.json`. Compact synchronization and cumulative PMF analysis
  have no 30 ps cap. The frozen runtime inspector retains the baseline
  controller pointer; consult the extension controls for its execution identity.
  A separate derived-only guardian uses the validated trajectory-time precision
  repair for stages09–14. Approximate additional runtime is one day after
  baseline completion, based on about four wall-hours per 5 ps block.
  The 60 ps target is a sampling budget, not a convergence criterion.

- **PMF reporting audit, 26 November:** the earlier PMF uncertainty bands were
  generated by block-bootstrap MBAR refits; the relay checkpoint had not run
  uncertainty estimation. This was a reporting omission, not zero uncertainty.
  New source-bound analysis recovers exactly the original 3,150 selected frames,
  path projections, bin edges and E–S reference, then calculates PyMBAR histogram
  covariance and **200 within-window completed-5 ps block-bootstrap MBAR refits**.
  All 43 occupied bins have support in all 200 replicates. The largest sampled
  PMF value is 22.154 kcal/mol relative to the fixed plotting reference, with
  pointwise bootstrap limits 20.642–23.326 kcal/mol; this is **not an accepted
  activation barrier or a confidence interval for the maximum**. The source
  acceptance gates remain unchanged. In particular, the bin near s=0.3125 has
  81 selected frames but bin-local weight ESS ≈1.02; s=0.1125 has five frames
  and weight ESS ≈2.56. More raw frames do not necessarily give effective
  reweighted sampling. Use the shared plotting/uncertainty entry points and
  `PRESENTATION_AND_FIGURE_STYLE.md`, with unconnected overlap-pair markers.
  The reported 31→57 line was a plot-order artifact from an extra state30–32
  comparison appended after regular pairs; it was not an MBAR network edge.
  Preserve the current 30 ps campaign and reassess local support, time stability
  and orthogonal mixing before deciding on targeted extensions or new windows.
  A separate, **non-adopted selection sensitivity check** concatenated the four
  production blocks/window before estimating equilibration/correlation: median
  statistical inefficiency increased from 3.08 to 6.73, maximum from 11.73 to
  33.23, and selected frames decreased from 3,150 to 2,840. The current blockwise
  estimator can miss correlations across segment boundaries. Before accepting
  a final barrier, analyze each continuous window with one documented production
  burn-in and a correlation estimate spanning its complete retained series;
  compare estimator/reference stability while retaining the original checkpoint.
  Evidence: `full_path_300K/analysis/pmf_uncertainty_20261126/selection_sensitivity.json`.
  Standard figures are now published as `plots/relay_20ps_standard_20261126_*`
  and stable `plots/pmf_latest.png`, `plots/mbar_overlap_latest.png`, and
  `plots/neighbor_overlap_latest.png`. The source-frozen, derived-only tmux
  session `relay_pmf_standard_figures_20261126` automatically applies the same
  200-refit uncertainty and plotting protocol at subsequent completed 5 ps
  checkpoints for up to 72 hours; it does not inspect or restart native workers.

- **Actual sampling at 10:03 EST:** all 57 relay windows have completed four
  accepted 5 ps production blocks: **20 ps/window, 1.14 ns aggregate**. Sampling
  in the fifth block continues. After the repair below, cumulative MBAR now
  includes all 20 ps/window. The original internal gap filled by 10 ps, but
  the 15/20 ps reactant minimum lacks occupied neighboring bins and effective
  sampling/overlap gates still fail. No converged barrier is accepted.
- **Analysis recovery is running:** later blocks failed derived trajectory-time
  validation in the original worker snapshot. The corrected guardian had exited
  while inspecting live queue claims. A tested, source-frozen `--analysis-only`
  guardian now regenerates derived receipts without inspecting or launching
  native workers. It preserves failed derived artifacts and native hashes.
  Session: `relay_us_analysis_only_20261126` on Azalea; receipts under
  `full_path_300K/relay_us_heat_recovery_v2_20261125/postprocessing_analysis_only_20261126/`.
  Existing compact sync and MBAR watcher consumed repaired receipts and
  generated shared 10, 15 and 20 ps checkpoints. No trajectory, Hamiltonian,
  bias, sampling gate or native worker source was changed.
- **US movie delivered:**
  `full_path_300K/vmd_progress/current/view_path.vmd`
  contains 1,140 frames (20/window) spanning four accepted blocks, **20 ps per
  window**, with 19,967 stripped atoms and all 98 QM atoms retained. This is
  window-ordered biased sampling, not a continuous reactive trajectory. The curated
  laptop pull includes the complete portable bundle, excluding job/source
  internals. Coordinates, mapping, checksums and Tcl policy were checked;
  graphical inspection remains on the laptop. The old 5 ps path redirects to
  `current`; original 5 ps data are preserved in central clutter. The all-blue
  display was traced to VMD Type coloring by the first character of the common
  `NYL12_` prefix. Policy `20261126-typeprefix3` uses unique leading color
  buckets; 21 other current viewers and their display watcher were refreshed.
- **20 ps sampling decision:** extend the current windows through the already
  approved 30 ps first. The four sub-threshold neighboring 3CV histogram pairs
  are activation windows 13–14, 14–15, 21–22 and 22–23 (intersections 0.091,
  0.081, 0.167, 0.167). Corresponding directed MBAR edge minima remain positive
  (approximately 0.16–0.23); there is no missing connectivity link. Many early
  windows still have only 13–28 decorrelated frames. Both endpoint minima need
  better bracketing. Reassess 25/30 ps block stability and physical-CV support;
  selectively extend weak windows, and consider endpoint flanks or interior
  bridges only where the distributions justify them. The histogram grids differ
  between checkpoints, so their numerical intersections are not a matched-grid
  time-trend test. Do not relax the existing gates or infer a missing window
  solely from a sparse 3CV histogram. Readable 20 ps PMF and overlap plots are
  under `full_path_300K/plots/` with plotting source and source hashes.
- **Plateau releases completed:** ten new bias-free optimizations (images
  11–19 and 21; image 20 reuses window 033) converged to addition-like structures.
  Images 11–15 have C–OG1 ≈1.81 Å and images 16–21 ≈1.74 Å; C–N remains intact
  and HG1 resides on Nalpha. Their maximum gradients are below 9.1×10⁻⁶ Eh/bohr.
  Neither group is certified as a local minimum. Window033's full Hessian has
  unresolved soft-mode signs, and its subsequent repair did not pass the
  finite-difference step-stability gate. Do not repeat a full Hessian or delay
  the current exploratory PMF pending a discrete TI1 claim.
- **Other completed work:** hydration reconnaissance reached 51 ps with acyl
  connectivity retained, but none of 1,021 frames contained a near-attack
  water (closest O–C 4.939 Å). Review water access and the report's unverified
  Nalpha protonation annotation before extending that branch. The WT/mask03/
  mask07 MD replicas and technical Boltz-2/OpenMM control are complete; matched
  MD analyses already exist. The forward relay Davidson calculation found a
  robust negative direction but did not resolve the full saddle index; one
  downhill recovery and both reverse repairs failed their scientific gates.
- **Mechanism-agent:** the read-only RAG pilot completed 20 requests and passed
  citation/provenance checks; factual review remains pending. Its result does
  not qualify scientific planning or job submission. The 20-case model
  benchmark is already complete; the earlier Nano startup failure is a
  superseded attempt, not a reason to restart it. No inference server was
  found on Marigold; Lavender's live check was blocked by SSH host-key
  verification, which was not bypassed.
- **Immediate priorities:** complete derived receipt recovery and inspect the
  updated PMF/overlap/chemistry before changing windows; provide a stripped US
  movie; finish the matched Mechanism-A/B common-coordinate analysis. Review
  completed hydration and design-control outputs before launching more science.
  Force-converged candidate geometries remain distinct from verified minima
  and first-order saddles.

Native release records: `12_nterminal_relay/03_plateau_releases/`.
Current task inventory: `provenance/current_task_status_20261126.json`.
The dated sampling-support audit and historical checkpoints below remain valid
for the data horizons they explicitly describe.

## Current checkpoint — 25 November 2026, 22:49 EST

- **Sampling:** continue the existing 57-window relay campaign through at least 30 ps/window. Automatic MBAR repeats at shared 5 ps boundaries. The first 5 ps overlap graph is connected, but local overlap, effective sampling and chemical/stationary-point validation remain unresolved.
- **Gap audit:** both visible gaps have full-trajectory visits. At s=0.1125, 21 frames become 3 after equilibration removal and 0 after decorrelation; at s=0.6875, 17 decorrelated samples were incorrectly assigned zero mass by weighted-histogram cancellation. Stable per-bin accumulation removes that numerical zero; the highest sampled bin is 20.61 kcal/mol, with **no accepted global barrier**. Original reports are archived; native inputs/sampling are unchanged.
- **Next sampling decision:** inspect 10 ps and 15 ps support in the actual CVs and the projection. Extend existing windows first; add targeted bridges only for persistent weak independent coverage/overlap after checking channel/projection behavior. Keep the fixed grid and predeclared gates.
- **Model policy:** retain GPT/PI scientific leadership; Nemotron remains a read-only evidence/citation assistant. We have no matched GPT-6.1 Sol Ultra versus Nemotron benchmark. The live RAG pilot and factual-count audit continue independently of science.

See [the sampling-support audit](RELAY_US_SAMPLING_SUPPORT_20261125.md) for exact counts, numerical correction, historical comparison, overlap figure and evidence paths.

## Human-guided mechanism-agent handoff — 25 November 2026

The PI requests mechanism-agent to continue the study after November 30 as a
**human-in-the-loop** system. Its governing staged plan is now
`/data/jpc/mechanism-agent/docs/ROADMAP.md`: finish/grade Nemotron vs Mistral;
qualify a better model for a small read-only retrieval pilot; then validate PI
decisions, persistent state and typed scientific adapters before enabling
bounded execution. Scope includes literature, Nyl12/PA6 mechanism work,
transition-region/verified-TS templates, symmetry-aware design and validation,
PI-reviewed experimental shortlists, and evidence-linked manuscript drafts.
Access the existing repositories through an artifact index and analysis tools,
preserving negative results, provenance and scientific limitations; do not
bulk-copy trajectories or treat every installed program as an enabled tool.
Existing US/MBAR science and the technical redesign control remain independent
of model qualification. The November 30 target is a tested handoff, not a claim
that autonomous science is already implemented or a deadline for a desired barrier.

## Operational corrections and authorized control — 25 November, 19:24 EST

- **Relay-US timing correction:** native inputs use `dt=0.0005 ps` and
  `nstlim=10000`: a 5 ps block requires 10,000 steps. Earlier status reports
  incorrectly used 1 fs to convert step counts to ps. At 19:20 EST the
  slowest window's native `mdinfo` reported 2,700/10,000 steps (1.35 ps),
  with about 2 h 55 min remaining at its observed rate. The corrected first
  all-window boundary is approximately **22:15 EST**, followed by receipt
  acceptance, synchronization and automatic MBAR. The earlier 20:15 forecast
  is superseded. The 65-minute production-start stagger and similar SCC/rates
  remain valid. No MD inputs, timesteps or running jobs were changed.
- **Automatic PMF:** all 57 workers and both sync/analysis services were
  verified live. The watcher requires all 57 accepted production receipts
  through a shared 5 ps horizon, then repeats at 10, 15, ... ps. A preliminary
  plot remains distinct from converged sampling.
- **GPU 0 technical redesign control is running:** WT and the accepted
  LigandMPNN sequence S50N/V111T/A187S/V188T receive two matched-seed Boltz-2
  predictions, followed automatically by existing parameterized Amber
  preparation, OpenMM CUDA minimization and raw/relaxed four-site geometry.
  No new sequence sampling, TS claims, activity claims or MD. GPU 1 and US
  remain untouched. Use `design/nyl12_ligandmpnn_pilot/pipeline_control/` and
  its `current_execution.json`. The initial waiting adapter used the wrong
  launcher-status field; it was interrupted before scientific postprocessing,
  retained, corrected to native `exit_status`, tested and requeued.
- **Nano benchmark:** the previous Mistral startup failed during FlashInfer
  autotuning, before any test cases; the controller then stopped Nemotron.
  No confirmed OOM or chemistry comparison exists. One bounded retry disables
  that autotuning option for both models, preserving weights and test cases.
  The retry controller `mechanism_nano_autotune_off_20261125` and both owned
  Nano server jobs were launched at approximately 19:31 EST. No completed cases
  or measured inference ETA yet; plan for overnight and grade chemistry before
  selecting a read-only assistant. The mechanism-agent receipt records live state.

Evidence: `provenance/relay_us_timing_correction_20261125.json` retains the
correction to the earlier timing receipt; the original observation artifacts
remain on disk. Control checks: 14 focused tests, 90% combined coverage, each
new module at least 81%; Ruff passes. This timestamp supersedes live-state
claims in the earlier checkpoint below.

## Verified execution checkpoint — 25 November, 18:56 EST

**US timing audit, 19:04 EST:** the production starts span 64 min 51 s
(17:09:56–18:14:47 EST), explaining almost all of the 2,048–4,848-step progress
spread (1.024–2.424 ps with the verified 0.5 fs timestep).
Current rates differ by only ~7% (1.379–1.477 s/step), with 20–21 SCC
iterations and no unresolved warning markers. Azalea has 28 physical cores /
56 logical CPUs; 57 workers share these resources, but the observed sampling
gap reflects staggered starts rather than a persistently slow worker. The
initial 20:12–20:17 EST forecast used the wrong timestep and is superseded by
the correction above. Receipt validation, synchronization and MBAR must finish
after the true boundary before a preliminary plot exists. Evidence:
`provenance/relay_us_timing_audit_20261125.json` and its correction receipt.

- **Azalea relay US:** 57/57 windows have accepted the 5 ps equilibration and
  are running their first 5 ps production block. Native production progress is
  1,704–4,496 steps/window (0.852–2.248 ps with the corrected conversion);
  zero accepted production blocks and zero reported live
  worker failures at this observation. Heating is no longer the active stage.
- **First new PMF:** not available yet. The cumulative MBAR watcher is active
  and waits for all 57 accepted first production blocks, then updates after
  each shared 5 ps increment. MD process health is not sampling convergence.
- **Violet:** one CPU Sander hydration run (`extend_06`, native time 34.9 ps)
  continues. GPU 1 runs the separate polyamide-amber-ff WT RC1/RC2 parameter
  comparison: the first RC1 replica is in 20 ns production, with six 20 ns runs
  total queued. GPU 0 is idle. US compact-sync/MBAR, VMD, MD-analysis and Nano
  controller processes are monitors, not additional simulation workers.
- **Marigold/lavender:** both verified checkpoints are loaded; first-use CUDA
  kernel compilation continues. No model inference requests or completed
  benchmark cases were observed yet. The 32 passing software tests do not
  establish Nemotron/Mistral scientific performance. Controller retry01 is
  active, with completion/error email and owned-service shutdown scripted.
- **Next decisions:** inspect the first complete relay PMF with overlap,
  chemistry and uncertainty diagnostics; grade paired local-model responses
  only after actual inference. Keep the matched Mechanism-A/B comparison and
  30 November technical design pilot as the existing priorities. No new
  scientific campaign or parameter change was launched by this audit.

Evidence: `provenance/four_host_status_20261125_evening.json`, the US
`current_execution.json`/live runtime summary, and mechanism-agent's
`benchmarks/nemotron_mistral_nano_20261125/controller_retry01_job/`.

## Matched mechanism comparison and redesign target: 25–30 November 2026

The PI requests a defensible comparison of two proposed acylation mechanisms
and a working LigandMPNN → Boltz-2 → OpenMM redesign pipeline by **30 November
2026**. The technical pilot is a deadline target; converged free energies,
verified first-order saddles and increased catalytic efficiency remain outcomes
to establish. Current QM/MM, MD and US jobs continue unchanged. This roadmap
update does not launch or silently extend any campaign.

**Mechanism A — direct proton transfer:** Thr233 OG1 attacks the scissile
carbonyl, with transfer of its hydroxyl proton directly to the substrate amide
nitrogen and C–N cleavage to form the acylenzyme. **Mechanism B — N-terminal
proton relay:** the Thr233 alpha-amino group first accepts the hydroxyl proton
and subsequently supplies a proton to the leaving nitrogen during acylation.
Neither description requires a stable, discrete tetrahedral intermediate.

### Common protocol and reference state

- Match the substrate length/register and direction, enzyme assembly,
  proton inventory/microstate, solvent/ions, atom mapping, force-field release,
  QM atom identities/link boundaries/charge, DFTB3 parameters/electrostatics,
  physical and electronic temperatures, mobile atoms and fixed exterior.
  Identical atom counts alone do not establish a match. Compare the actual
  frozen-coordinate and input hashes before reusing a reference.
- Use one equilibrated E–S ensemble and common basin definitions. If a
  relay-ready conformation is selected, quantify its preorganization cost or
  label its barrier conditional; do not give it a separately favorable zero.
  Prefer a shared accepted E–S endpoint and frozen environment for the bounded
  comparison, then assess sensitivity to additional independent starts.
- Apply the same heating/equilibration, timestep, thermostat, bias-potential
  convention, output cadence, SCC/force validity rules and production-block
  analysis. Start with 30 ps/window in retained 5 ps blocks, and extend based
  on the same convergence gates. Equal window counts or equal elapsed times
  are not required: coverage and statistical precision are the matching goals.
- Record a common set of acylation/proton-transfer coordinates and all relevant
  component distances/proton coordination in both campaigns. Window centers
  can follow different channels; their exact native biases must be evaluated
  and removed consistently. Keep both elementary reactions asynchronous.
- Reconstruct a common multidimensional free-energy description where sampled
  support permits, then compare channels using the same coordinate measure,
  binning/metric, normalization and E–S basin reference. A marginal in
  `q_acyl` and a density along normalized path arclength are different
  observables; their maxima cannot simply be subtracted. Separately zeroed
  path PMFs also cannot establish a route preference without compatible
  reactant/channel normalization and any necessary restraint corrections.
- Monitor all exchangeable protons, not only HG1. If direct-seeded windows
  follow the N-terminal relay, report channel interconversion rather than
  calling their marginal a pure Mechanism-A PMF. Do not introduce artificial
  proton-identity walls merely to preserve a pathway label. With a common
  Hamiltonian and fully equilibrated common observables, seed-dependent global
  PMFs should agree; mechanism differences must be resolved as physical
  channels and bottlenecks, not as an unexplained disagreement between runs.
- Require connected overlap, correlation-aware effective sample sizes,
  block/independent-start stability, stationary orthogonal observables,
  chemical retention and bootstrap uncertainty. Use MBAR consistently, with
  WHAM as an estimator cross-check when the chosen coordinates support it.
  Report uncertainty in the route comparison, including reference uncertainty;
  sparse matched DFT/MM checks remain necessary for fast-QM interpretation.

The existing Mechanism-A-seeded campaign is valuable evidence, but its
approximately 47.05 kcal/mol `q_acyl` marginal is nonconverged and mixes direct
and relay configurations (low target ESS and drift flags). It is not a
matched pure-direct activation free energy. First reanalyze its retained data
in the common observables and audit physical identities. Reuse valid evidence;
prepare a matched direct-channel control only for the missing support. Keep
the present relay campaign running. No high barrier is discarded or lowered
to meet the deadline.

**Compatibility audit, 25 November:** retained direct and relay records have
identical topology, QM/mobile selections, 3ob parameter trees, native QM/MM
electrostatics and 300 K physical/0 K electronic temperatures. Five direct
seeds compared with the relay source have exactly matching cell parameters;
100,757 frozen atoms differ by at most 3.81×10⁻⁶ Å per Cartesian component
(RMS 1.56×10⁻⁶ Å), consistent with restart serialization precision. No fit or
periodic realignment was used. The direct calculation used `scfconv=1e-10`
versus `1e-12` for the relay. Check this numerical/exterior precision effect
on representative common configurations before any pooling. The next bounded
step is common-observable reanalysis and support assessment, not an automatic
full direct-route rerun. Matching physical inputs does not fix mixed channels,
reactant normalization, low ESS or lack of convergence.
Source-bound comparison evidence is retained in
`qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/direct_relay_compatibility_20261125.json`
and included by the existing laptop pull script.

### Work sequence and measurable deliverables

| Target dates | Priority | Deliverable and acceptance |
| --- | --- | --- |
| 25–26 Nov | Relay sampling and comparison definition | Retain automatic 5 ps MBAR checkpoints; audit overlap/chemistry as soon as all windows have production. Freeze the direct/relay match table, common reference and channel observables. Prepare a tested continuation before requesting sampling beyond the current six-block MD limit. |
| 26–27 Nov | Direct-channel comparison | Reanalyze existing direct-seeded frames first. Prepare missing matched seeds/inputs and resource plan using Violet/Azalea headroom. Launch a control only after topology, frozen exterior, CV/bias parity and restart checks pass. |
| 27–29 Nov | Mechanistic qualification | Compare common-reference free-energy channels with uncertainty; resolve coverage/hysteresis failures. In parallel, finish bounded TS-curvature/downhill tests and select representative barrier-region structures without claiming an unverified structure is a true TS. |
| 25–29 Nov, parallel | Design-pipeline commissioning | Verify reactive-atom mapping, protected-site manifest and separate alpha/beta sequence ties; run a small existing-structure/WT control through LigandMPNN, Boltz-2, OpenMM and geometric reporting. This is an engineering test, not evidence for a new catalytic design. |
| By 30 Nov | End-to-end pilot and decision | Deliver reproducible scripts, current slide deck, matched-PMF status table and a complete technical pilot. Start a bounded scientific design campaign when the WT PMF/template gates below pass; otherwise deliver a runnable pipeline with explicit scientific holds and the next resolving calculation. |

Use the mechanism with supported chemistry and a qualified WT free-energy
profile for the first scientific design pilot. Prefer a verified PES saddle
as the chemical anchor/control and three to five representative near-barrier
conformations for ensemble testing. A PMF maximum is not automatically a
dynamical TS; retain that distinction and use appropriate connectivity or
committor diagnostics. If only a qualified transition-region ensemble is
available, the previously approved provisional ensemble pilot may proceed
with that label; do not invent a verified saddle to meet 30 November.

The design pilot preserves catalytic/maturation context and assembly symmetry
using the reviewed protected-site/alpha-beta tie manifests. Boltz-2 predictions
and OpenMM relaxation test structural plausibility, not activation barriers.
Use validated parameters for noncovalent E–S relaxation. Do not relax a TS or
covalent intermediate with the E–S bond topology/charges: retain explicitly
declared motif constraints for structural tests or use QM/MM refinement.
Evaluate whole sequences across E–S and reactive-state templates, with WT and
known-active controls; native-Thr similarity remains a diagnostic. Reserve
matched mutant/WT barrier calculations and experiments for a small shortlist.
No new force-field-development, large substrate-modeling or assay-temperature
campaign is added to this deadline's critical path.

### Presentation narrative

Explain experiment and ultimate catalytic-efficiency goal first, then the
parameterization/MD/protonation/QM-MM workflow and the two proposed mechanisms.
Define geometric criteria before showing structural classifications. Refer to
both substrate directions, and show the position of the scissile amide N
relative to Thr233's amino group rather than relying on the word “flipped”.
Use Mechanism A/B consistently; include pending-result placeholders and keep
PES maxima, preliminary PMFs and verified TS barriers distinct. The current
presentation filename and visible date must reflect the actual build date.

The PI prioritizes a converged free-energy profile for the current modeled
E–S → N-terminal relay/addition-like region → acylenzyme reaction. A discrete
TI1 minimum and full PES TS certification are not prerequisites for sampling.
Those diagnostics continue in parallel; the umbrella Hamiltonian, mapped
chemistry and sampling validation are the critical path.

**Exploratory-template clarification, 2026-11-26:** full-path PMF convergence
and a verified PES saddle are not prerequisites for a small provisional
LigandMPNN template pilot. Distinct chemically plausible US barrier-region
frames can be used directly with explicit provenance and provisional labels.
Quantitative barrier/mechanistic claims and larger design commitments retain
the qualification gates below. A verified saddle remains a desirable control,
not a requirement to minimize every thermal parent. See the design roadmap
clarification for fixed side-chain context and hydrogen/charge feature limits.
No new design campaign was launched.

**Redesign decision, 2026-11-25:** obtain a stable, converged and chemically
credible WT PMF before the next mechanism-guided sequence campaign. A reasonable
barrier is a result to evaluate, not a numerical target to force. Afterward,
compare one verified PES saddle with three to five actual representative
finite-temperature conformations near each relevant PMF bottleneck. Keep the
saddle as a chemical anchor and single-template control; retain ensemble
parents without minimizing away their diversity. Near a PMF maximum is not
automatically a dynamical TS: use bounded bias-free shooting/committor checks
where that claim matters. A stable TI1 is optional. The proposed matched-budget
pilot, explicit alpha/beta sequence ties, cross-state evaluation and primary
precedents are specified in `design_ranking_and_iteration_20261022.md`, section
“PES saddles and transition-region ensembles”. No new design campaign is launched
by this decision; US/MBAR remains the compute priority.
If saddle certification remains unresolved after the PMF/ensemble is qualified,
defer the single-saddle control rather than block a provisional ensemble pilot.

**Current US recovery checkpoint, November 25, 13:26 EST:** all 24 initial
Sander heating jobs completed normally (2,000 steps each, exit 0, resolved SCC),
but the workflow wrongly required an equilibrated temperature after the 1 ps
25→300 K ramp. Their final temperatures were 172.79–182.21 K. All 48,000 steps
and checkpoints were preserved and independently adopted as heating data;
none counts as production. Heating now has an explicit ramp-only acceptance
role; the existing 5 ps constant-300 K equilibration follows from those exact
checkpoints, with the stricter equilibrium/production gates retained.

The recovered campaign now has **57 actual Sander workers / 57 windows / zero
worker failures** at the stated observation. Azalea has 28 physical cores and
56 logical CPUs. All windows can progress concurrently with slight CPU time
sharing; no window was removed just to match hardware. A source-frozen staged
28→40→57 worker benchmark and reversible all-logical/physical-core affinity
comparison measure aggregate throughput while preserving MD/RNG/checkpoints.
Sampling convergence remains unestablished. The stable status files are
`12_nterminal_relay/umbrella_sampling/full_path_300K/current_execution.json`
and `runtime_summary.json`; these supersede the earlier failed-run summary
and are included by the curated pull. Exact recovery plan, retained failures,
acceptances and resource results are under `full_path_300K/relay_us_heat_recovery_v2_20261125/`
and `full_path_300K/recovery_operations_v2_20261125/`.

**Earlier launch checkpoint, November 25, 12:38 EST:** 24 one-thread Sander
workers were verified running on Azalea, with zero worker failures at that
observation. The full schedule has **57 states: 31 activation + 26 delivery**;
these include bridges needed to keep each biased-axis increment below 0.15 Å.
Inputs are frozen in `12_nterminal_relay/umbrella_sampling/full_path_300K/`.
The authoritative process snapshot is
`full_path_300K/relay_us_azalea_20261125/runtime_summary.json`.
All windows receive 1 ps heating, 5 ps equilibration, and their first 5 ps
production before any production extension. Six retained production blocks
give the initial 30 ps/window target; duration alone does not establish convergence.
Both-chart native compatibility and end-to-end MD checks passed. The real
native-data → independent DUMPAVE-CV parity → compact-copy → strict-loader
chain also passed for both charts using the actual analysis request. Native
restraint CVs and trajectory times are aligned using the stage's restart time;
rounding bounds are explicit. The terminal frame lacks a native DUMPAVE match
and is retained with parity marked unevaluated, rather than censored or fabricated.
All failed derived attempts are preserved. A source-frozen 120-hour Azalea
guardian revalidates even previously completed derived records, repairs only
derived outputs, and resumes clean analysis-held workers without repeating
accepted MD. A 120-hour compact sync runs every 60 seconds; raw trajectories
remain on the execution host. Extensions wait if diagnostic acceptance fails.
The replacement seven-day analysis watcher (`nyl12_relay_us_analysis_20261125_c03`)
waits for matched accepted 5 ps production increments across all 57 states,
then retains each cumulative MBAR report and a common-reference progression
plot under `full_path_300K/analysis/cumulative_progression_20261125/`. Disconnected
overlap cannot produce a global PMF. No accepted production receipts existed
at replacement, so no new PMF exists yet. The watcher can analyze horizons
beyond 30 ps, but the current MD runner stops at six blocks; a tested native
continuation must be prepared before extending that allocation. See the stable
`12_nterminal_relay/umbrella_sampling/README.md` for the exact scope and gates.

1. Finish the generalized native CV/analysis adapter and runnable inputs.
   Activation windows bias `(q_acyl, q_act)`; delivery windows bias
   `(q_acyl, q_delivery)`. All three features must be evaluated for every
   sampled frame and all cross-window biases. Both chart families share one
   physical Hamiltonian and overlapping relay-state sampling.
2. Cover the **whole accepted path**, including E–S and acylenzyme, before
   extending already sampled windows. Place centers no more than 0.15 Å apart
   in either chart axis, adding bridge windows at the coarse proton-transfer
   gap. The earlier 18-window local-grid preparation is a protocol fixture,
   not a reason to delay full-path coverage.
3. Run short native compatibility/bias-parity checks, one per chart, before
   dispatch. Prefer idle Azalea CPU capacity; add Violet workers only within
   measured headroom, preserving current jobs. Use native host-specific Amber
   binaries/environments, one thread per window, source/input hashes and tmux.
4. First accumulate an accepted 5 ps production block in every window.
   Run overlap/MBAR diagnostics once coverage exists; clearly mark gaps or
   disconnected components. Then extend in retained 5 ps blocks toward at
   least 30 ps/window and further as required by the eight existing convergence
   gates. No nominal duration alone establishes convergence.
5. Evaluate proton inventory/coordination, heavy-atom chemistry and orthogonal
   environment distributions throughout. Use a common statistical normalization
   across both chart families; independently zeroed PMFs must not be added.
   A final barrier also needs a defensible bound-state/conformational reference,
   method sensitivity and appropriate uncertainty. A desired numerical barrier
   is not an acceptance criterion.

This campaign stays at **300 K** with the existing 104,922-atom single-substrate
system, 98 real QM atoms/eight links/charge −2, original 4,165-atom 12 Å mobile
region/frozen exterior and zero electronic temperature. The separately prepared
348.15 K full-system WT branch is not substituted silently. The present goal
is acylation in the modeled bound complex; deacylation/product release remain
separate mechanistic work.

Current preparation/execution records belong under
`qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/`.
The source-bound combined movie is `12_nterminal_relay/vmd_progress/view_path.vmd`.
Detailed convergence and interpretation requirements remain in
`docs/qmmm_umbrella_sampling_and_pmf_protocol_20261109.md`.

## Supporting checkpoint — retained endpoints and TS work, 2026-11-25

**Combined inspection movie:** `qmmm/pa6_flipped_model30/12_nterminal_relay/`
`vmd_progress/view_path.vmd` contains 39 frames / 19,967 solute atoms, E–S→shared
relay/TI1-like seed once→acylenzyme. All frame coordinates and periodic cells and
15 bundle checksums passed, including an independent copied bundle. The pull
script includes the root bundle; original 39-frame plot numbering is unchanged.
Reverse windows 012→013 contains a coarsely sampled HG transfer/C–OG change, with smooth
local torsions/RMSD and no branch-switch evidence. A bounded bidirectional
4-window refinement per direction is being prepared under the same branch.

**Current follow-up:** plotted images 11–21 are being released without reaction
bias under the same 12 Å mobile-region protocol. Image20 is the identical,
already force-converged window 033 seed and is reused; ten other images are
running concurrently in one-CPU tmux sessions, up to 10,000 cycles each, with
independent force/SCC/exterior checks. Inputs and plot-index mappings are in
`12_nterminal_relay/03_plateau_releases/prepared.json`. Compare chemical
retention and distinct endpoints before further curvature calculations; these
neighboring path starts do not estimate equilibrium basin populations.

The forward candidate's conditional real-QM Hessian has one negative mode
(−0.06495 Eh/(Å² amu)); its largest displacement is transferring HG1. Both
finite-difference directional curvatures are negative (−0.10515/−0.10541
Eh/Å²). The inherited relative precision gate of 1e−4 failed (observed 0.00238),
and that failed receipt is preserved. A separately documented negative-sign
gate permits exploratory downhill connectivity only when both estimates are
negative and differ by<1% of their smaller magnitude, followed by actual
energy and outward-slope decrease. This does not certify full-active index one
or relax the original precision qualification for that purpose. Both-sign
all-mobile downhill jobs passed actual energy and outward-slope checks and are
minimizing; a bounded 80-HVP full-active Davidson diagnostic is also running.
Both guarded reverse dimer repairs stopped safely (attack: energy-rise limit;
PT: SCC guard). Two original-seed partial-QM screens are authorized as
nonstationary search-direction diagnostics before further reverse attempts.
Final stationary-force/TS verification gates are unchanged.

**SCC guard correction:** an intermediate “not converged after70cycles” message
may be followed by valid convergence (the retained window 033 release includes
70→96-cycle recoveries). The first guard incorrectly rejected these recoverable
warnings. Old aborted attempts are preserved. A new frozen policy allows
bounded recovery through500actualSCCiterations, still rejects fatal or
unresolved returned electronic evaluations, and retains the same zero electronic
temperature, SCF/density tolerances, energy-rise limit and finalforce criterion.
The minus downhill retry starts from the original independently validated
0.05 Å displacement, not an interrupted electronic trial. This is warning-
classification repair, not a change in the Hamiltonian or physical convergence.

**TI1 terminology:** window 033 is a retained, force-converged, partially
addition-like candidate (C–OG1=1.740 Å), not a verified stable TI1 minimum.
Its source was a biased QM/MM MD frame, followed by bias-free optimization;
the matched window 032 seed lost addition character under the same protocol.
Conformational/environmental differences are implicated, but no individual
stabilizing interaction or force-field change has been shown causal. See
`docs/neb_us_nterminal_relay_release_20261123.md`, section “Why window 033 retained
addition character”, and the updated Results/Methods.

Both exploratory branches reached geometrically appropriate endpoints and
retained them after removing the reaction bias, within the unchanged 12 Å
movable model/frozen exterior. Independent maximum gradients are 4.79e−6
Eh/bohr (E–S) and 9.88e−6 Eh/bohr (acylenzyme). E–S C–OG1/C–Nleave are
2.963/1.380 Å; acylenzyme values are 1.327/3.019 Å. HG1 is on OG1 in E–S
and on the leaving nitrogen in acylenzyme; Nalpha retains two hydrogens in
both endpoints. Atom/proton inventory, topology and exterior are shared.

The combined driven PES maximum is 16.31 kcal/mol above the regenerated E–S
reference; the forward dimer candidate gives 16.312 kcal/mol from independent
same-Hamiltonian energies. Acylation endpoint ΔE is −4.15 kcal/mol. These are
**potential energies**, not a PMF or an established biochemical barrier. The
forward candidate meets its dimer force tolerance (maximum 3.40e−5 Eh/bohr)
but still needs curvature/index and downhill-connectivity verification. The
reverse TS attempt was stopped after thousands of SCC failures/resets, large
energy/curvature excursions and return to E–S-like geometry; its last checkpoint
and diagnostics are retained. This does not invalidate the accepted E–S endpoint. Window033 minimum status remains
unresolved. No new full-shell Hessian or US campaign has been launched.

Updated common-reference plot/table/source hashes:
`qmmm/pa6_flipped_model30/12_nterminal_relay/plots/connected_relay_pes.*`.
Current branch movies remain `BRANCH/vmd_progress/view_path.vmd`.

Hydration completed all accepted blocks through 26 ps (521 saved frames),
with no annotated near-attack water; closest remains 5.422 Å at 9.7 ps.
The conditional extension stopped before new native sampling because its
analysis manifest hashed a lifecycle status file that changed during final
postprocessing. Repair is limited to stable status/manifest handling; native
data are accepted and retained. Do not claim that the 51 ps continuation ran.

Four WT/mask_07 local inspection structures and corrected PML scripts are
available under `md_followup/analysis/structural_followup_20261124/`. Coordinate
and atom-map readback passed; PyMOL was unavailable for a rendered GUI check.
Further extraction of all eight existing CSVs found no sampled Tier A/B
observation with Nalpha–Nleave ≤3.5 Å; conditional medians are 5.99–6.34 Å.
This is a distance diagnostic, not proof against an induced-fit or water-
assisted relay. It strengthens the need to establish the reactant-ensemble/
preorganization reference before treating a low conditional PES/PMF barrier as
an overall enzyme barrier. Native-Thr similarity remains diagnostic only.

## Previous checkpoint — approved relay exploration, 2026-11-24

This dated block supersedes historical dispatch gates below. The PI explicitly
approved exploring window 033 without full-minimum verification and preparing a
348.15 K WT branch. Window033 is force-converged (maximum gradient
7.94e-6 Eh/bohr), with tetrahedral-like geometry, but is **not a verified local
minimum or chemically accepted TI1**. The reused QM-block Hessian audit did not
resolve all signs under its uncertainty criteria; no more full Hessian work is
queued. Original failed minimum-verification receipts remain intact.

1. **Primary: N-terminal-relay acylation.** Two reviewed exploratory branches
   under `qmmm/pa6_flipped_model30/12_nterminal_relay/` connect window 033 toward
   E–S (`01_activation_reverse`) and toward acylenzyme through N-terminal proton
   delivery (`02_acylation_forward`). Explicit exploratory dispatch records
   supersede the old minimum gate only for these branches. The seed has
   Nalpha–HG1=1.050 Å, HG1–Nleave=1.975 Å and a 179.4° donor–H–acceptor angle.
   Use the identical DFTB3 Hamiltonian, 98 QM atoms/charge −2, 12 Å movable
   region and fixed exterior. Each asynchronous driving window is followed by
   chemical/force checks and updated movie/energy/distance plots. Candidate
   endpoints receive unbiased optimization; eligible interior maxima may seed
   dimer searches. Neither branch completion nor force convergence certifies
   a TS. Actual running/failed status is in each branch's retained job records.
   Both branches were dispatched at 22:29 EST on November 24 in tmux
   `nyl12_w033_exploratory_reverse_20261124` and
   `nyl12_w033_exploratory_forward_20261124`, one nice CPU each. Inspect
   each branch's `exploratory_job/`, `branch_status.json`, `plots/` and
   `vmd_progress/view_path.vmd`; derived products update after completed windows.
2. **Secondary: bounded deacylation hydration.** The existing 300 K run continues
   toward 26 ps. A detached conditional queue waits for all accepted blocks and
   a complete water analysis, then adds five 5 ps blocks, to **51 ps total**, only
   if no near-attack water is found. Any annotated candidate holds the queue for
   review; failures stop explicitly. No classical acylenzyme parameters or water
   insertion are introduced. Current accepted analysis: 221 frames through
   11 ps, closest Ow–acyl-C 5.422 Å at 9.7 ps, angle 47.59°; no combined
   ≤3.5 Å/95–125° annotation. The current native segment extends beyond this
   accepted analysis. Movie: `13_deacylation_pilot/hydration_reconnaissance/`
   `vmd_progress/view_hydration.vmd`. Native MD acceptance and repaired export
   receipts remain separate; original failed refresh attempts are retained.
3. **Matched MD finished.** All eight 100 ns WT/mask_03/mask_07/V111C runs passed
   completion and final-temperature gates; both independent replicas are
   analyzed at four sites. Mask_07 Tier A is higher than WT in both replicas,
   but its A-or-B advantage is not reproduced in replica 2. Known-active V111C
   scores lower than WT: these are geometry-retention diagnostics, not activity
   predictions. Native-Thr similarity remains diagnostic only. Final table and
   plot: `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/`
   `analysis/matched_replica_comparison.{json,csv,md,png}`.
4. **Temperature: prepared, not launched.** Two independent full-system WT E–S
   endpoints are staged under the same `md_followup/temperature_348K/` directory:
   200 ps 300→348.15 K NVT, 1 ns NPT and 2 ns unrestrained NPT, then review.
   Forty-two staged artifacts passed hash checks. Keep ongoing comparisons and
   hydration at 300 K. The four-substrate design-MD topology is distinct from
   the single-substrate legacy QM/MM system; map and equilibrate a selected
   mechanistic state explicitly before any assay-temperature PMF. Do not heat
   only a local frozen-shell model and claim full-system thermal equilibration.

**Current PI follow-up (November 24):** further WT/mask_07 structural
evaluation and an automatic accepted-window connectivity audit are being
prepared. Both relay branches accepted window 001 and continued into window 002.
The required endpoint is an auditable E–S → relay → acylenzyme connection; the
next PMF coordinates are selected from that chemistry, not from a desired
barrier. Detailed preregistered gates and candidate local coordinates are in
`docs/qmmm_umbrella_sampling_and_pmf_protocol_20261109.md`. No new umbrella
campaign is dispatched before the chemistry/mapping/coordinate gates pass.

**Barrier strategy and stop rules.** Seek a chemically credible pathway, not an
artificially low number or an obligatory TI1 minimum. First inspect forward and
reverse relay paths, proton inventory and distance continuity. Use successful
paths to choose proton-transfer and heavy-atom progress coordinates for new
QM/MM umbrella sampling, with connected overlap, effective sampling, block and
independent-start checks. Track all exchangeable N/O hydrogens. Missing slow
proton/environment coordinates can distort a projected PMF (Rosta et al., 2009).
Apply sparse DFT/MM checks to new E–S/TS-region/product configurations; a
correction requires energy-distribution overlap. Account for protonation-state
population costs before comparing different microstates. A verified TS still
requires first-order curvature and downhill connectivity.

At 348.15 K, hypothetical ΔG‡=25 kcal/mol gives k≈1.47e−3 s⁻¹ (κ=1), whereas
28 kcal/mol gives 1.92e−5 s⁻¹. These are elementary bound-step illustrations,
not rate estimates from crude-lysate conversion. See
`analysis/nyl12_mechanism_priority_20261124/eyring_temperature_reference.csv`.
Neither raising temperature nor more sampling guarantees a lower ΔG‡. The old
35.81 kcal/mol **potential-energy** path already received umbrella sampling;
its latest apparent PMF barrier was about 47.05 kcal/mol. Do not repeat that
high-barrier direct route as the main campaign. Legacy ASH274 NEB may continue
as authorized, but remains a secondary unconverged diagnostic.

**Remaining controls.** If no attacking water is found by 51 ps, inspect
retained-product occlusion/solvent access before further extension; any product-
separation control must preserve atom balance. Keep the matched MD parameter
release fixed. The newer polyamide force-field release remains a research
candidate with documented accuracy limitations; consider one bounded WT
sensitivity test before strong quantitative design/PMF claims, rather than
silently replacing the force field. Crude MS assays do not determine kcat or a
unique barrier; L4→L2 remains a modeling hypothesis. Mechanism-agent stays
stopped, with tested email alerts ready for later authorized runs.

## Historical checkpoints — retained evidence, not current dispatch instructions

## LIVE STATUS — completed hydration pilot, 2026-11-24

The retained-product deacylation hydration reconnaissance completed natively:
1 ps heating plus 5 ps 300 K NVT DFTB3/MM, 0.5 fs steps, original 12 Å movable
region and frozen exterior. Final temperature 300.40 K; native
process, SCC history, finite-state and frozen-boundary checks passed. Acyl
C–OG1 = 1.310 Å and C–leaving-N = 3.250 Å at the endpoint.
Nearest existing water O is still 5.786 Å from acyl C (initial 6.246 Å),
so the final structure has no close attacking water. This short MM-water
hydration check is not a deacylation pathway or evidence of hydration convergence.
Retained report: `qmmm/pa6_flipped_model30/13_deacylation_pilot/`
`hydration_reconnaissance/hydration_endpoint_status_20261124.json`.

Next: inspect the hydration trajectory for transient accessible waters before
choosing a QM-water model; reoptimize any expanded-QM starting state before
driving hydrolysis. Do not extend blindly or infer that water cannot enter
from this 6 ps pilot. The full-window 033 minimum remains unqualified after
the cutoff-column repair failed; neither conditional N-terminal relay path
is running. Diagnose chemical-subspace curvature separately from unresolved
full-shell numerical qualification rather than repeat the entire Hessian.
No new pathway is dispatched on an unresolved gate. Four second-replica GPU
MD runs and their postcompletion analysis watchers continue. Legacy ASH274
NEB is still nonconverged (current provisional profile maximum ~46.7 kcal/mol).
Two Gaussian terminal-pair diagnostics on Violet belong to polyamide-amber-ff,
not the Nyl12 deacylation pilot.

Mechanism-agent readiness: the PI confirmed email receipt on November 24.
Completion/error/PI-decision alerts are enabled by default for future run CLI
commands (152 tests passed). Its idle project Ollama server on Azalea and
SSH tunnel on Violet are stopped; installed models, environments and ongoing
scientific jobs are preserved. No review/controller is left running.

## CURRENT — native-Thr diagnostic-only policy and replica-2 analysis, 2026-11-24

The PI approved native-Thr reference similarity as a **diagnostic, not an
eligibility filter or ranking veto**. Original displacements, torsions,
thresholds, native_pass and historical joint_pass flags remain available.
Future eligibility uses original catalytic Tier A/B geometry with unchanged
chemical integrity, protonation, atom mapping and relevant process/force
checks. This does not certify activity or lower a mechanistic/TS criterion.
Native-reference failures alone no longer block MM→QM/MM admission or
candidate extraction; relay candidates retain their chemical ranking.

Replica 1 findings: mean Tier A/B occupancy WT/mask_03/mask_07 is
45.7/23.3/69.5%. The old joint native/Tier A/B fractions 0.05/0.30/0.25%
are descriptive only. Within Tier A/B frames, OG1/Nalpha/psi fail together;
chi1 is usually retained. Static minimized WT OG1 shifts 0.757–0.880 Å
barely exceed 0.75 Å, whereas MD median shifts are larger. Alternative
scaffold fits and unchanged circular torsion definitions were audited.
The native crystal reference is not a validated activity classifier.

The same four-site, fixed-scissile-amide/PBC analysis is queued for all
four second replicas at their execution hosts: Violet focused_mask03 and
V111C; Azalea WT and broader_mask07. Detached sessions verified running:
`nyl12_replica02_analysis_violet_20261124` and
`nyl12_replica02_analysis_azalea_20261124`. All four cases are pending MD
completion at dispatch. Analyze only completed, validated
100 ns production data, at 100 ps stride with all-frame denominators and
first/second-half summaries. Native similarity is reported separately.
One CPU analysis runs at a time per host; MD/QM jobs remain untouched.
Outputs: existing `md_followup/analysis/replica_02/DATASET/`; persistent
queue/process/scientific receipts under `replica_02/queues/HOST/` establish
actual dispatch status. No second-replica outcome is claimed before analysis.
Compare independent velocity replicas explicitly; four sites and successive
frames are nested/correlated, not independent replicate experiments.

## CURRENT — window 033 cutoff diagnosis and localized repair, 2026-11-24

**Matched MD preliminary analysis is now available:** replica 1 WT/mask03/
mask07 were sampled every 100 ps at all four fixed sites (1,000 frames/site).
Mean Tier A/B occupancies are 45.7/23.3/69.5%; joint native-Thr+TierA/B fractions
are 0.05/0.30/0.25%. Broader-mask geometry retention is encouraging, but this
is one trajectory per sequence with starting-pose confounding and almost no
joint gate passes, not evidence of improved catalysis. The four-site adapter
is implemented/tested and these three analyses completed. Remaining
replica-2 analyses are now queued as described above. Compact source-bound table:
`design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/analysis/preliminary_comparison.json`.

**Latest repair outcome (November 24):** the localized repair finished but
54 of 383 columns failed the unchanged adjacent-step stability test. It stopped
before deriving a new full matrix; no minimum was qualified and neither relay
path has started. This is a numerical qualification failure, not evidence of a
negative chemical mode. Original failed receipts and all new columns are retained
in `07_cutoff_repair/`. No further Hessian evaluations are automatically queued.
The hydration pilot and matched MD continue independently. Earlier launch
statements below describe their recorded state at submission.

The full 12495-DOF Hessian finished natively (exit 0). All symmetrized
matrix eigenvalues are positive, lowest 2.89147e-5 Eh/(Å² amu), but the
unchanged uncertainty criterion leaves 2339 signs unresolved. Window033 is
therefore **not formally minimum-verified**. The conditional relay watcher
stopped at the failed scientific gate; neither new relay branch started.
Do not interpret the initial-only movies as completed reaction paths.

A bounded native follow-up supports a concrete cause: finite displacements
crossed the 10 Å nonbonded cutoff. Implicated MM pairs lie at 9.999925 and
10.000087 Å. TYR115 HH xy nonreciprocity fell from .0040027 to 8.67e-11,
and WAT10834 H1 yz from .0036062 to 9.70e-10, with smaller displacements.
All six tested columns pass the unchanged finer-step relative stability gate.
This supports a numerical cutoff artifact; it is not a minimum classification.
See window 033 `minimum_verification/numerical_diagnosis_20261124/`.

One localized repair is authorized: up to 383 implicated columns (3.1% of the
full matrix), at h=3.125e-5 and 1.5625e-5 Å, with at most 1,532 native gradient
calls and eight nice CPU workers. Original failed matrices/receipts remain.
Every replacement must pass the original 1e-4 relative step test, followed by
complete-matrix inertia with the same sign/uncertainty/residual criteria.
No threshold, protonation, Hamiltonian, movable shell or frozen exterior is
changed. No relay dispatch is authorized on an unresolved result.

The localized repair is now actually running in tmux
`nyl12_w033_cutoff_repair_20261124`, separate `07_cutoff_repair/`. Its nine
focused tests passed (93% script coverage); original failure receipts and relay
stop remain unchanged. Hydration is running in tmux
`nyl12_deacylation_hydration_20261124`, under `13_deacylation_pilot/`
`hydration_reconnaissance/`; four focused tests passed (90% combined coverage).
Native SCC/integration steps advance. A two-gradient same-coordinate diagnostic
passed: all 4,165 active analytic gradients are identical at retained precision
with/without native belly pruning. QM energies also match; the MD total-energy
offset is not a chemical reaction energy. Neither calculation is yet scientifically
complete; inspect its receipts rather than inferring success from submission.

Matched MD snapshot on November 24 morning: broader_mask07 replica 01,
focused_mask03 replica 01 and V111C replica 01 completed their 100 ns production
protocols. WT replica 01 was near 99 ns on Violet; V111C replica 02 had begun
production there. Azalea runs WT replica 02 (~26.5 ns at 08:40 ET) and broader
mask07 replica 02 (~29.8 ns). Violet's **retry_01** GPU0 queue is live and will
start focused_mask03 replica 02 after WT; its retained original failure is
superseded and must not be mistaken for a blocked active queue. Progress values
are snapshots, not completion receipts. No new geometry/activity conclusion
is drawn before the four-site streaming analysis is implemented and run.

Deacylation: the retained-product force-accepted acylenzyme from the old direct
route has no water O within 6 Å of its acyl carbon. Nearest WAT23286 O is 6.246 Å
from C and 8.571 Å from neutral Nalpha, so it is not reaction-ready. The approved
bounded follow-up is1 ps heating plus5 ps NVT hydration reconnaissance using
existing solvent, original 98-real-atom QM/charge −2, and the 4,165-atom movable
12 Å region. No water is inserted, no reaction/position bias is applied, and
the complete amine product is retained. Water is MM during this hydration-only
stage; a candidate must be promoted to QM and reoptimized before hydrolysis or
water-mediated proton transfer can be tested. This old endpoint has its own
frozen exterior: do not splice its energies to window 033. Audit and preparation
are under `qmmm/pa6_flipped_model30/13_deacylation_pilot/`. Check native job
receipts before calling the pilot running. The existing GPU MD and legacy
ASH274 NEB remain unchanged; the four-site MD analysis adapter is still pending.

## CURRENT — retained window 033 candidate; full Hessian running, 2026-11-23

### Approved one-week finish line — 2026-11-23

The PI has approximately one more week of hands-on project time. The full
scientific endgame remains experimentally improved Nyl12, supported by an
adequately validated catalytic model. Finishing a calculation, retaining TI1,
or writing a manuscript is not equivalent to that success.

The PI approved the one-week plan and conditional pathway/TS continuation.
For this week the primary question is whether an N-terminal proton
relay can connect the modeled PA6 E–S complex to acylenzyme through a supported
acylation pathway. A discrete TI1 is optional. A publishable mechanistic claim
requires chemically mapped endpoints and proton inventory, qualified minima,
index-one saddles with chemically relevant modes and appropriate downhill
connectivity, and a bounded Hamiltonian/method-sensitivity assessment. Reference
energies must use consistent Hamiltonians, boundary coordinates and constraints;
independently released US frames must not become a spurious barrier by total-
energy subtraction. PES and PMF claims remain distinct. L4-to-two-L2 cleavage
and its relevance to a longer polymer remain model hypotheses.

| Deadline | Bounded deliverable and decision |
|---|---|
| Days 1–2 | Finish window 033 curvature classification. If unstable, follow the resolved negative mode in both directions and classify the endpoints; do not force a TI1 well or repeat full matrices without a demonstrated numerical reason. Finish the existing eight MD trajectories when their queues complete. |
| Days 2–4 | Pursue one N-terminal-relay acylation branch with compatible endpoints and explicit proton tracking. Refine/verify its saddle(s) using the existing fast Hamiltonian; use already available DFT benchmarks and at most a bounded set of representative new checks. By the end of day 3 require a connected chemically interpretable path and a credible stationary-point verification route to justify further mechanism computation. |
| Days 4–5 | Analyze all four occupied sites in both independent velocity replicas per sequence. Report block/replicate variation, initial-pose bias, native-Thr failures and retention denominators. Select zero to a few experimental candidates, not a forced winner. |
| Days 6–7 | Freeze claims, figures, native inputs/receipts and manuscript/SI; perform a reproducibility and uncertainty audit. Explicitly label unverified or incomplete calculations. A review-ready package is achievable; publication or biological improvement is not guaranteed. |

No new broad model generation, protonation survey, LLM deployment, force-field
development or long umbrella campaign is included this week. Additional
US is justified only after a connected chemically valid path exists and sampling
can be assessed by overlap, effective sample size, block stability and relevant
orthogonal chemistry. US is not mandatory for a paper explicitly restricted to
potential-energy mechanisms; quantitative kinetic/free-energy claims need
appropriate free-energy evidence. No smearing, gates or force tolerances are
changed to meet the calendar. Existing jobs continue. The new conditional
continuation is `qmmm/pa6_flipped_model30/12_nterminal_relay/`, tmux
`nyl12_w033_relay_after_minimum_20261123`: after qualified full-shell minimum
acceptance it runs two one-CPU branches from the same window 033 exterior,
back toward E–S and forward toward acylenzyme. Accepted unbiased endpoints and
interior physical-PES maxima trigger native dimer TS searches. Movies, energies
and proton/bond-distance plots refresh after each window. Failure or unresolved
minimum curvature prevents dispatch. First-order TS verification remains a
subsequent scientific gate. See the campaign README and live watcher receipt.

### Proposed bounded deacylation pilot — not launched

The PI asked whether difficulty retaining TI1 reflects poor adaptation to PA6
and whether to include acylenzyme hydrolysis. Poor differential stabilization
of a reactive state is plausible, but loss of a TI1 well does not establish it:
substrate pose, protonation, hydration, Hamiltonian and environment also matter,
and hydrolysis need not have a discrete TI1. The force-converged window 033
candidate is currently awaiting curvature classification. Low polymer conversion
does not by itself identify an elementary rate-limiting step.

Include deacylation in the complete-cycle hypothesis and eventual design
assessment. Andjelkovic et al. (2023; manuscript R41) found deacylation rate
limiting for human asparaginase III, with the Thr N-terminal amine activating
water and subsequently returning a proton to the catalytic oxygen. This is
a testable Ntn-hydrolase analogy, not an established Nyl12 mechanism. Lauko
et al. (2025; DOI 10.1126/science.adu2454) also show why an acylation-competent
design can require further improvement of deacylation.

Keep the accepted one-week acylation path/TS work as the priority. The bounded
addition proposed here is hydration/acid-base geometry screening and one fast
QM/MM deacylation feasibility path from an accepted acylenzyme, preferably the
endpoint generated under the current window 033 boundary conditions. Do not
start a second large PMF campaign or imply that a short hydration simulation
has converged. Existing accepted acylenzyme coordinates can support a preliminary
geometry audit, but energies from different frozen exteriors must not be spliced.

1. Retain the full amine-ended PA6 fragment initially; establish the complete
   proton inventory, endpoint chemistry and product contacts. A retained
   product can affect water access, electrostatics and proton transfer.
2. Screen water access and donor/acceptor geometry. A participating water must
   enter the QM region; the current 98-real-atom QM region contains no water.
   Use the unchanged 12 Å movable region. Short QM/MM MD avoids inventing
   classical acylenzyme parameters; it requires explicit preparation and does
   not automatically provide adequate equilibration. The first base hypothesis
   is neutral Thr233 Nalpha; Asp272/Asp274 or additional-water alternatives
   require supporting geometry/protonation evidence.
3. Test water activation/attack followed by proton return to OG1 and ester
   cleavage, allowing asynchronous changes. Retain any supported TI2 minimum
   without requiring one. Qualify endpoints and saddle connectivity using the
   same scientific gates as acylation before a mechanistic claim.
4. A product-released acylenzyme is a separate model. Deleting the amine changes
   composition and possibly charge, QM boundaries and counterions. Its raw PES
   cannot be joined to the retained-product profile. If product release becomes
   necessary, define a separate dissociation/free-energy leg with consistent
   standard state, concentration, solvation and restraint corrections (R48).
   Do not add a generic gas-phase translational-entropy bonus to a solvated
   enzyme profile. A retained-product path is conditional on that bound state.

Design should improve activation free energies relative to each corresponding
basin, rather than maximize TI1 or acylenzyme binding. Initially protect the
catalytic/maturation residues and tetramer interfaces and vary the existing
neighboring design positions to tune oxyanion donors, proton-relay geometry
and water access. Extra positive charge could also over-stabilize E–S, alter
base pKas or obstruct substrate/water, so no specific substitution is justified
from electrostatics alone. Use a single symmetry-consistent sequence across
the E–S, acylation TS, acylenzyme/water, deacylation TS and relevant product/open
states. LigandMPNN supports weighted tied decoding; alternatively generate
from individual states and cross-score the same sequences on every state.
Do not use residuewise majority consensus or interpret model likelihood as
an activation barrier. Missing TSs currently limit the strength of this design
assessment. The complete cycle is a stronger eventual narrative; it is not a
guaranteed additional deliverable within this week.

The preferred paper centers on a defensible acylation result and its limitations,
using earlier unsuccessful direct routes as controls. Cite the separate PA6
parameter release within its validated scope; do not imply that every modular
polyamide/state has been validated. Treat sequence-design/MD as exploratory or
a later paper unless a coherent mechanistic connection is demonstrated. If no
supported pathway emerges by the decision deadline, do not claim the mechanism
is solved or presume a collection of negative calculations is publishable.
Assess whether validated parameter/structural results support a distinct,
novel standalone contribution; otherwise finish an honest handoff/data package.
The most useful experimental follow-up includes clean, expression/maturation-
controlled activity measurements and direct testing of the modeled soluble
substrate/cleavage hypothesis. Current crude-lysate MS results cannot alone
establish intrinsic catalytic improvement.

The matched MD follow-up uses mask_07 seed 04, mask_03 seed 03, V111C seed 02
and WT seed 04, two independent 100 ns replicas per sequence. All four initial
sites per structure have negative projected phi under consistent scissile-amide
mapping; WT seed 04 has four Tier A sites but zero native-Thr passes. Starting
poses differ, and no orientation-holding restraint is applied. This tests a
possible Thr N-terminal relay to the leaving amide nitrogen. It does not assume
closer N–N proximity: WT negative-phi endpoints average 4.65 Å versus 3.61 Å
for positive-phi endpoints, so direct transfer from these flipped starts needs
rearrangement or possible water mediation. A shorter distance alone does not
prove a mechanism. Historical poses are retained; no additional positive-pose
MD is launched. Mask lists, canonical mutations and
selection evidence: `docs/ligandmpnn_mask_site_rationale_20261123.md`.
Actual job status is in `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/queues/`;
process completion, geometry retention and sampling convergence remain separate.
Violet WT/broader replicas and Azalea focused/V111C replicas have passed heating,
density and weak-restraint stages and entered unrestrained stages; second replicas
are serially queued with one replica of every sequence on each host. Azalea's
generic CUDA startup failures were preserved and corrected using the previously
validated Blackwell build. Read `AMBER_WORKSTATION_RUNTIME.md` before any future
GPU launch; the campaign host launcher now rejects incorrect executable hashes
and host-specific configs before launch. The original 12 Å Hessian continues.

All three unrestrained 12 Å DFTB3/MM releases independently force-converged.
Window033 retains distorted tetrahedral-like geometry: C–OG1 1.740 Å,
C–N 1.490 Å, C–O 1.309 Å, Nalpha–HG1 1.050 Å; maximum force
7.94e-6 Eh/bohr. Window032 is noncovalent with HG1 on Nalpha
(C–OG1 2.661 Å), while028 returned to E–S-like geometry with HG1 on OG1
(C–OG1 3.286 Å). None has yet been certified as a minimum by curvature.

The original window 033 controller correctly stopped before a full Hessian:
exact repeated gradients/serialization parity passed, but .001/.0005 Å
finite-difference columns failed step stability. The bounded follow-up completed
72 native analytic-gradient evaluations. All nine sampled coordinates pass the
unchanged 1e-4 relative screen across three adjacent finer-step comparisons.
At h=.000125 Å, the problematic Thr234 HB yz reciprocity difference falls to
5.08e-9 Eh/(Å² amu), versus5.89e-3 at h=.001 Å. These are precision diagnostics,
not minimum evidence by themselves.

**Full 12495-DOF Hessian is actually running** in tmux
`nyl12_relay_w033_full_hessian_20261123`, using h=.000125 Å (largest tested
step on the stable plateau), eight nice single-CPU workers, zero electronic
temperature and the unchanged 4165-atom movable shell/fixed exterior.
Results/progress: window 033 `minimum_verification/06_full_hessian/`;
current interpretation pointer: `minimum_verification/verification_summary.json`.
Require qualified zero-negative-mode inertia before a minimum claim; retain
unresolved signs honestly. No new reaction path has been started.

Final optimizer VMD bundles retain the canonical paths. Endpoint-inclusive
downsampling and an explicit accepted-restart frame correct the previous
omission of the exact endpoint; full raw trajectories/restarts are retained.
See `neb_us_nterminal_relay_release_20261123.md` for measured export checks.

## Earlier checkpoint — calibrated Thr gate and relay releases, 2026-11-23

The PI approved calibration and three unbiased DFTB3/MM releases of the
N-terminal branch sampled in existing NEB-US data. No thresholds, proton
inventory, Hamiltonian, QM partition or movable shell changed.

1. **Gate calibration completed in tmux.** All eight cryo 9Z0B sites, both
   room-temperature 9Z0A sites and four derived n12c sites pass the original
   gate. The original sites span two diffraction datasets, not ten independent
   experiments. Maximum OG1/N shifts are 0.1232/0.0841 Å versus fixed limits
   0.75/1.25 Å. Known-active K59D still has zero joint Tier A/B+native passes
   after relaxation; structural compatibility is not activity validation.
   See `native_thr_gate_calibration_20261123.md`, including reference-fit
   sensitivity, MD scaffold diagnostics and the n12cH A228/P228 edit audit.
2. **Relay releases on Violet.** Sessions are
   `nyl12_relay_ti1_window032_20261123`,
   `nyl12_relay_ti1_window033_20261123` and
   `nyl12_relay_early_window028_20261123`. Two distinct tetrahedral-like
   seeds and one earlier-attack seed have actual Nalpha–HG1 bonds with H1/H2
   also retained on Nalpha. All retain the original 4165 movable atoms in
   the 12 Å region and fixed exterior; no positional or reaction-coordinate
   potentials act in the movable region. QM charge −2, topology, 98 real QM
   atoms and eight link atoms are unchanged. Initial finite energy/gradient
   and SCC checks passed; force convergence and chemistry remain to assess.
   See `neb_us_nterminal_relay_release_20261123.md`. The campaign is
   `umbrella_sampling/nterminal_relay_releases_20261123/` under the existing
   direct-path root. Verified stripped VMD bundles retain explicit atom maps.
   Window032 subsequently converged with maximum force 8.50e-6 Eh/bohr:
   C–OG1 2.661 Å, C–N 1.371 Å, C–O 1.262 Å, HG1–Nalpha 1.079 Å.
   It is a force-accepted noncovalent NH3+/alkoxide configuration, not the
   original neutral-N/OG1–H E–S state. Its Hessian has not been evaluated;
   no extra window 032 curvature campaign was launched. Other seeds continue.
3. **Do not duplicate the 35.81-path campaign.** Its 35.80895 kcal/mol value
   is NEB ΔE, not ΔG‡. That path family already supplied the refined
   46.41-kcal/mol candidate and completed US campaign with apparent PMF
   height 47.0475 kcal/mol. Reference geometries and observables differ.
   Finish the bounded existing-candidate Hessian diagnosis before claiming a
   first-order TS; no duplicate TS/US run or sampling extension is launched.
4. **PI-authorized offline follow-up (23 November):** leave all three original
   optimizations running to their convergence checks. Persistent tmux session
   `nyl12_relay_minimum_w033_20261123` passed all 40 prerequisite hashes and
   is waiting in window 033's `minimum_verification/00_queue/`. It starts
   window 033 verification as soon as that seed converges;
   it does not wait for the other two. If tetrahedral-like geometry survives,
   freeze the accepted restart and verify repeated analytic gradients,
   serialization parity and finite-difference step stability before computing
   the full 12495-DOF active-region Hessian (at most eight low-priority CPU
   workers, no GPU). Require zero resolved negative eigenvalues and no unresolved
   signs under the recorded numerical screen before accepting a local minimum
   in this fixed-exterior model. No QM-only Hessian or smaller movable shell
   substitutes for that check. Numerical ambiguity is reported as unverified;
   neither the skew diagnostic nor force convergence alone proves a minimum.
   No new reaction-path campaign starts automatically. The original 10000-cycle
   per-seed budgets remain in force; an exhausted budget is a continuation
   checkpoint, not convergence. See `relay_minimum_queue_20261123.md`.

   **Decision after release:** classify all proton occupancies, C–OG1, C–N,
   C–O and force convergence. Retention warrants further stationary-point
   classification; reversion or onward acylation instead constrains the path.
   A discrete TI1 minimum is optional. HG1→Nalpha followed by another N–H
   donating to the leaving N can be legitimate relay chemistry.
5. **NylC MCCE reruns completed, still unqualified.** 5XYG/7YU0/7YU1 native
   completion records retain `assigned_states:false` and
   `sampling_convergence:false`. Bounds, poor target fits and the known
   proline-template inconsistency require matched-input curve/parameter QA.
   This remains mechanism-agent's bounded task. Muse means **Glimmer 30B**;
   no model installation or benchmark was launched this turn.
6. **Progress briefing:** editable PowerPoint, PDF, slide JSON, source hashes
   and key figures are in `docs/presentations/nyl12_progress_20261123/`.
   The deck includes the current L4-only modular FF release, Eyring table
   and explicitly marked VMD snapshot placeholders. The laptop pull includes
   the briefing, calibration and new relay inspection bundles.

## Earlier checkpoint — executable priorities, 2026-11-23

**Scientific conditions remain unchanged:** no new protonation assignment,
smearing, QM partition, or reduction of the 12 Angstrom movable shell.
The structural-design pilot remains exploratory; no new design MD is queued.

1. **NylC pKa controls are running.** Repaired 5XYG/7YU0/7YU1 inputs passed
   native MCCE preparation and conformer generation; three independent Violet
   tmux jobs are in DelPhi step3. Their accepted coordinates were reused byte
   for byte. The REMARK failure was our input-parser issue; allowed metadata
   are now preserved and audited, while unknown/invalid inputs still fail.
   Native completion must be followed by target-curve/chain/parameter QA;
   no protonation states are assigned automatically.
2. **TS electronic gate corrected; endpoints independently accepted.** Amber
   recovered 5/9 inner SCC trials in the saved minus/plus descents. The old
   regex misclassified them as terminal failures. Fresh native gradients
   accept E-S and acylenzyme endpoints. Three fresh TS gradients are identical;
   maximum force9.24e-5 Eh/bohr meets the original1e-4 criterion. The candidate
   lies46.41 kcal/mol above the relaxed E-S endpoint; the route remains too
   unfavorable to establish the desired viable mechanism. No first-order TS
   or formal IRC is yet certified.
3. **Next TS calculation: a bounded12-evaluation Hessian diagnostic, not yet
   launched.** At the unchanged candidate recompute Thr930 HB x/y/z columns
   at h=6.25e-5 and3.125e-5 Angstrom. Require step stability and improved
   mixed-derivative reciprocity before forming any separately retained
   derived Hessian and reassessing the unresolved soft-mode signs. Retain
   the original matrix and all4165movable atoms. Do not reduce the global
   uncertainty threshold simply to obtain acceptance. A new full matrix
   needs a demonstrated precision plateau first.
4. **Mechanism-agent next:** a human-reviewed, tool-backed NylC pKa QA case.
   Audit these MCCE curves and run/compare PROPKA on the same prepared inputs;
   distinguish mature Thr267 NTR (5XYG C/D) from continuous precursor chains
   (7YU0/7YU1 A/B). Return a source-linked PI decision packet. Keep GPU model
   servers off between requests; no new model benchmark or autonomous
   scientific dispatch is queued.

Detailed evidence, tolerances, source diagnosis and caveats are in
`nyl12_stationary_validation_loop_20261121.md`. Current NylC job pointers are
`analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/{5XYG,7YU0,7YU1}/current.json`.
The native five-point check is under
`.../validation_loop_20261121/diagnostic_repeat_20261123/`;
`interpretation.json` explicitly distinguishes original from tighter diagnostic
force criteria. All earlier queues below are dated history when inconsistent
with this current section.

## Structural-design checkpoint — 2026-11-22

The bounded 19-versus-29-position mask comparison is complete: 55/55 model
seeds force-converged (52 new, three reused WT). Joint Tier A/B + native-Thr
sites are broader 12/80, focused 8/80, WT 1/20. For the negative-phi pose,
counts are 11/80, 3/80 and 0/20. Four selected sequences per arm and five
seeds per sequence provide descriptive evidence only. Broader mask_07 has
6/20 joint sites in 3/5 seeds; focused mask_03 has 3/20 in 2/5, all negative
phi. K59D has zero joint passes despite known activity in noisy lysate/MS
assays; filters remain unvalidated as activity classifiers. Next bounded
structural decision: inspect these candidates against WT and V111C before
selecting any matched MD follow-up. Do not widen the campaign or infer TS
stabilization automatically. Exact sources and reporting-repair receipts are
in `design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/analysis/`.

## CURRENT — completed structural pilot and matched mask comparison, 2026-11-22

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

# Nyl12–PA6 mechanism: decision roadmap

Status: working mechanism subproject plan, updated 2026-11-22. The overarching
`MASTER_PLAN_NYL12_COVALENT_REDESIGN.md` remains controlling for redesign;
this document is its **single forward-looking mechanism execution plan**.
The append-only action journal, run manifests and manuscript source ledger hold
execution history. Earlier dated mechanism plans are provenance, not queues.

## Closeout — historical neutral Thr/Asp272/Asp274 direct-route NEB-seeded US (2026-11-22)

All 63 windows completed six accepted 5 ps blocks (30 ps/window): eight Violet
workers exited 0 on November 22 at 10:43–10:47 UTC and 24 Azalea workers exited
0 on November 21 at 11:55–20:49 UTC. The latest receipts are in
`umbrella_sampling/dispatch/coverage_first_20261120/{violet,azalea}/` within
the campaign root stated below.
Only local analysis PID 4100242 and sync PID 4098079 were stopped with SIGTERM;
both were absent afterward, as recorded in
`umbrella_sampling/dispatch/stop_20261122/violet_monitor_stop_receipt.json`.
No campaign Sander remained active on either host. Completion is not sampling
convergence.

The refreshed apparent PMF maximum is 47.0475 kcal/mol versus 46.9414 in the
same E-S density bin. Overlap is connected, target ESS is 7.29, and windows
028/030/032 have drift flags. This PMF is not converged and is not a validated
barrier. The scientific decision is to stop further investment in this forced
direct route; preserve raw trajectories, restarts, source snapshots and the
negative finding, with no extension queued. This applies only to the historical
neutral-Thr, ASP272/ASP274 direct-NEB-seeded US under
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/`.
Other work, redesign and the separate ASH274 NEBs are unaffected.

The verified movie is at
`umbrella_sampling/vmd_dynamic/view_umbrella.vmd`: 30 frames/window,
1,890 frames and 19,967 atoms. Every selected frame reproduces the four
recorded reaction distances within 0.00005 Å; bundle checksums pass.
Window concatenation is not a continuous dynamical reaction trajectory.

The PI's model-capability and research-efficiency questions are addressed in
[`mechanism_agent_model_and_strategy_20261120.md`](mechanism_agent_model_and_strategy_20261120.md).
It records replacement candidates, hardware/license constraints and lessons
from published scientific agents. These are advisory recommendations; no model
deployment, new mechanism branch or change to the active queue is implied.

## Current structural validation alongside mechanism work

PI-approved on 2026-11-22: 45 Boltz-2 models (WT + eight unique LigandMPNN
variants, five seeds each), four zwitterionic PA6 L4 copies per tetramer, then
fully unrestrained solvated RESP/GAFF2/ff19SB relaxation and all-four-site
geometry comparison. Violet/tulip predict; geranium relaxes under renewed PI
host authorization. Existing QM/MM/US jobs are untouched. Details/live counts:
`design/nyl12_ligandmpnn_pilot/boltz_validation/README.md` and `status.json`.
This advances the exploratory E-S redesign lane; it does not resolve or replace
the reaction-mechanism work, certify a TS, or authorize new candidate dynamics.

## Active queue and next gates

### Authorized — WT/known-active filter calibration and design planning (2026-11-21)

**Execution update:** PI approved a single-parent, 24-sample LigandMPNN pilot
using WT Model23, restricted to violet/azalea while current QM/MM/US/MD
continue. WT40/41 failed occupied-site native-Thr checks; only WT23 passed the
bounded screen of seven existing Tier A/B WT models. The accepted preparation
manifest is `design/nyl12_ligandmpnn_pilot/accepted_parent/run/preparation_manifest.json`.
Completed: 24 samples, eight unique sequences, 16 duplicates, 4-6 substitutions each.
All sequence-integrity and independent tensor/FASTA checks passed; results are
in `accepted_parent/run/results/` and the campaign README.
Retain Model23’s earlier QM/MM deterioration as a negative result. This is an
exploratory single-parent E-S sequence pilot, not validated-TS design; candidate
prediction/relaxation, matched MD and activity remain unvalidated. See
`docs/nyl12_ligandmpnn_pilot.md` for the fixed mask and sampling budget.

The PI authorized a parallel exploratory engineering lane. Start by identifying
the exact measured WT/mutant labels and reanalyzing existing static/MD outputs
under their frozen filters. Keep every current QM/MM and US calculation running.
The initial authorization was retrospective analysis and planning. The later
approved LigandMPNN pilot above is now complete; new predictions, MD and local-
backbone generation remain outside that completed pilot. See the updated master plan
and `plm_and_filter_calibration_20261121.md`.

Keep experimental PA6/polymer labels separate from PA66-only, unmeasured and
uncertain variants. Absence of activity data is not inactivity. Match comparison
strata by ligand/register, predictor, relaxation, FF family and protonation;
do not pool mismatched protocols into a mutation effect. Retain missing fields
and exclusions explicitly. Report current positive-control retention and
individual gate failures before considering revised rules. No cutoff tuning
or claims about specificity/AUC without appropriate independent labels.

The analysis must retain immutable input/source snapshots, experimental-label
evidence, environments, hashes and process/scientific acceptance. Run persistent
analysis in a low-priority tmux job using existing dependencies. Use the results
to plan 24–48 fixed-length, symmetry-tied LigandMPNN variants plus controls,
protecting maturation, proposed catalytic machinery and essential interfaces.
The final editable positions and ranking policy depend on calibration; do not
invent a successful outcome or treat Tier occupancy as measured activity.

PI assay-quality clarification: these are crude-lysate, mass-spectrometry-based
product measurements, with considerably discordant duplicate experiments for
**K59E**. Primary workbook:
`/data/jpc/enzyme-substrate-modeling/config/cb150_nyl12_ssl_sites_1_4_pa6_pa66_pmaxa.xlsx`.
Retain both experiments and establish replicate type, backgrounds, WT controls
and normalization before interpreting a fold difference. Do not select one
run, silently average discordance, or treat a derived positive signal as clean
ground truth. Existing single-value activity comparisons are provisional;
the completed structural count/gate audit does not validate an activity ranker.
Workbook K59E PA6 L2 values are 0 and 1.1465832774; their arithmetic mean is
the legacy 0.5732916387. The two entries remain separate in the source audit.
The manuscript reports independently re-screened hits and triplicate MS
measurements with internal-standard/product calibration, but source rows cannot
yet be assigned to those replication levels. The PI believes the updated Figure
7 uses the same data; experimentalist confirmation is deferred a few days.
Keep that linkage provisional without blocking structural audit or planning.

The first static audit has completed in tmux with 900/900 prespecified models,
zero Tier A/B parity mismatches and exact historical count reproduction.
Outputs: `analysis/nyl12_filter_calibration/{report.json,tables/,plots/}`.
WT PA6 A/B retention is 7/50; V111C 1/50; K59D 18/50. Every variant retains
at least one A/B model. The Burgi-Dunitz angle is the dominant rejection gate,
not a validated activity predictor. All thresholds remain unchanged.
The approved specification is `docs/nyl12_ligandmpnn_pilot.md` and the
accepted-parent manifest above: 24 raw samples from Model23 at temperatures
0.1/0.2, with two independent RNG seeds and six samples per job. Sequence ties
are explicit across A-H with 0.25 weights; fixed-length design uses the declared
19-position mask and protects maturation/catalytic context. No inference report
is available yet. These static results do not certify the newer RESP/GAFF2
mutant-MD filters or an experimental activity ordering.

Exploratory E-S-based sequence engineering no longer waits for a stable TI1 or
verified TS. TS-directed claims still require verified stationary points; TI1
is optional if the actual pathway has no corresponding minimum.

The PI's NylB/poly(MAXA)-to-MAXA-selective concept is recorded as a separate
future hypothesis, with product/cleavage mapping and patent-scope review needed
before execution. Scaffold-constrained Nyl12 active-site rebuilding is a later
comparison, not a switch away from the present calibration-first sequence lane.
Preserve the mature scaffold and assembly rather than physically deleting
uncharacterized catalytic contributions. Mechanism-agent remains a bounded
evidence/automation aid: autonomous scientific planning is unqualified by the
current benchmark; no new model deployment is implied.

### Queued — final matched WT occupancy and relay analysis (2026-11-21)

The PI requested that the remaining WT MD finish and analysis run automatically.
Tmux `nyl12_final_matched_wt_analysis_20261121` now waits on all twelve accepted
100 ns replicas under `10_wt_protonation_md/analysis_after_md/`. Initial live
preflight accepted eight and correctly left four pending. It checks actual stage,
terminal-job and 100 ns receipts plus production-topology/source hashes; a 10 ns
completion marker alone is insufficient. Upstream failures stop analysis, never
resubmit or alter MD. Polling is once per minute with a 48-hour waiting limit.

After completion, two low-priority single-thread CPU stages run in parallel:
all-frame occupancy/native analysis for the matched twelve replicas, and the
four newly completed relay screens. The eight previously screened matched
replicas are reused with verified hashes. Historical unmatched ensembles are
excluded from the three-state comparison. Grouping explicitly distinguishes
neutral and cationic Thr N, not just Asp residue names.

Report Tier A, exclusive Tier B, inclusive A-or-B, signed phi and native-Thr
intersections. Summarize four independent replicas per state with sample SD/range
and 10 ns block stability; do not treat correlated frames or four enzyme sites
as independent replicas. Join unique-frame relay counts to the **all-frame**
native measurements rather than a twenty-frame native shortlist. Fresh trajectory
hashes must match the reused relay hashes before comparison.

Outputs remain inside the existing MD campaign: `plots/tier_occupancy.png`,
`tables/occupancy_by_{replica,state}.csv`, `occupancy/summary.json` and
`relay/comparison/`. Read `waiting/status.json` for progress/errors. These are
fixed-state structural occupancies, not measured pKas, catalytic rates or a
validated mechanism. Current MD and NEB-seeded US remain unchanged.

The PI also authorized geranium and tulip on November 21. Initial default-key
authentication failed; the explicitly authorized aster key subsequently worked
for both hosts. Use `~/.ssh/id_ed25519_aster` with explicit identity selection.
Geranium has 48 logical CPUs and tulip 24; both report about 503 GiB RAM and
two idle RTX 6000 Ada GPUs with 49,140 MiB each at inspection. Amber Sander and
pmemd executables were located on geranium; tulip's Amber location remains
unresolved. See `provenance/workstation_resources_20261121.json`. No transfer,
installation, native compatibility run or job migration was performed. Recheck
contention and the exact environment before using extra hosts for a specific
calculation; existing MD/US and post-MD analysis remain on their current hosts.
PI resource restriction: use geranium and tulip simultaneously only for work
expected to finish within a few hours. Before every launch, check other users'
sessions, CPU/GPU workloads and queued work on each proposed host. An earlier
idle check is not permission to launch later. Ask the PI if availability or
expected duration is unclear; preserve these shared resources for other users.

### Authorized next step — leaving-group proton-relay screen (2026-11-21)

The PI approved screening existing MD and selecting the best proton-delivery
configurations before launching another reaction-path campaign. Canonical new
analysis: `qmmm/pa6_flipped_model30/11_leaving_group_relay_screen/`.
Keep the current NEB-seeded US and all other scientific jobs running unchanged.
Nyl12 remains primary; no new NylC or polymer-parameter campaign is implied.

The two initial workers completed fourteen accepted 100 ns replicas. A second
queue completed two original-pose WT 100 ns replicas, whose cached signed phi
is negative throughout. Combined coverage is 1.6 μs and 32,000 coarse frames
at 50 ps spacing. The original-pose files save every 50 ps, so their explicit
execution stride is one; the fourteen flipped-cohort files save every 10 ps
and use stride five. Incomplete or unavailable replicas remain excluded.

The first selection/refinement completed successfully but mostly identifies
zero-intermediate OG1→leaving-N direct-route **controls**, not newly discovered
relays. Do not treat `selected/selection.json` as evidence for a new mechanism.
The assisted-only review completed exit zero in tmux
`nyl12_assisted_relay_candidates_20261121`, with local source-frame refinement.
Final outputs: `assisted_review_refined/{selection.json,summary.md}` and
`candidate_readiness.json`. There is **one unique strict assisted Tier B frame**:
matched NH2/ASP272/ASP274 replica 03, frame 9386, with a
Lys156→water29282→water4517→leaving-N network. It fails current native-Thr
displacement gates (OG1 2.396 Å, Nalpha 2.958 Å). A Tier A single-water
Thr-OG1→water9397→leaving-N nearpass in ASH272 replica 1 frame 7316 misses
the H–N cutoff by only 0.00964 Å, but also fails native displacement/psi gates.

Twenty diverse shortlisted frames were audited against native Thr; none passed
all those criteria. Three unchanged restarts were exported for diagnostic
inspection, **not accepted as new QM/MM seeds**. Automated "joint qualified"
means hydrogen-bond legs plus Tier A/B only; the explicit readiness record
adds native/chemical acceptance. Do not count the repeated coarse/refined
observation of frame 9386 as two independent candidates. The native audit is
bounded, so this does not prove that every original frame lacks a viable relay.
Review these geometries and native criteria before another path or MD campaign;
no automatic relaxation or altered-protonation calculation was launched.

Before any selection ran, review found that per-frame route deduplication
preceded full-network priority annotation. The corrected selector is retained
separately in `selection_revision_02/source_snapshot/` with a regression test.
Only the original orchestration parent was paused; neither data-reading worker
was interrupted. The replacement retired that parent after both workers finished
exit zero. Its expected retirement status 143 is execution provenance, not an
MD/US or coarse-screen crash. Corrected baseline selection/refinement exited zero.

The chemical question is whether an observed configuration supplies a plausible
complete proton network to the **leaving amide nitrogen**, together with useful
attack/oxyanion geometry. A persistent TI1 is optional. A distance-difference
target, favorable Tier score, or hydrogen-bond path alone does not demonstrate
proton transfer, a stable intermediate, or an acceptable activation barrier.

Execution is deliberately bounded:

1. Inventory completed WT production trajectories and per-topology proton
   inventories. Reuse original atom maps and geometry implementations; never
   transfer atom indices across protonation-dependent topologies unchecked.
   Include the historical neutral/charged-N and ASH272 ensembles and completed
   matched ASP/ASP, ASH274 and cationic-N replicas where available. Exclude
   actively written segments and clearly label any partial completed coverage.
2. Screen the bound active site for zero-, one- and two-water directed paths
   from Thr OG1, catalytic Nalpha, protonated Asp272/Asp274 and Lys156 to
   leaving N, using actual topology-owned hydrogens and periodic images.
   Standard structural criteria are D–A ≤3.5 Å, H–A ≤2.5 Å and D–H–A ≥120°.
   Keep complete directed networks separate from heavy-atom proximity and
   orientation nearpasses. Neutral Nalpha donation is conditional on prior
   activation; amide-N acceptance is a reaction-preparation hypothesis.
   Never invent a proton, rotate H atoms, or infer a neutral Lys microstate.
3. Use a 50 ps coarse screen when saved frames are 10 ps apart, then inspect
   original saved frames around selected events. Report d1/d2/d3, Bürgi–Dunitz
   angle, **signed** Herschlag phi, substrate N geometry and native Thr
   displacement/torsion checks. Screen both pose signs; rank chemical pathways
   within each orientation instead of requiring the historical direct-route
   sign universally. Rank strict directed paths with Tier A/B geometry first;
   retain clearly labeled nearpasses and native-gate failures for diagnosis.
4. Select at most two primary structural candidates plus a small, diverse
   reserve set, avoiding adjacent correlated frames. Export unchanged source
   snapshots with topology/frame/time/water identities and full provenance.
   No new MM, QM/MM or altered-protonation launch is part of this screen.
5. Review the selected complete acid–base networks. The next separately
   specified path must allow attack, C–N cleavage and individual proton-site
   occupancies to vary independently. Do not lock Nalpha in NH3+ while asking
   it to donate onward, or confuse proton labels with chemical occupancy.
   Use sparse matched-embedding DFT and one targeted QM-region sensitivity
   test only on genuinely discriminating new configurations.

Queue CPU-only analysis in tmux with immutable source/config snapshots and
input hashes. Use a small number of low-priority single-thread workers; do not
consume GPUs or stop US to make room. Completion, structural acceptance and
mechanistic validation remain separate records. The completed screen, including
negative results, determines the next chemistry decision rather than another
unbounded collection of reaction-path attempts.

### Latest NEB-US MBAR PMF refresh — completed 2026-11-22

The fresh MBAR calculation completed with exit 0 using all 63 windows and six
accepted 5 ps blocks per window (30 ps/window, 1,890 ps total; 37,800 raw and
7,871 thinned frames). Sampled overlap is connected and target ESS is 7.2928;
all 50 connected 5 ps block bootstraps completed. The apparent q_acyl maximum
is 47.0475 kcal/mol relative to the same fixed E-S density bin as the previous
46.9414 kcal/mol result (+0.1062 kcal/mol). The fixed-peak-bin bootstrap 95%
interval is 45.4070–47.8696 kcal/mol; this is not comprehensive PMF uncertainty.
Windows 028, 030 and 032 are flagged for drift.

Frozen source/input hashes, same-bin numerical-range regression, independent
2D marginal and publication checks passed. Detailed outputs are in
`qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/umbrella_sampling/analysis/pmf_refresh_20261122/`.
The stable sibling `analysis/pmf_latest/` now points to its `compact_export/`;
the prior export is preserved.

Sampling convergence remains false, and this is not a validated activation
barrier. No sampling jobs were stopped or changed, and no additional sampling
extension was inferred from this analysis.

### Prior checkpoint — preliminary NEB-seeded PMF and continued sampling (2026-11-21)

Read [the frozen preliminary findings](neb_us_preliminary_findings_20261121.md).
All 63 windows have production data (1,190 ps at the frozen checkpoint);
connected overlap permits a preliminary q_acyl marginal with a 46.50 kcal/mol
maximum relative to the specified E–S density bin. Target ESS is only 28.28;
sampling convergence and a mechanism-specific activation free energy remain
unestablished. The source 0 K regular NEB maximum is 35.80895 kcal/mol at
image 8, with the acylenzyme endpoint at +4.92846 kcal/mol.

HG1 can bind Thr N while tracking both native distance-difference umbrellas;
window 031 has mean q_H 0.4417 Å versus target 0.4568 Å despite HG1 binding
Thr N. This is observed CV degeneracy, not evidence that the same restraints
need strengthening. Keep all branches and the unchanged 12 Å movable region.
The PI requested continuation: 31 workers were active across violet/azalea at
the live audit. Finish the existing minimum 30 ps/window, then assess block
stability, ESS, overlap and branch exchange before any new campaign or claim.


### Current — completed Hessian and targeted validation, November 21, 2026

The full matrix completed with one chemically relevant negative mode in the
stored symmetric Hessian, but numerical consistency and downhill connectivity
remain unverified. Read the completed-outcome section in
`nyl12_stationary_validation_loop_20261121.md`. A bounded two-worker directional/
localized-asymmetry diagnostic is now running with conditional candidate-mode
descents; the original failed gate is retained. Reuse this matrix rather than
repeat the full 12,495-column calculation. The diagnostic pointer is alongside
the original validation current.json, named diagnostic_current.json.

Both local model runs completed. Paired semantic review shows a narrow Mistral
critical-case improvement but more major scientific errors; neither model is
qualified for scientific planning or dispatch. Original response text, format
compliance and scientific grading remain separate in the private paired review.


### Authorized — qualified model, stationary-point loop and filter calibration

The PI approved official Apache-2.0 Mistral Small 4 setup on azalea and a
matched GPT-OSS comparison. The private agent repository holds the frozen
14-case benchmark: six historical development cases and eight newly authored
held-out-in-workflow cases, including positive TS/PMF controls. Auditor rubrics
are excluded from inference. Source, precision and tokenizer/runtime identities
are retained; downloading/serving is not scientific qualification. Scientific
jobs are preserved. See the private `benchmarks/mistral_small4_20261120/` and
`docs/QUALIFIED_SCIENTIFIC_CONTROL.md` for actual status and promotion gates.

The numerical validation loop reuses the existing force-converged, zero-
smearing direct-route Nyl12 candidate with all 4,165 movable atoms (12 Å).
The native DL-FIND gradient writer rounded to seven decimal places. An isolated
precision-only library patch is being verified before a complete movable-region
Hessian, explicit index/sign uncertainty and conditional both-sign downhill
minimizations/endpoint checks. Existing installed engines and Hamiltonian are
unchanged. Read its source-frozen receipts before calling the TS verified;
this route's unfavorable energy is not being treated as a viable mechanism.

The next design deliverable is retrospective recovery of known-active controls,
not a larger candidate campaign. Read
[`plm_and_filter_calibration_20261121.md`](plm_and_filter_calibration_20261121.md).
Keep L4-to-two-L2 as a model reaction and polymer-powder assay labels distinct.
Reuse PA6 parameter work and matched protonation ensembles. Calibrate geometry
and confidence as surrogates before rejecting variants; preserve hard chemical
validity checks. Optional ESM models do not yet have installation or inference
approval. No extra mechanism branch or narrowed movable region was introduced.

### Active — complementary pKa qualification and NylC controls (2026-11-20)

The Nyl12 MCCE pilot completed technically, but its quantitative pKas remain
unqualified. Read `pka_trust_and_active_site_audit_20261120.md`: matched PROPKA
inputs, saturated/poor titration curves, a contradictory His257 fit label and
a proline-template inconsistency are recorded. No MD/QM microstate was changed.

The PI now explicitly requested repaired mature5XYG and uncleaved7YU0/7YU1
MCCE runs. Initial gap reconstruction's soft geometry was inadequate; full
Amber ff19SB/GBn2 refinement of modeled atoms/H preserves original X-ray heavy
coordinates and applies strict bond/clash/sequence/true-terminal gates. Mature
5XYG has separate alpha/beta chains; both mutant precursors retain A/B1–355
with internal267 and their deposited substitutions. No backmutation is implied.

Three independent persistent controllers are launched under
`analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/repaired_comparison_20261120/`.
They wait for exact preparation checkpoints, then refine/audit and run native
MCCE only on passing models. Use each structure's `current.json`/README/status,
not a historical run name. Earlier failed attempts and one builder source-history
gap are retained; authoritative precursor jobs replay from immutable sources.
The pinned native parameter set is a diagnostic comparison, not yet validated
quantitative catalytic pKas. No US/GPU MD process was changed.

After completion, review full curves and matched preparation/electrostatic
assumptions. Do not promote MCCE output merely because the process exits zero.
The QM/MM and US mechanism objectives remain the project's primary goals.

### Completed — first read-only mechanism-agent pilot (2026-11-20)

The separate private `/data/jpc/mechanism-agent` prototype now uses official
local gpt-oss:120b on one azalea GPU with two CPU inference threads. Five
source-hashed evidence cases completed in tmux and produced a PI review packet.
Runtime residency and structural checks passed. Unsupported model claims
were identified; scientific/model capability remains ungraded.
No scientific submission or shell tools are exposed. The six-role ensemble
is still a design; this is its single-model baseline. PI review, then a broader
held-out evaluation, precede any consideration of execution adapters. Existing
US/MD jobs and the mechanism priorities below remain unchanged.

### Launched — matched WT protonation ensembles (2026-11-19)

The user authorized all three WT conditions. The common-box campaign in
`qmmm/pa6_flipped_model30/10_wt_protonation_md/` now runs unrestrained GPU
MM before four independent 100 ns replicas per state, sequentially allocated
to two free Violet GPUs. All three actual preparation audits passed; full
MD acceptance and occupancy results remain pending. Read its README and
frozen orchestration status. These simulations test fixed-microstate
structural ensembles; they do not allow proton exchange or establish pKas.
Proton movement from cationic Thr N to either Asp would reach a proton
assignment already considered, although a different protein/water ensemble
can affect access to the reactive geometry. Do not equate this with finding
a new stoichiometric intermediate.

### Active — distributed NEB-seeded US, approved and launched 2026-11-20

The PI explicitly approved transfer, a native compatibility test and up to
24 CPU workers on azalea. **32 native Sander calculations are verified
running: eight on violet and 24 on azalea**, with positive integration steps.
Static ownership covers all 63 windows once: 16 locally (eight adopted plus
eight subsequent windows), 47 remotely. Native violet jobs were not restarted.
The canonical `umbrella_sampling/` directory is identical on both hosts.
Read its `status.json`, `distributed_launch_acceptance_20261120.json`,
`dispatch/plan.json` and `dispatch/remote_status/azalea/` for durable state.

Azalea uses its installed architecture-native CPU Sander and analysis Python;
no executable binary was copied. The 100-step native test completed normally,
with initial EPtot matching violet at printed precision (difference 0.0000
kcal/mol; allowed tolerance 0.001). All 1,046 transferred input/source hashes
passed. Native method, topology, charges, seeds, restraints, 12 Å/4,165-atom
movable selection and at least 30 ps production/window are unchanged. The
initial24-controller attempt failed before MD due to an Amber libtinfo/tmux
ABI conflict; failed artifacts were preserved. The tested launcher isolates
tmux's library environment while retaining scientific job library paths.
The remote retry generation is execution provenance, not a new reaction path.

At 17:18 EST on November 20, 32/63 windows had three accepted 5 ps production
blocks each and were running block four. The remaining 31 windows were queued;
the product-side gap includes windows 045–052, which had not started. Sampled
overlap remains disconnected. **A full-path PMF/barrier is not yet estimable.**
Current plots are `analysis/latest/{sample_coverage,histogram_overlap,mbar_overlap}.png`.
Middle-window HG1 occupation of Thr Nalpha is retained as competing chemistry;
do not label a resulting projected surface a verified direct-route mechanism.

The PI authorized **coverage-first scheduling**. The tested, separately frozen
controller in `dispatch/coverage_first_20261120/` now waits for each current
native block to finish normally, then starts the missing windows. No native
Sander process was signaled or restarted; all 32 original native PID/start
identities and positive integration steps were verified after handoff.
Only the old Python queue/window controllers are replaced. Their original
receipts and all native data remain; intentional controller termination does
not indicate native-MD failure. Use the new controller/worker status files
and top-level `status.json` pointer rather than the retired controller jobs.
An independent tmux gate holds extensions until all 63 windows have accepted,
diagnostic-ready first 5 ps blocks, then the workers resume balanced extensions
to the unchanged minimum 30 ps/window. No scientific inputs or seeds changed.

A seven-day CPU-only MBAR/plot watcher and seven-day compact results mirror
continue in tmux. The mirror publishes validated production observables/bias
matrices and endpoint views on violet, with acceptance metadata last; raw
trajectories stay on their execution host. The laptop pull follows this same
layout and excludes raw trajectories, transfer staging and archived attempts.
Measured production throughput is about 1.34 ps/hour locally and 2.80 ps/hour on
azalea. The frozen `timing_update_20261120.json` describes the superseded
whole-window queue. Revise its estimates after the newly prioritized windows
complete heating/equilibration. Full first-block coverage still does not
guarantee connected overlap or sampling convergence.

Next gates: verify finite/SCC-complete blocks, proton inventories and competing
branches; obtain coverage across all 63 windows; review connected overlap,
autocorrelation/ESS and block stability before a PMF interpretation. Extend
sampling selectively under the approved protocol. Any new mechanism or
Hamiltonian decision remains a PI discussion. Keep existing NEB/MD/parameter
work intact and use completed results as evidence rather than another pivot.

### Launched — umbrella sampling from the converged direct NEB

The 35.80895 kcal/mol historical band is force-converged at cycle 477,
with ASP272/ASP274, neutral Thr N, 98 real QM atoms and charge -2. It
connects E–S to acylenzyme by direct OG1-proton delivery to leaving N.
A regular-band maximum is not a verified TS, but the images can seed US.
The current ASH274+water band is a different model/seed, not an acid-only
comparison; its original HG1 has already exchanged onto water.

NEB's Cartesian path spans the declared movable QM/MM atoms; it does not
force an umbrella coordinate choice. The authorized native-Sander campaign uses independent
chemical CVs: q_acyl=d(C,N_leaving)-d(C,OG1), and
q_H=d(OG1,Hgamma)-d(N_leaving,Hgamma). Each window receives separate
harmonic terms about the image-derived centers. Monitor all four component
distances and exchangeable proton identities so degeneracy or competing
proton transfer is not hidden by a summed coordinate. Seed a local pilot
from the existing band, then add windows/orthogonal coverage based on overlap
and hysteresis; do not interpret an unsampled narrow ribbon as a complete
two-dimensional free-energy surface.

Alternatively, a path CV can express progress s and distance z from the
reference path using these reactive distances. This needs a reviewed native
implementation or explicit dependency approval; no PLUMED installation is
implied. Avoid restraining every QM Cartesian coordinate to each image.
NEB springs are removed. All applied umbrella/orthogonal biases must be
represented in MBAR reduced potentials; an extra confining wall that cannot
be unweighted defines a conditional PMF and must be reported. Preserve the
reference Hamiltonian and boundary conditions for a controlled PES/PMF
comparison. Production is at least 30 ps/window, with extension based on
overlap, autocorrelation/ESS and block/replica stability, not elapsed time alone.
Eight CPU workers are now launched in the original band's `umbrella_sampling/`
child (see its README, prepared hash and native retry02 receipt). There are
63 target windows with ≤0.10 Å adjacent two-CV spacing, starting from actual
NEB restarts; no Cartesian seed interpolation was used. Each receives
1 ps heating, 5 ps equilibration and six 5 ps production blocks. The exact
4,165-atom movable boundary is preserved, making the PMF conditional on the
frozen exterior. Initial pre-MD filename-buffer failures were retained and
repaired using short local names without changing physics or seeds.
CPU-only MBAR diagnostics update after blocks; disconnected overlap blocks
a global surface. At the launch checkpoint no production was complete.
A lower free-energy barrier is not assumed, and a stable TI1 is not required
for concerted acylation.

Method context: Branduardi, Gervasio and Parrinello, JCP 126, 054103 (2007),
DOI 10.1063/1.2432340; Lonsdale, Harvey and Mulholland, Chem. Soc. Rev. (2012),
DOI 10.1039/C2CS15297E. These support CV/bias methodology, not a Nyl12 mechanism.


### Cationic Thr233 comparison and NEB checkpoint — 2026-11-19

Cationic Thr233 N termini with ASP272/ASP274 should be retained as a competing
WT hypothesis, not assumed superior or substituted across all mutants. Existing
PROPKA site-E estimates (8.93–9.70) motivate this comparison but omit PA6 and
are conformation/coupling sensitive. NH3+ cannot accept HG1 as a neutral
amine base would; its reactive pathway must involve prior deprotonation or
another base/relay. A minority neutral population could still carry catalytic
flux. No new microstate or MD campaign is authorized by documenting this option.

Two charged-N WT Model30 replicas already completed 100 ns with the same
PA6 RESP/GAFF2 release and ff19SB/TIP3P. Their Tier A counts are 487/10,000
and 601/10,000 (4.87%, 6.01%), versus baseline 37 and 238 and ASH274 1,677
and 0. The charged joint Tier/native screen yielded only 0 and 1 frame, with
no temporally supported seed; its one diagnostic frame subsequently failed
post-MM geometry. No charged QM/MM followed. Thus simple Tier A persistence
is more consistent between those two charged replicas, but improved catalytic
competence is unestablished. Reuse the cached diagnostics before new runs.

A bounded WT comparison would include neutral-N/ASP274 baseline, neutral-N/
ASH274 and cationic-N/ASP274, keeping Asp272 deprotonated. The last two have
the same dry proton count and formal charge, making proton location a useful
controlled variable. Existing boxes are not identical (118,829 versus 104,934
atoms; Na/Cl 123/119 versus 110/106), and dry partial charges differ by 0.004 e
from the neutral-Thr recipe's -0.001 e per copy. Preserve these differences;
do not call the historical simulations an identical-solvent proton-relocation
experiment. Validate terminal charge bookkeeping before new matched builds.
Original cached evidence: `docs/md_model30_charged_screen_20261111.md` and
`qmmm/pa6_flipped_model30/04_md_screen_20261111/`.

At 2026-11-20 02:12:36 UTC (Nov 19 local), the extended CG NEB reached cycle
260: physical profile maximum 43.642 kcal/mol, maximum force 8.093 times
target, RMS force 0.182 times target, and zero geometric image resets. It
exactly replayed the 79-cycle pilot before advancing. Uniform24 L-BFGS was
at cycle 426 / 51.405 kcal/mol / 59.842 times max-force target; uniform48 at
219 / 47.187 / 23.572; adaptive24 at 428 / 44.553 / 44.484. Original NEB
was at 734 / 49.799 / 71.076. All five engines were verified running.
None passes the maximum-force criterion; these are provisional PES maxima,
not verified barriers or PMFs. No calculation settings changed for this audit.

Actual native geometric image-reset counts were 50 for original, 5 for
uniform24, 11 for uniform48, 9 for adaptive24 and 0 for either CG. The old
monitor's `reset_count` counts BFGS memory resets and must not be interpreted
as this geometric-reset count. The explicit corrected audit and frozen plots
are under the existing `04_neb/plots/status_refresh_20261120_021148/`.


### Active — long CG identical-start replay (2026-11-19)

The user authorized a much larger CG budget, including starting over if needed.
The finished 79-cycle pilot has no usable optimizer checkpoint. The supported
same-start replay is running in `04_uniform_24_cg/long_run/`, under the same
`04_neb/spacing_tests/` campaign. Session: `nyl12_neb_cg_long_20261119`.
Only iteration ceilings changed: DL-FIND maxcycle 250→32,000 and Amber maxcyc
100,000→300,000. The resulting 224,000-evaluation budget is 128× the pilot.
All physical settings, thirteen control files and tolerances are unchanged;
input equivalence was checked by inverse replacement of those exact tokens.
The pilot and all other NEB jobs remain unchanged.

The first native profile reproduces all 24 pilot energies and Max/RMS gradient
values exactly. An independent tmux checker compares through cycle 79 and
records whether the replay subsequently advances beyond it. Full overlap is
pending; allow about 65 minutes to retrace the pilot, not to converge.
Current plots and movies publish under `long_run/plots/` and
`long_run/progress/latest/vmd/`; recursive laptop pull filters include them.
Read `long_run/README.md` and `replay_verification/last.json` before reporting
this as an exact optimizer restart or a scientifically converged path.


### Running — four parallel NEB numerical-spacing trials (2026-11-19)

**Early progress checkpoint:** uniform24/adaptive24 L-BFGS reached cycle80,
dense48 reached40, and the short CG diagnostic reached78. At approximately
comparable energy/gradient work, all four new bands have zero automatic image
resets; the original band has31. Profile maxima are39.8,55.1,42.0 and46.4
kcal/mol respectively; these are provisional physical-image maxima, not
validated barriers. Maximum-gradient/target ratios are113,81,64 and21; all
four RMS-gradient criteria pass, but no band has converged. CG currently has
the smallest maximum gradient despite a higher profile maximum; this is
promising numerical evidence, not a basis to rank physical mechanisms.
Its short1,750-evaluation budget is nearly exhausted by design. Continue the
planned matched-work checks before selecting a production continuation.
Snapshot: `04_neb/spacing_tests/plots/early_comparison_status_20261119.json`.


The user authorized four violet trials under the canonical
`09_water_relay_acylation/04_neb/spacing_tests/` directory. The retained original
NEB continues as context; it is not a matched restart control.

| Trial | Regular images | 13-control seed | Optimizer and scope |
|---|---:|---|---|
| `01_uniform_24_lbfgs` | 24 | Uniform QM arc length | Primary numerical baseline |
| `02_uniform_48_lbfgs` | 48 | Same uniform controls | Primary image-density comparison |
| `03_adaptive_24_lbfgs` | 24 | Double weight near source images 9–10 | Primary seed-detail comparison |
| `04_uniform_24_cg` | 24 | Same uniform controls | Short curvature-history diagnostic |

All four reuse the immutable cycle-268 band (`04_neb/progress/update_000075`),
periodic DFTB3/3ob at electronic temperature 0 K, QM102/54 real H/eight link
atoms/charge −1, ASP272/ASH274 and neutral Thr Nalpha. Frozen endpoints,
exterior, force tolerances and the complete-residue 12 Å/4,662-atom movable
region are unchanged. No chemical or positional restraint was added. Common
`maxstep=0.02` bohr limits each Cartesian step component. Native initialization
still uses piecewise Cartesian linear interpolation: these trials change its
control distribution or image count, not its interpolation formula.

The installed native code repairs bends **above 90° after the optimizer step
cap**, interpolating neighboring images; reducing `maxstep` does not cap those
repairs. In the original band, reset-associated maxima jumped
35.08→82.74 kcal/mol (cycles 244→245) and 32.83→75.15 (269→270). This is a
concrete numerical diagnosis, not evidence that a lower transient maximum
represents progress. L-BFGS remains primary; CG only tests the role of its
curvature history. CG has a 1,750-evaluation budget (`maxcycle=250` in this
adapter), roughly 79 full 24-image cycles; exhaustion is not convergence.
The three L-BFGS trials retain 112,000 evaluations each.

The 48-image initialization has one explicitly reviewed shared reacting H:
H20318 in image 15 is 1.521 Å from leaving N and 1.590 Å from OG1. The strict
retention gate remains false; `inputs/uniform/shared_proton_review.json` pins
this exact unconverged seed. All other images and spectator checks pass; the
exception neither accepts a stationary state nor adds a chemical restraint.
Per-trial `startup_precision_review.json` accepts native atom-map agreement at
1e-6 Å; original stricter 1e-7 Å failures remain preserved. Startup mapping
acceptance is separate from force convergence and chemical acceptance.

All four portable movies have published under `<trial>/progress/latest/vmd/`.
Native engines and individual coherent exporters persist in recorded tmux
sessions; comparison session `nyl12_neb_spacing_comparison_20261119` writes
`spacing_tests/plots/comparison.png` and source-hashed JSON/CSV/checkpoints.
Compare Max/RMS gradient-to-target ratios, reset counts/timing, coherent
chemical-coordinate gaps and elapsed work at 50/100/250 cycles where budgets
allow. With frozen endpoints, 24/48-image bands require 22/46 interior native
evaluations per full cycle; compare actual `Calculating the energy` counts,
not cycle counts or maxima alone. No validated barrier or improved convergence
has yet been established. Curvature and downhill connectivity remain later gates.

Authoritative current details: `04_neb/spacing_tests/README.md`, per-trial
configs/reviews/job records, and the latest section of
`docs/standalone_handoff_20261020.md`. Earlier queue entries below retain the
execution history; they do not supersede this active comparison.

### Running context — original ASP272/ASH274 NEB (2026-11-19)

The user authorized NEB. Current canonical stage is
`09_water_relay_acylation/04_neb/` within the ASH274 campaign. Native session
`nyl12_ash274_acyl_neb_20261119` and coherent-band export session
`nyl12_ash274_neb_progress_20261119` are launched. The24regular-image band is
seeded by13controls from38optimized stage endpoints plus12window025 trace
checkpoints. Exact unbiased endpoints/Hlabels/native placement were checked;
no controls dropped. Earlier draft seed-order error was rejected before engine
launch and retained under preparation_history. QM102/−1/54realH/8links, zero
Kelvin electronic temperature and4,662mobile atoms in12Å remain unchanged.
Maximum16000cycles implies112000energy evaluations in this installed adapter;
no climbing image yet. There are no extra reaction/positional restraints.

Watch native process/scientific records and `04_neb/progress/latest/`, not
unusable native numbered restarts. Curvature/downhill verification remains
separate. VMD's PRO1312display error is fixed in `progress/latest/`; original
coordinates unchanged. Laptop pull includes both latest bundles and NEB mdin.


### Overnight outcome — product retained, reverse scan complete (2026-11-19)

The unbiased acyl-enzyme candidate force-converged and retained the ester and
balanced proton inventory. All36reversewindows and a final unbiased E–S
release force-converged. The product is +11.06 kcal/mol relative to returned
E–S; window 020 is the sampled physical maximum at +50.89 kcal/mol (bias
0.000108kcal/mol). Neither curvature nor a continuous MEP is verified.
Window24→25 drops18.59kcal/mol with a proton-occupancy change; inspect/refine
that interval before treating the profile as a continuous reaction path.
Thr Nalpha remains NH2 and the tracked transferring H is original waterH20318;
originalHG1 remains on water, so no Nalpha-water relay is established.
**Next bounded work:** retain the product endpoint, inspect/refine that jump,
and test window 020 as a saddle candidate with curvature and downhill checks.
No endpoint/reverse-scan jobs remain active; the original NEB and the four
parallel numerical trials above are separate running work. Source/summary:
`09_water_relay_acylation/overnight_summary_20261119.json` relative to ASH274;
final inspection under the same campaign's `progress/latest/`.


### Authorized execution — balanced acyl endpoint before another path

The user authorized the revised sequence on 2026-11-18. The bounded screen of
all six existing 100 ns replicas is complete: 60,000 correlated frames include
18,958 Tier A/B frames with positive signed phi, but only 51 also satisfy the
unchanged native-Thr gates. None contains the specified preorganized one-water
Nalpha–water–leaving-N bridge. Counts by microstate are in
`09_water_relay_acylation/00_screen/bounded_screen_counts.csv`, relative to the
ASH274 campaign. This excludes a ready bridge in that screened subset; it does
not exclude water reorganization, multiple waters or other conformations.
No MD snapshot or alternative microstate was silently promoted.

The independent endpoint test uses the **already accepted water-QM E–S source**,
not a claimed positive result from that screen. It preserves QM102/54 real H,
charge −1, eight link H, the original WAT1576 oxygen, periodic DFTB3/3ob at
zero electronic temperature, and the complete-residue 12 Å/4,662-atom movable
selection. There are no added positional restraints; the exterior remains
frozen as in the accepted calculation. Product construction conserves every
atom/proton and the ASP272/ASH274 inventory. The source E–S MM topology is
**not** used as a classical force field for the acyl-enzyme: reacting atoms are
QM throughout. Repeating E–S MM minimization is unnecessary for this accepted
source.

The first ideal planar-ester construction failed contact checks: imposing a
120° ester angle rotates an entire polymer arm and creates severe overlaps.
The source approach angle OG1–C–acyl-Cadj is 83.97°. The bounded 285-position leaving-fragment search plus a 24-value acyl torsion
sweep yielded zero severe contacts, with maximum heavy/H moves 2.172/2.366 Å.
The selected +15° torsion about bonded C19920–C19919 retains the acute source-like
83.97° ester angle for subsequent optimization. This is a coordinate hypothesis,
not an optimized endpoint or a reaction path. TIP3P's topological H–H constraint is omitted only
from the builder's chemical graph; the topology itself is unchanged. Failed
construction artifacts and the contact diagnostic are retained.

Once a seed passes the contact/inventory checks, the scripted tmux sequence is:

1. Product-seed QM/MM preoptimization with explicit temporary harmonic terms
   for ester formation, amide separation and leaving-amine proton occupancy.
2. Remove every driving term and optimize for up to 2,500 accepted steps.
   Require the unchanged force gate and balanced acyl-enzyme/amine retention;
   report every real-QM-H owner and all spectator proton inventories.
3. Only after that gate, run 36 reverse windows toward E–S with the single
   coordinate `r(OG1,C) − r(C,Nleave) + 1 Å*[CN(OG1,H) − CN(Nleave,H)]`.
   All 54 real QM H contribute; attack, cleavage and proton motion can vary
   independently orthogonal to this scalar. Nalpha occupancy is unrestrained.
   No stable TI1 is required. This scalar can still be chemically degenerate,
   so a coordinate value or energy maximum alone never establishes success.
4. Release the reverse endpoint without bias, audit connectivity and inventory,
   then nominate a chemically valid saddle for separate curvature/downhill
   verification. Production PMFs and redesign remain downstream gates.

Entry points are `scripts/prepare_relay_acyl_campaign.py` and
`scripts/run_relay_acyl_campaign.py`; one logical campaign is
`qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/`.
The campaign `protocol.json`, stage records and `campaign_status.json` will
distinguish running, force convergence, chemical retention and process failure.
The live accepted-state movie/plots use `progress/latest/`, and the laptop pull
includes compact inspection artifacts. Verify launcher records before calling
a queued/prepared calculation running. The native workflow is launched in tmux `nyl12_acyl_endpoint_20261118`;
the accepted-state monitor is `nyl12_acyl_progress_20261118`. Inspect
`02_acyl/job/` and stage scientific records for current progress.

### Controlling revision — productive hydrolysis, not a required TI1 minimum

The user's mechanistic review request motivates this revision (2026-11-18).
**Success is a chemically balanced, kinetically credible product-connected
pathway and verified relevant transition state(s). A stable TI1 is optional.**
Literature supports both stepwise Ntn chemistry (Andjelkovic 2023, R41) and
amide acylation without a metastable tetrahedral intermediate in the calculated
path (Palermo 2015 FAAH, R48). “Concerted addition” can still yield a TI, so
identify the actual coupled bond changes (Nutho 2019, R46). Negoro's mature
NylC Fig. 10 is a BapA-guided proposal; precursor autocleavage evidence must
remain separate. Nitrogen inversion/proton access are explicit design concerns
(R44,R45), not merely oxyanion geometry.

Both outstanding calculations completed normally. The unbiased D release
force-converged back to E–S after 249 steps: OG1–C 3.211 Å, Thr-NH2/OG1H,
with Thr/water proton-label exchange. DFT and xTB embedded-QM gradients at D
both oppose additional attack while retaining the specified proton occupancy.
This is a negative result for this configuration/branch, not proof that no TI1
or catalytically viable pathway exists. D was still only partially added before
release. No new barrier, minimum or broadly validated Hamiltonian follows.

The missing chemical test is **leaving-group protonation plus C–N scission**.
Holding q_occ near 2.10 favors Thr-NH3+ and penalizes the normal net donation
of its extra proton onward unless another proton replenishes it. The pilot
therefore was not a complete test of an N-terminal general-base/general-acid
cycle. Coordinate audit also finds Nalpha–leaving-N 5.80 Å in E–S and 5.23 Å
at D; direct donation is not arranged. At D, QM-water O–leaving-N is 3.35 Å
and its nearest H is 2.60 Å. A reorganized water relay is a concrete hypothesis,
not an established path. Asp272 oxygens and the ASH274 proton are too distant
for direct donation in these snapshots. Nyl12 Phe113 also lacks the phenolic
OH of the corresponding proposed NylC residue Tyr146; its backbone NH remains.
Do not silently add Tyr chemistry or neutralize Lys156 to obtain a desired path.

**Bounded next sequence (proposed; no new jobs launched by this review):**

1. Screen existing equilibrated structures for a credible complete proton
   network to the **leaving amide N**, retaining useful carbonyl/oxyanion
   geometry. Prioritize an Nalpha–water–leaving-N relay if geometry supports it;
   inspect Asp-assisted alternatives only when they supply a reachable relay.
   This is distinct from requiring the earlier fixed-pocket OG1–Nalpha ring.
   Keep at most one primary route and one justified control.
2. Prepare and relax a chemically balanced acyl-enzyme/amine endpoint in the
   same atom mapping and explicitly declared microstate. Thr covalent catalysis
   is the best family-guided hypothesis, not a consequence of net hydrolysis
   stoichiometry; direct activated-water attack remains a lower-priority
   alternative if the covalent route is not supported.
3. Connect E–S to that endpoint with independently evolving attack, C–N
   cleavage and proton-transfer coordinates. Permit the catalytic Nalpha and
   water occupancies to change; do not lock in a TI1 protonation pattern or
   demand a minimum at a prescribed window. Allow nitrogen pyramidalization,
   inversion and polymer reorganization. Classify any TI only if it emerges.
4. Verify the chemically relevant saddle and downhill endpoints, then assess
   deacylation before advancing multistate design or production PMFs. A
   persistent TI1 is not required for LigandMPNN/redesign; a credible TS,
   reactant and product/intermediate ensemble is required.

Do not reopen another TI1-only grid or repeat Hamiltonian swaps on the same
geometrically incomplete branch. Percent polymer conversion in crude lysate
cannot determine kcat or a unique chemical barrier: adsorption/accessibility,
productive-state population, enzyme amount and turnover steps are confounded.
All completed negative results remain evidence. See manuscript R44–R48 and
`analysis/nyl12_mechanism_evidence_20261118/` for literature/coordinate provenance.

### Previous decision — complete pilot, test D release and electronic response

All four points are force-converged. The proton-label-invariant coordinate
retains the intended N3/OG0 state without water proton exchange; nevertheless
attack reaches only 1.854 Å and costs +42.763 kcal/mol relative to original
E–S. Carbonyl C–O is 1.292 Å, so D fails the established ≤1.75 Å attack and
≥1.30 Å tetrahedral C–O screens. There is no accepted TI1 or validated barrier.

The endpoint audit also identifies loss of one oxyanion contact: substrate O
to Phe113 N grows 3.161→4.110 Å, while Asn186 ND2–O is 2.837 Å. This motivates
testing both geometry persistence and the local electronic response rather
than increasing the bias or starting another scan. Evidence:
`occupancy_callback_pilot/endpoint_audit/{audit.json,summary.csv,supplement.json}`.

The user-authorized next tests are launched:

1. Remove **both** driving biases from D; optimize in the unchanged 12 Å
   active region with native gradient L-BFGS and the same force gate. Session
   `occupancy_D_unbiased_release`; output `pilot/point_D/unbiased_release/`.
   Audit all proton ownership and original bonds on completion. A returned
   E–S endpoint is a negative result for this starting geometry/microstate,
   not proof that the enzyme lacks a viable hydrolysis mechanism.
2. One exact-D ωB97X-D/def2-SVP/GFN2-xTB force comparison in the existing
   finite embedding, under `d_reference/` (session `occupancy_D_method_reference`).
   Inspect projected attack/proton-transfer and oxyanion-contact derivatives.
   Do not compare embedded-QM-only slopes directly to native total QM/MM
   slopes or promote a Hamiltonian on a single result.

**Next decision:** review both results before another path/sampling campaign.
TI1-like force convergence would nominate a stationary-point candidate;
minimum/TS verification requires separate curvature/connectivity evidence.

### Current repair — native gradient-based continuation, 2026-11-18

The line-search failure is explained quantitatively by ordinary MM pair
11272/11330 crossing the 10 Å nonbonded cutoff. The resulting ≈0.001 kcal/mol
energy jump exceeds the predicted smooth energy change by about 32,000-fold;
pair energy and force predictions match within native print precision. This
is not an SCC failure, a changed QM membership or evidence against activation.
`occupancy_callback_pilot/line_search_diagnosis/` retains reproducible evidence.

The user-authorized fix resumes from the final accepted A snapshot using the
installed DL-FIND gradient-based L-BFGS routine, with 0.05 Å maximum per-atom
steps and explicit history resets on nonpositive curvature. No dependency,
potential, protonation, spring, QM/active selection or convergence criterion
changed. Session: `occupancy_attack_pilot_continuation`; script:
`occupancy_callback_pilot/resume_pilot.py`. The original failed terminal records
are preserved in named history folders; native calls and accepted frames append
in place. Points B–D still require A force convergence. The same diagnostic
scope applies: no automatic unbiased release or claim of TI1/TS/PMF.

Repair verification succeeded: point A reached maximum biased gradient
0.109355 kcal/mol/Å in three additional steps (iteration 161), retaining
three Thr-N protons, deprotonated OG1 and intact QM water/ASH274. Point B is
now active. The accepted-step plot/movie monitor was restarted against the
same progress folder. This establishes successful numerical continuation,
not stability of the biased activated geometry after release.

### Latest execution checkpoint — point A incomplete, 2026-11-18

The four-point pilot stopped at A after 158 optimizer iterations/284 native
evaluations, with L-BFGS-B status 2 (`ABNORMAL` line-search termination).
Native/SCC checks passed, but maximum active biased gradient 0.277869
kcal/mol/Å exceeded the 0.12 criterion. Points B–D remain unstarted.
`pilot/completion.json` is `incomplete`; launcher exit 0 is process completion,
not scientific acceptance. Retain this attempt. Diagnose the final line-search
energy/gradient behavior before a recorded continuation from its final accepted
snapshot; do not loosen the tolerance or change chemistry silently.

### Current decision — verified callback and four-point pilot, 2026-11-18

The analytic occupancy bias now reaches the actual optimizer through
`ActiveBiasCallback`, using native Sander single-geometry energies and analytic
gradients. The two-iteration functional proof passed all seven gates: a known
0.20 kcal/mol energy addition, nonzero added force, correct gradient composition,
native energy conversion, optimizer execution, and requested geometry mapping.
It used four native evaluations; its iteration-limit exit was intentional and
does not establish stationary-point convergence. Evidence:
`08_activation_attack_grid/occupancy_callback_pilot/callback_proof/`.
Native NFE/PMD remains rejected for this DL-FIND minimization backend.

The authorized four-point diagnostic is launched in tmux
`occupancy_attack_pilot`, under the same grid's `occupancy_callback_pilot/pilot/`.
Starting at q_004_attack_001, it targets q_occ=2.10 independently of OG1–C
targets 3.1251784, 2.60, 2.10 and 1.60 Å (points A–D, serial). Harmonic
potentials are 0.5*kq*(q−q0)^2 + 0.5*kr*(r−r0)^2, with kq=100
kcal/mol per coordination-unit squared and kr=200 kcal/mol/Å². These smooth
occupancies are not formal protonation labels. Both component counts, all
54 real QM H owners, original QM bonds and attack distance must be inspected.
QM102/charge −1, 4,662 movable atoms, periodic electrostatics and 0 K electronic
temperature are retained. There are no additional positional restraints or
bounds. Each point permits 1,500 L-BFGS-B iterations/3,000 function calls and
must reach maximum active biased-gradient component ≤0.12 kcal/mol/Å before
continuation. Process failure or an unconverged point stops the pilot explicitly.
This is a bounded diagnostic, not a converged reaction path, PMF or barrier.
Inspection products are `occupancy_callback_pilot/pilot/progress/accepted_progress.png`
and `progress/vmd/view_path.vmd`, refreshed from accepted steps and selected
accepted restart snapshots. Plot ΔE is relative to the first pilot physical
energy, not the original E–S reference. The laptop pull preserves these paths
and excludes full native-evaluation artifacts.

Saved reference-gradient analysis is complete, without new QM evaluations:
`reference_checks/gradient_analysis/`. At the three sampled states, DFT and
GFN2-xTB embedded-QM gradients both oppose shortening OG1–C at fixed q_occ;
the strongest opposing slopes are about 60 and 55 kcal/mol/Å. Some small
proton-transfer slopes have opposite signs. GFN2-xTB therefore is not promoted
as globally validated. Native DFTB3 *total* QM/MM gradients are recorded
separately and cannot be compared directly to embedded-QM-only gradients.

**Next gate:** inspect the pilot's actual occupancy/connectivity and force
acceptance. Only a chemically suitable force-converged candidate warrants a
separately recorded unbiased release or further matched method checks. Do not
restart the labelled-H grid or umbrella sampling from a nominal coordinate
value alone. Earlier checkpoints below remain execution history.

### Previous checkpoint — occupancy audit and unbiased release, 2026-11-18

The authorized audit of all 25 points is complete: 54 real QM hydrogens per
point (1,350 assignments; fictitious link H excluded), with nearest-heavy
ownership, ambiguity thresholds and smooth coordination counts. Rows 1–3
keep HG1 nearest OG1, row 4 places it nearest Thr N, and row 5 exchanges it
with a water proton. No detected hydrogen escaped the real QM heavy-atom
region. Threshold-based assignments do not establish electronic bond orders
or protonation free energies. Evidence is in `08_activation_attack_grid/proton_occupancy_audit/`.

The separately authorized one-candidate diagnostic release completed with
normal engine/SCC/force convergence in
`08_activation_attack_grid/diagnostic_activation_release/`. It started from
q_004_attack_001 and removed both reaction restraints while retaining the
QM102/−1 setup and 4662-atom movable shell. It did **not** retain activation:
HG1 moved to water, water H1 moved to OG1, and Thr N retained only its original
two hydrogens. The endpoint is therefore E–S chemistry, 4.7757 kcal/mol lower
in native potential energy than the original E–S reference. Its geometry is
Tier B (d1/d2/d3 3.205/3.103/2.968 Å; attack angle 117.39°; signed phi
+102.52°). Do not promote it automatically over the original Tier A seed or
call it a stable activated intermediate. The helper's ligand_rmsd=0 field
is not a measured displacement from the original E–S reference. A 13-frame
stripped movie retaining WAT1576 is in `diagnostic_activation_release/vmd/`.

The proton-label-invariant definition is
`q_occ = CN(Thr N,H) − CN(Thr OG1,H)`, where each count is
`sum_H 1/(1+(r/1.25 Å)^6)` over the 54 real QM hydrogens. Both separate
counts must also be monitored: their difference alone is not a complete
chemical-state assignment. The reusable analytic value/gradient implementation
is `src/nyl12_covalent/proton_coordination.py`; 12 focused tests passed with
100% line coverage, including hydrogen permutation invariance and a custom
derivative check on synthetic coordinates. It agrees with the independent
audit for all 25 points to 1.2×10⁻¹⁵. These smooth counts are not integer
valences, formal charges or pKas.

**Backend preflight failed:** the installed Sander/DL-FIND minimization did
not initialize native NFE/PMD occupancy restraints. A deliberately nonzero
0.20 kcal/mol bias gave zero energy contribution, with no PMD/CV output.
See `occupancy_coordinate_preflight/attempt_2/scientific_acceptance_corrected.json`;
the original 0.40 expectation is retained with its half-k correction. No
Amber installation was changed. Before another occupancy-driven path, wire
the tested analytic bias into a supported energy/gradient callback and prove
that it contributes to the optimizer's energy and forces, or explicitly repair
the native initialization with equivalent verification. Merely adding native
`&pmd` inputs to current DL-FIND minimizations is invalid. No new driven path
or sampling campaign has launched.

The three-state finite-embedding reference calculations also finished. Relative
energies (E–S zero) for the partially activated/attack structures are
18.56/39.23 kcal/mol with DFTB3, 13.84/33.45 with ωB97X-D/def2-SVP, and
19.95/27.20 with GFN2-xTB. These energies alone do not validate a fast
Hamiltonian or locate a minimum; reaction-direction gradient interpretation
remains pending. Preserve the two current priorities: correct coordinate
integration, and interpretation of the saved reference gradients. Do not
restart the labelled-H grid or umbrella sampling.

### Execution result — completed grid, 2026-11-18

All 25 reported points and all five row sweeps completed. No point passed the
declared TI1 geometry gate. For rows 1–4, the last attack points reached
OG1–C = 1.887, 1.886, 1.860 and 1.783 Å, with physical relative energies
+28.34, +32.71, +36.62 and +54.58 kcal/mol. Energy increased along attack
in every row; these are constrained PES observations, not barriers or proof
that no TI1 basin exists. The 1.50 Å attack target was not achieved.

The largest-q row is a **proton-label exchange branch**: Thr HG1 binds to
WAT1576 O (0.986 Å in its first grid point), while water H1 binds Thr OG1
(0.980 Å). Thr remains an alcohol and water remains water; Thr N is not
protonated by this exchange. Its −4.54 kcal/mol first point is a rearranged
reactant geometry, not a favorable activated intermediate. The labelled-H
coordinate therefore did not sample complete intended N-terminal activation
over its full nominal range. Any follow-up must distinguish chemical proton
occupancy from atom labels, e.g. through smooth O–H/N–H coordination variables;
do not merely extend this degenerate coordinate or increase its spring.

The E–S three-method reference finished normally, but selected-state exports
stopped because Amber truncates its Gaussian executable field to 20
characters. A real Sander export preflight verified the repair: a short local
wrapper plus the required Gaussian template. Session
`ash274_reference_continuation` now runs only the two unfinished activated/
attack references concurrently, reusing completed E–S and DFTB3 results.
Original failures remain under each state's
`extern_export_failed_pathlength/`. No TI1-like state was available for the
third selection. No unbiased release has run; the reference-review gate is
still pending. The next decision is to inspect the completed matched relative
energies/gradients and then select whether either actual candidate warrants
release. Do not launch a replacement pathway campaign before that decision.

### Controlling decision — after completed water-relay scan, 2026-11-17

**Current state: no active calculation in this path; no validated TI1 or
first-order acylation TS.** The 25-window water-to-Asp272 scan completed with
force-converged windows and continuous monitored reactive distances. Its final
physical potential energy is +54.0235 kcal/mol relative to its own E–S seed
(bias-inclusive +58.1547; bias +4.1312). OG1–C reached 1.852 Å despite a
1.650 Å restraint target. The endpoint is `endpoint_negative_or_incomplete`;
the conditional unbiased release did not run. Both intended proton transfers
occurred: Thr HG1 moved to water and water H2 moved to Asp272. ASH274 retained
HD2, leaving **both Asp residues protonated and Thr N neutral**. This particular
charge redistribution is energetically unfavorable along the sampled scan.
It does not establish the absence of a tetrahedral basin elsewhere on the PES.
Evidence: `qmmm/pa6_flipped_model30/07_ash274_comparison/07_water_to_asp272_ti1_like/`
(`protocol.json`, `window_025/acceptance.json`, and `energy_ledger.csv`).

**Immediate objective:** discriminate an inadequate prescribed path from
absence of a persistent tetrahedral intermediate in this microstate. A stable
TI1 is a hypothesis, not a prerequisite imposed on the mechanism. Experimental
hydrolysis establishes overall chemistry but does not identify a microscopic
route, intermediate lifetime, or the correctness of the modeled binding pose.

1. **One bounded alternative on the same accepted E–S seed.** Prepare a coarse
   5 × 5 relaxed PES map with independent coordinates
   `q_N = d(Thr929 OG1,HG1) − d(Thr929 N,HG1)` and
   `r_attack = d(Thr929 OG1,substrate carbonyl C)`.
   Suggested initial ranges are the measured E–S q_N (~−1.78 Å) to +1.05 Å,
   and r_attack 3.10 to 1.50 Å. Freeze the exact grid and restraints in a
   reviewed protocol before execution. This tests direct N-terminal activation
   and different timings of attack without tying the two targets together.
   Retain the same 102 QM atoms, −1 total QM charge, water, ASH274/ASP272
   starting microstate, topology, and 12 Å movable shell. Treat any unexpected
   proton relocation as a distinct chemical branch requiring review. Preserve
   full proton accounting; water and the Asp pair may reorganize. Use the
   grid to find candidate valleys, not a barrier or a proof of a minimum.
2. **Validate only informative points.** At most four representative E–S,
   activated, attack, and tetrahedral-like geometries receive matched-embedding
   DFTB3/GFN2-xTB/ωB97X-D checks, reusing previous results where identical.
   Inspect reaction-direction forces as well as relative energies. Keep the
   point-charge membership and boundary treatment consistent in cross-method
   comparisons. The existing benchmark found DFTB3's relative TI1-like error
   +7.42 kcal/mol versus DFT and xTB's direct-PT error +8.02; neither is a
   validated universal replacement. Stop fast-method production if the new
   reaction-region comparison disagrees qualitatively with the reference.
3. **Release candidates, then connect endpoints.** Release reaction restraints
   on at most two plausible candidate basins. Require force convergence and
   retained chemistry. If retained, refine the valley with a string/NEB seeded
   from the map and then verify its saddle(s). If candidates collapse, one
   bounded restrained QM/MM equilibration followed by independent releases may
   test environmental reorganization; declare sampling and restraints in
   advance. A constrained basin alone is not a stable TI1.
4. **Explicit decision if TI1 remains absent.** Seek an E–S → acylenzyme path
   allowing attack, leaving-group protonation, and C–N cleavage to couple;
   do not insist that it pause in TI1 or reuse only the strained direct
   four-center route. Use the observed proton network to select the relay.
   If the grid/reference checks still leave no credible route, make one
   targeted QM-region/protonation sensitivity test (e.g. the documented
   Lys156/Gly273 expansion), or obtain external QM/MM review of the setup.
   Avoid another broad collection of paths with untested assumptions.

**Work limits:** this is the next proposed calculation sequence; updating this
document does not mean these calculations have launched. Keep the old direct
route and its unresolved Hessian index as a comparator, defer further costly
full-shell eigensolver attempts and umbrella sampling until a credible path is
selected, and preserve all prior negative results. Do not shrink the movable
shell. Refine the grid only where a valley or discontinuity gives a specific
reason. A second starting conformation is a robustness check after the first
mechanism passes, not another simultaneous campaign. Success means connected,
chemically assigned endpoints, stationary gradients, a reproducible index-one
TS and downhill connectivity, followed by selected DFT checks. PMFs follow
that gate. No computational plan guarantees those results for the present
model; each step must resolve a named uncertainty.

The historical updates below remain evidence. This dated decision supersedes
their outstanding next-run instructions where they conflict.

**Execution status — 2026-11-18:** The bounded activation/attack grid is
**RUNNING** as five independent tmux row workers. Its 25 reported points use
135 bounded optimization stages. The same QM102/−1 model, 4662-atom movable
shell and periodic Hamiltonian are retained. The E–S ωB97X-D/def2-SVP force
reference queue is now running with up to two eight-core DFT workers; it uses
a separately declared nonperiodic full-PC sensitivity model. Gaussian runs
directly on the exported 110-atom QM/link geometry and 104,824 adjusted charges
to avoid its field-query centre limit. The earlier Sander/EXTERN exit-0 result
failed Gaussian's termination/gradient gates and is rejected. The progress
monitor refreshes a measured-coordinate energy plot and per-row VMD bundles.
The candidate-release watcher requires explicit comparison review before
running at most two releases. Launcher records, frozen inputs and logs are under
`08_activation_attack_grid/job_row_N/` and `reference_checks/`; running status
does not establish force convergence or a retained intermediate.

**Update — corrected-shell ASH274 E–S accepted, 2026-11-17.** The matched
ASP272/ASH274 seed now has a force-converged, accepted Tier A DFTB3/MM
endpoint in `qmmm/pa6_flipped_model30/07_ash274_comparison/06_qmmm_optimization/replica_1_tier_A_seed/corrected_shell/`.
It retains ASH274 HD2 and neutral Thr233 N, with d1/d2/d3 =
3.07665/3.06706/2.94462 Å, attack angle 116.807°, and signed phi +97.925°.
The corrected movable shell contains 4,650 atoms; water 1576 is now movable
and its oxygen shifted 0.1465 Å. This corrects the inherited shell that froze
111 residues within 12 Å; the earlier accepted optimization was an incomplete
relaxation and is superseded. The final audit
`corrected_shell/water_1576_full_hbond_audit_final.json` finds
Thr233 OG1–HG1···water O at 167.9° and H···O 1.811 Å, and water
O–H2···Asp272 OD1 at 173.4° and H···OD1 1.533 Å. It does not establish a
bridge to Thr233 N or substrate leaving N. This local geometry is not evidence
of water occupancy or reactivity. A one-water-in-QM (102 atoms, charge −1)
reoptimization completed in tmux `ash274_water_qm_reopt` in the sibling
`water_qm_reoptimization/`. It started from the corrected-shell converged
endpoint, preserved the −1 QM charge and eight link cuts, and expanded the
12 Å movable region to 4,662 atoms to include four boundary waters. Sander
exited normally and DL-FIND converged at cycle 248. The endpoint retained
Tier A geometry (d1/d2/d3 = 3.100/3.063/2.959 Å; attack angle 116.44°;
signed phi +98.75°) and the declared Thr, water and Asp274 covalent proton
inventory. Thr HG1···water O is 1.932 Å at 165.5°, and water H2···Asp272
OD1 is 1.730 Å at 172.7°; see `scientific_water_assessment.json`. There is
still no bridge to the substrate leaving N (water O···N = 4.107 Å). This is
a stable *local geometry*, not evidence of water occupancy, proton transfer,
or a viable barrier. Next test one bounded water-assisted Thr activation and
OG1 attack route toward TI1, with explicit proton/geometry gates and a
separate unrestrained endpoint release. Compare with the direct route only
after identifying chemically continuous states; do not broaden the campaign.

1. **Asp-pair microstate comparison:** all six matched 100 ns, 300 K
   replicas completed. Tier A counts per 10,000 frames were 37/238 for
   ASP272/ASP274, 3886/1724 for ASH272/ASP274, and 1677/0 for
   ASP272/ASH274. Site-1 median closest Asp-pair O–O distances were
   4.89/4.90, 2.56/2.56, and 2.73/2.65 Å, respectively. The large
   between-replica variation, including zero Tier A frames in one ASH274
   replica, prevents ranking protonation microstates from frame counts
   alone. Next inspect ligand register, the 10 ns block series, and
   proton/water geometry before choosing at most one follow-up snapshot
   per microstate. Fixed-charge MD cannot exchange the Asp proton.
   A three-frame, protein-only PROPKA 3.5.1 screen is recorded in
   `docs/propka_asp_pair_pilot_20261116.md`; it does not select a
   protonation microstate. Assess ligand/water sensitivity before choosing
   the next QM/MM E–S model.
   The matched 10 ns block and water-oxygen analysis is in
   `docs/asp_pair_md_water_screen_20261116.md`. Water proximity is common
   in one ASH274 replica but the inspected candidate lacks a complete
   hydrogen-bond relay.
   That two-leg check has now completed: only 1, 6 and 32 Tier A frames
   in the three screened replicas pass both hydrogen bonds to Thr N,
   with at most three consecutive 10 ps samples for one water. None
   passes both legs to the leaving-group N while Tier A. Treat an
   ordered water relay as unproven and screen the remaining substrate
   register/protonation gates before choosing QM/MM seeds.
   Compare Asp306/Asp308 pKas from mature NylC 5XYG with **uncleaved**
   precursor controls 7YU0 and 7YU1, explicitly representing the new
   Thr267 N terminus only in the mature structure. These precursors carry
   maturation-blocking substitutions; differences cannot be assigned solely
   to cleavage. The first comparison completed with process exit 0, but a
   malformed `TER` record prevented recognition of mature Thr267's free
   N terminus. An atom-identical corrected input restores N+267 and removes
   its warning. Corrected PROPKA predictions for mature 5XYG are N+267
   6.35–6.55 (default) or 4.59–4.60 (coupled alternative), Asp306
   2.19–2.21, and Asp308 5.49–6.49 (default) or 7.24–8.45 (coupled
   alternative). The earlier mature pKas must not be used. The corrected
   audit is `analysis/negoro_2023_protonation/terminal_pka_20261116/before_after.json`;
   the corrected scripted run passed end to end under
   `analysis/negoro_2023_protonation/crystal_propka_20261116/attempt3/`.
   the precursor comparison is under
   `analysis/negoro_2023_protonation/crystal_propka_20261116/attempt2/`.
   Precursor predictions vary by chain and coupled-state treatment. These
   are structure-dependent predictions, not measured protonation assignments.
   Compare standard and coupled-site PROPKA predictions, then use an
   isolated MCCE4 trial on a small, matched structure set if installation and
   input validation succeed. Neither method establishes protonation by itself.
   The isolated mature-5XYG MCCE4/DelPhi trial finished, but its pKa bounds
   (both Asp <0; Thr267 N terminus >14) are **scientifically rejected**:
   MCCE generated 13 charged termini at internal missing-coordinate
   boundaries. See `analysis/negoro_2023_protonation/mcce4_alpha_trial_20261116/scientific_acceptance.json`.
   Do not use these bounds to assign an Asp state or run the mutant precursor
   controls with the same input preparation. Next preflight one complete
   Nyl12 Tier A protein model with MCCE step 1 only and require exactly its
   true mature-chain termini; launch a full pKa calculation only if that
   structural gate passes. Independently compare one matched Tier A
   ASP272/ASH274 E–S seed with the previous ASP/ASP and ASH272/ASP274
   simulations using identical QM/MM atom mapping, charge and relaxation
   gates. This is a stability test, not an absolute pKa determination.
   The complete Nyl12 Tier A MCCE4 step-1 preflight also **failed**:
   it retained all eight mature protein chains but assigned zero NTR/CTR
   groups and emitted 2,840 `Error!` diagnostics. No full pKa run was
   launched. See `analysis/nyl12_mcce_step1_preflight_20261116/neutral_replica_02_frame284_attempt2/scientific_acceptance.json`.
   Stop MCCE4 troubleshooting for this mechanism campaign unless a tested
   structural-input fix becomes available. Use the matched microstate
   stability comparison and, if needed, a targeted free-energy calculation
   to estimate proton occupancy; static PROPKA remains a sensitivity screen.
   Defer constant-pH MD until these cheaper screens identify a specific,
   consequential Asp/Thr-N microstate ambiguity and a validated protocol can
   handle the tetramer, termini and nonstandard PA6 substrate.
2. **Direct-route saddle:** 0 K DFTB3/MM dimer refinement from the
   converged NEB maximum finished at cycle 99 without SCC failure and with
   native force convergence. Its geometry matches the earlier downhill-tested
   candidate within 0.006 Å RMSD over the movable region. The bounded
   block-LOBPCG probe finished, but its approximate eigenvalues
   (−0.227, +0.033 Eh/(Å² amu)) have residual norms 0.051/0.061,
   far above 10⁻⁴. The Hessian index is unresolved. Require exactly one
   reproducible negative mode and two-sided endpoint connectivity before
   a true-TS claim. Do not report the
   35.81 kcal/mol NEB image maximum as a validated activation barrier.
   A bounded operator check now shows exact repeatability, <1% step
   sensitivity and 0.37% bilinear asymmetry; see
   `docs/ts_hessian_operator_diagnosis_20261116.md`. Gross force noise
   is unlikely to explain the stalled solver. Audit preconditioning or
   focused mode seeds before another full-region eigenvalue search. The next
   bounded test should recycle approximate Ritz vectors, improve conditioning
   without shrinking the 12 Å movable region, and compare residual-versus-
   gradient-evaluation curves against the previous probe. Accept a first-order
   saddle only after a reproducible single negative mode with converged
   residual and two-sided downhill connectivity from this same TS structure.
   A first warm-start pilot reused the two saved **partial** Ritz vectors and
   completed 38 HVPs/76 native gradients. Residuals improved to
   4.13×10⁻⁴ (negative candidate) and 2.20×10⁻² (next candidate), still
   above the 10⁻⁴ gate. Thus the Hessian index remains unverified; a
   better-conditioned follow-up is needed before a true-TS claim. The pilot,
   failures and acceptance record are under the direct-route TS tree's
   `mode_probe_block_warmstart_shortpath_20261116/` directory.
   A six-vector lowest-algebraic pilot then exhausted its 60-HVP budget
   before returning Ritz pairs; process exit 0 does not imply scientific
   acceptance. Its inputs and failure record are under
   `mode_probe_block6_warmstart_20261116/`. Require the eigensolver to save
   intermediate Ritz/residual checkpoints at its budget cap, then run one
   budgeted wider-block continuation. If the second low mode still fails
   to converge, diagnose a preconditioner or alternative eigensolver before
   spending more gradients; do not call the current structure a verified TS.
   The final six-vector run returned one negative and five near-zero positive
   *estimates*, but all full-Hessian residuals were 0.0076–0.0156 versus
   the 10⁻⁴ gate; index one is unresolved. Test numerical restart/damping or
   a positive damped preconditioner on synthetic clustered modes first,
   then allow one bounded native pilot if it actually improves convergence.
   A scalar Hessian shift cannot untangle eigenvectors because it preserves
   their spacings; adding physical restraints would change the PES and
   cannot certify the unrestricted saddle. Prioritize the triad-assisted
   mechanism if this high-barrier direct route remains difficult to verify.
   A 12-gradient precision gate on the unchanged system found exact repeated
   HVPs, but 0.025 vs 0.0125 Å√amu products differed by 8.44×10⁻⁴ and
   1.79×10⁻³ Eh/(Å² amu), and bilinear asymmetry worsened at the smaller
   step. The nominal 10⁻⁴ residual target is not numerically supported by
   this central-difference operator. Do not relax the gate arbitrarily;
   test a higher-order finite-difference HVP or analytic/mixed Hessian method
   on a bounded set of directions and calibrate the tolerance before another
   full eigensolver run. Evidence is in the direct-route TS tree's
   `hvp_precision_gate_20261116/attempt2/`.

Full commands, hashes and exact paths are at the top of
`docs/standalone_handoff_20261020.md`. Do not open a new mechanism campaign
until these gates are assessed.

## Question and boundaries

Determine a chemically credible, reproducible route for WT Nyl12 cleavage of
PA6 L4 at the L2-producing bond and obtain stationary E–S, intermediate,
transition-state and acylenzyme structures useful for multi-state redesign.
Keep the tetramer and explicit solvent. Use the flipped, positive signed
Herschlag-phi pose as the first test, while retaining the opposite pose as a
bounded alternative. Treat substrate hydrolysis separately from NylC/Nyl12
autocleavage. Do not expand to PA66, poly(MAXA), mutants or a full polymer
binding free energy in this mechanism campaign without a new decision record.

The working chemical comparison keeps the Thr233 N terminus neutral and tests
ASP272/ASP274, ASH272/ASP274, and ASP272/ASH274. The Negoro NylC protonation
proposal is an inference, not established Nyl12 protonation. Test a direct four-membered
OG1-HG1-to-leaving-N route against one triad/water-assisted route. Add another
route only when geometry, conserved residues, or a failed chemical gate gives
specific evidence. The experimental assay at 75 °C and uncertain Tris pH at
that temperature cannot be equated to the 300 K MD elementary rate.
The initial neutral Thr233 N terminus is also a hypothesis. If it is NH3+,
it lacks a lone pair and cannot accept Thr HG1 as a general base in that
microstate. Then test a carboxylate- or water-assisted Thr activation, or
transient N-terminal deprotonation, with an explicit proton/charge inventory.
An NH3+ group might instead donate a proton to the departing substrate N
during tetrahedral-state collapse; favorable electrostatics alone do not
establish this catalytic role. Do not carry a neutral-N reaction coordinate
unchanged into an NH3+ calculation.

## Current evidence and what it does not prove

* Two ASH272/neutral-N 100 ns replicas completed. The geometry screen found
  3,886 and 1,724 Tier A frames versus 37 and 238 in the earlier ASP272
  replicas. The three campaigns share a static model30 protein/substrate pose
  and 300 K protocol, but solvent placements and random seeds differ. This
  is a promising association, not a causal protonation or pKa result.
* The first direct-route NEB stopped unconverged at 500 cycles with a
  ~70.6 kcal/mol outlier. A later 24-regular-image, trace-seeded NEB converged
  at cycle 477 and has a 35.81 kcal/mol image maximum relative to its first
  image. Its endpoint identities and local image continuity pass; a full
  first-order saddle is unverified. The older dimer point is ~46.4 kcal/mol
  above a separately minimized E–S reference, so these are not comparable
  activation barriers or PMFs.
* All four ASH272 charge-audited DFTB3/MM E–S optimizations converged with
  99 real QM atoms, eight links and charge −1. All four are geometric Tier A,
  but each transfers the prepared Asp272 proton to Asp274. Treat this as a
  microstate observation under the tested Hamiltonian, not proof of in-vivo
  protonation. Do not label their endpoints ASH272-protonated E–S minima.

## Execution gates, in order

1. **Resolve the E–S microstate comparison.** The six matched 100 ns
   replicas, Tier A/10 ns block analysis, and initial pKa and water screens
   are complete. Confirm the PA6 register and Thr orientation in candidate
   frames, test whether including PA6 or alternative pKa methods change
   the Asp/Thr interpretation, and select at most one accepted
   representative per microstate for follow-up QM/MM. Do not multiply
   reaction paths by six or assign protonation from static PROPKA alone.
2. **Validate electronic setup.** Parameterize the QM/MM input by declared
   protonation and charge; fail on mismatch. Keep the predeclared 12 Å movable
   region initially. Check charge, link cuts, periodic embedding, explicit
   waters, and MM-only versus QM/MM endpoint geometry. Stop a candidate if it
   loses the PA6 register or catalytic geometry during unbiased relaxation.
3. **Choose two mechanisms, then screen cheaply.** From a chemically accepted
   E–S microstate, build chemically explicit acylenzyme and possible tetrahedral
   states; equilibrate the surrounding MM degrees of freedom if an endpoint
   collapses, with the QM state restrained and all restraints recorded. First
   test triad/water-assisted proton transfer and attack. Retain the direct
   route as a comparator until its TS and energy are verified. Use restrained driving
   only to construct seeds; it is not an energy barrier.
4. **Converge a path before TS refinement.** Optimize both endpoints without
   reaction-coordinate restraints. Run regular NEB/string with smooth,
   geometry-checked images and enough cycles to meet predeclared force and
   displacement tolerances; only then use climbing-image refinement. Keep
   image-wise energies, forming/breaking distances, proton inventory and VMD
   frames updated. If the band stalls, diagnose its highest-force image and
   missing coordinate, then change **one** element of the path setup. Do not
   equate a 500-cycle stop or shoulder with convergence or a new minimum.
5. **Verify stationary points.** Refine each candidate TS with DL-FIND or an
   appropriate mode-following method at the screening Hamiltonian. Require
   stationary gradient, exactly one relevant negative Hessian mode, and
   downhill/IRC connectivity to the intended minima. The 4,165-atom numerical
   Hessian eigensolver is currently inefficient; diagnose precision and
   eigensolver strategy before merely raising its evaluation cap. Perform
   sparse hybrid-DFT/MM checks of key geometries and relative energies. If a
   fast Hamiltonian changes state ordering or proton identity relative to
   DFT, do not promote its PMF to mechanistic evidence.
6. **Extend to deacylation only after acylation is credible.** Add water attack,
   the second tetrahedral state, product-forming transition state(s) and final
   products. For redesign, use a consistent ensemble of E–S, acylation TS,
   acylenzyme and deacylation TS when available; an early design pilot may use
   validated acylation states but must be labeled incomplete-cycle design.

## Free-energy decision

Umbrella sampling is a **second-phase mechanistic result**, not a substitute
for locating a viable path. Launch a small 2D pilot only after a fast method
passes the chemistry/DFT checks and the relevant proton and bond coordinates
are identified. Check neighboring-window overlap, equilibration, forward/
reverse hysteresis, block stability and replica variation; extend windows
based on these diagnostics rather than a fixed 30 ps rule. If these gates pass,
compute an acylation PMF, then consider deacylation. A PES TS and a PMF
maximum answer different questions and should be reported separately.

## Manuscript scope and structure

The mechanism paper should include the **released PA6 parameter set and its
validation as a concise Methods foundation**, with versioned citation to
`/data/jpc/polyamide-amber-ff`; detailed modular residue derivation, broad
polyamide coverage and parameter benchmarks belong in that separate force-
field paper/SI. Include WT MD for E–S pose/protonation selection and sampling
limits. Include mutant MD only if it answers a predeclared mechanistic test;
otherwise leave redesign validation for a subsequent paper. The central
Results should be competing QM/MM paths, stationary-point validation, the
role of protonation/water, and DFT sensitivity. Include umbrella sampling in
the main paper if the convergence gates pass and it changes or strengthens the
mechanistic conclusion; otherwise report a bounded pilot as a limitation, not
a headline free-energy barrier. A full chemical-cycle claim requires both
acylation and deacylation, not only TS1.

Maintain `manuscript/INTRODUCTION_DRAFT.md` for prior NylC/Nyl12 work,
experimental motivation and competing hypotheses; `METHODS_DRAFT.md` for
executed protocols (main-text concise description plus SI detail); and
`RESULTS_DRAFT.md` for claims linked to exact run manifests and negative
results. Keep references in `REFERENCES.md`, claims/artifacts in
`SOURCE_LEDGER.md`, and decisions here. Move only mature, accepted findings
from the run ledger into manuscript prose. Do not paste long run logs into the
manuscript.

## Student MD comparison arm — 2026-11-19

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

## Operating rules

### Current bounded water-relay diagnostic (2026-11-17)

The accepted ASH274/ASP272 E–S endpoint with WAT1576 in the 102-atom QM
region is the sole starting structure. A 25-window DFTB3/MM scan drives the
composite Thr929 OG1–HG1 / water / Asp272 OD1 proton relay and OG1 attack on
the substrate carbonyl. The remaining MM atoms in the established 12 Å shell
are movable; no global positional restraints are added. Each accepted window
updates a stripped VMD trajectory and energy/distance plots. The runner stops
on failed convergence, proton identity, or geometry continuity. Only an
accepted endpoint triggers an unbiased minimization. This route would yield
Asp272-COOH and neutral Thr N: call it an **alternative TI1-like state**, not
canonical TI1, whose required Thr beta N is NH3+. Neither a biased scan nor
an unconverged endpoint establishes a stable intermediate or transition state.

The first execution stopped at window 5 on a false continuity rejection of
a distant ASH274–water separation; its reactive geometry and proton inventory
passed. Revision 5 preserves that original decision, checks only the seven
declared reactive distances for continuity, and resumes in place from window
6 after validating all completed output hashes. The path remains a diagnostic
until the endpoint and unbiased release pass their separate gates.


One active reference system, at most two competing routes, and one next
scientific decision at a time. Each queued run has a one-page hypothesis,
atom mapping, charge/proton inventory, method, restraints, input hashes,
acceptance gates and stop condition. Script and run in tmux; check only
meaningful milestones. Archive failed attempts without renaming them as
successful. Use a short weekly dashboard: completed, running, blocked,
next decision, and the single figure or structure the PI should inspect.
Change this roadmap only when a dated decision record cites new data and the
specific gate it changes. Avoid new calculations merely because hardware is
idle.

Selected methodological context: Negoro et al., *FEBS J.* 2023,
DOI 10.1111/febs.16755; Kästner and Sherwood, *J. Chem. Phys.* 2008,
DOI 10.1063/1.2815812. The manuscript reference ledger controls final
bibliography and claim verification.


### Native-Thr gate diagnosis — 2026-11-24

The rare joint native-Thr/Tier A/B passes are dominated by simultaneous
OG1/Nalpha displacement and psi deviations, rather than chi1 rotamer changes.
Among Tier A/B frames, WT versus mask_07 failure rates are 99.89/98.67% for
OG1, 99.84/99.03% for Nalpha, 83.82/90.21% for psi and 5.47/3.71% for chi1.
Modest diagnostic relaxation (1.0 Å OG1, 1.5 Å N, 45-degree psi, 30-degree chi1)
still admits only 2/1829 and 14/2779 frames. The original static WT starts are
borderline OG1 failures (0.757–0.880 Å versus 0.75 Å); the MD distributions
show larger real shifts. Independent alpha- and beta-scaffold fits do not
remove the discrepancy in a bounded WT sample. The original torsion
definitions/circular differences were verified; N-terminal phi is not used.

No thresholds changed. Native resemblance should be retained as a continuous
structural diagnostic, not treated as a validated activity classifier. The
crystal and ligand-bound MD need not have identical ensembles, and known-active
control rejection already limits hard-filter interpretation. Matched-site
representative and rare closest-native overlays, initial controls, scalar
sensitivity and portable PyMOL scenes are in
`design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/md_followup/analysis/native_thr_md_gate_audit_20261124/`.
