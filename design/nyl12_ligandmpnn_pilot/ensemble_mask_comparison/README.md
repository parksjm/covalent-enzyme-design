# Matched focused versus broader sequence-design comparison

This PI-authorized test compares the original 19-position mask with a nested
29-position mask on two distinct, force-converged WT backbones. Both parents
have at least one occupied site passing Tier A/B and the native Thr233 gate.
The four active sites within a tetramer are nested observations. These two
minimized bound conformations sample limited E–S heterogeneity; they are not
an open/closed ensemble and do not establish substrate-access flexibility.

The sequence experiment has 2 masks × 2 backbones × 2 temperatures × 6 samples
= 48 raw designs (24 per arm), with four RNG seeds paired across mask arms.
All sequences, duplicates and lineage are retained. Alpha A–D and beta E–H
positions are explicitly tied; no indels are allowed. Catalytic/maturation
protections are unchanged. Ten new alpha positions are eligible: 28, 37, 38,
92, 96, 97, 160, 161, 162 and 164. Beta positions remain 39 and 41.

New positions require relative solvent accessibility ≥0.25 in every equivalent
copy of both parents and no other-chain heavy-atom contact ≤4 Å. The existing
focused mask remains unchanged. These are conservative design-mask heuristics,
not activity-trained filters or a claim that these sites control flexibility.
RDKit's bundled FreeSASA (Lee–Richards/Protor, 1.4 Å probe) computes exposure;
normalization follows [Tien et al., 2013](https://doi.org/10.1371/journal.pone.0080635).
No dependency installation was needed. Full atom/position decisions and source
hashes are in `preflight/`; the accepted generation input is
`prepared/generation_manifest_complete_runtime.json`.

The follow-up selects at most four distinct sequences per arm, balancing parent
and temperature strata through a deterministic sequence-diversity rule, plus
matched K59D and V111C positive controls. These controls have noisy lysate/MS
product evidence, not measured intrinsic turnover enhancements. K59E remains
unresolved and is not used as a negative control. Five matched Boltz seeds,
four PA6 L4 molecules, crystal-template/contact conditioning and the existing
RESP/GAFF2/ff19SB/TIP3P relaxation protocol are retained. Reuse of completed WT
controls requires scientific-input and runtime identity; the two earlier Tulip
WT predictions are regenerated on Violet to remove that software-build factor.
Maximum new prediction budget is 52 models; duplicates can reduce it.

Only Violet is assigned for new CPU generation and GPU validation. Existing
scientific jobs remain intact. New artifacts stay in this comparison directory.
No mutation-dependent WT total charge or counterion count may be assumed:
sequence-derived charge validation is required before new Amber preparation.

This is an E–S structural-diversity and filter-calibration pilot. More passing
models may reflect better modeled preorganization but do not establish lower
chemical barriers, transition-state stabilization or improved polymer turnover.
There is no verified TS conditioning, no fitted activity cutoff and no claim
that sequence symmetry guarantees stable tetramer interfaces.

Operational state belongs in job receipts/status files. The initial named-flag
CLI startup failure is preserved in `generation/job/`; it occurred before engine
execution. The first native attempt omitted the already installed frozen OpenFold support
path and failed before sampling; its receipts are also retained. The complete
runtime uses that same verified support with no installation; its persistent
job is `generation/job_complete_runtime/`.

Sequence generation completed: all eight native jobs exited successfully and
all 48 sequences passed the protection/symmetry audit (39 distinct sequences
across arms). The broader arm had 24/24 unique sequences; focused had 15/24.
A CSV publication error occurred after successful sampling; the original
wrapper exit 1 and source are retained. Corrected export and an independent
re-audit succeeded without resampling; `generation/results/acceptance.json`
is the explicit scientific handoff gate. This is separate from native engine
completion and from downstream structural acceptance.

Absolute MM potential energies are not used to rank different sequences or
net charges. The matched endpoints are geometry, assembly integrity, confidence
and structural diversity; chemical barrier changes require separate reactive
state validation and experimental testing.

The full follow-up ran in tmux
`nyl12_matched_masks_validation_20261122`. All 52 new predictions and
minimization processes have now finished; all 55 endpoints, including three
reused WT models, reached the force tolerance. The final plotting phase
rejected the new candidate names; reporting-only recovery is retained under
`validation/reporting_repair/`. No scientific calculations were repeated. The authoritative completion receipt is
`validation/COMPLETED.json`; the reporting-only job is
`validation/reporting_repair/job_final/`. The original outer job failure is
retained under `workflow_job_checked/`. It contains
52 new predictions and three eligible reused WT models (55 model seeds total).
Both Violet GPUs were assigned to prediction first, then to matched Amber
relaxation. The four-site comparison is now complete. The first
preflight hash-dictionary mismatch stopped before GPU work; its status and
source snapshot are preserved, and the corrected source is frozen separately
in `validation_source_snapshot_checked/`. All 52 job records and all five
hashes of each reused WT model passed the native preflight before resumption.

## Completed geometric comparison (2026-11-22)

Each candidate has five model seeds with four occupied sites per model. Joint
pass means Tier A/B geometry plus the native-Thr233 conformation gate. Model
seeds and sites are nested within four selected sequences per design arm;
pooled site percentages are descriptive and do not establish significance.

| Cohort | Joint sites | Models with a joint site | Tier A/B sites | Native-Thr sites |
| --- | --- | --- | --- | --- |
| WT | 1/20 | 1/5 | 10/20 | 5/20 |
| Focused 19-position mask | 8/80 | 4/20 | 24/80 | 21/80 |
| Broader 29-position mask | 12/80 | 7/20 | 20/80 | 43/80 |
| Known-active K59D control | 0/20 | 0/5 | 10/20 | 6/20 |
| Known-active V111C control | 7/20 | 3/5 | 10/20 | 16/20 |

Consult `validation/analysis/comparison_summary.json` for independently audited
counts and source hashes; all 55 minimizations are force-converged.
The broader arm has more joint passes, primarily alongside better native-Thr
retention, while its Tier A/B-only frequency is lower than the focused arm.
It is therefore not evidence for uniformly improved active-site geometry.

Signed Herschlag phi remains essential: restricting the joint passes to the
negative-phi orientation leaves broader 11/80, focused 3/80, WT 0/20 and V111C
7/20. Positive-phi joint passes must be reported separately instead of being
silently treated as the same mechanistic pose. No gate was retuned to obtain
these numbers. The WT baseline differs from the earlier mixed-runtime
45-model campaign because two WT seeds were regenerated on Violet.

Broader `mask_07` is the most repeatable designed candidate: 6/20 joint sites
in 3/5 seeds, all at negative phi. Focused `mask_03` gives 3/20 in 2/5 seeds,
also at negative phi. Focused `mask_05` has 4/20 jointly passing sites, all in
one seed and at positive phi. Its apparent advantage is therefore conditional
on seed and substrate orientation. K59D's zero joint passes despite activity
in noisy lysate/MS assays indicates that the joint gate is not validated as an
activity classifier. Preserve it as structural characterization rather than
an absolute rejection rule for enzyme activity.

Next decision gate: inspect representative negative-phi mask_07/mask_03,
WT and V111C structures, confidence and assembly contacts before any matched
MD follow-up. No new simulation campaign or mechanistic claim follows
automatically from these static geometric counts.
