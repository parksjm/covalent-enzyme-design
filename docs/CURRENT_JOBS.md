## Focused update — relay reactant coverage, 28 November 2026

Six additional DFTB3/MM umbrella histories started on Geranium at 15:56 UTC
(10:56 EST), in tmux `nyl12_es_boundary6_20261128`. Six native Sander PIDs
were verified against their exact working directories. They use three new
reactant-side centers and two E–S starting geometries per center, with 30 ps
production after 1 ps fresh-velocity equilibration at 300 K + 4 ps continued
equilibration. This is running sampling,
not completed or scientifically accepted sampling. Estimated wall time is
about 20 hours from measured prior throughput; no GPU is allocated.

All existing Violet/Azalea relay US jobs are preserved. The new branch is
`qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/activation_support_boundary_20261128/`.
Its hashed launch receipt is under
`analysis/reactant_coverage_plan_20261128/preparation/launch_evidence_20261128.json`
in the same campaign. Separate automatic checkpoint analysis is running in Violet tmux
`relay_boundary106_mbar_20261128`, with common new-production checkpoints at
10/15/20/25/30 ps. It compares 106 histories/95 biases with a frozen original
70 ps baseline and preserves the primary collector and current links. The
first new checkpoint retains one 5 ps chunk per new history; it is a point
and mixing diagnostic, without a qualified pooled confidence interval.
See `NYL12_MECHANISM_ROADMAP.md` for current scientific priorities. Entries
below are historical snapshots, not a fresh full-host status audit.

## CURRENT — resumed pKa controls, 2026-11-23

Three NylC MCCE jobs (5XYG,7YU0,7YU1) have passed preparation and conformer
generation and are running DelPhi step3 on Violet in their per-case
`rerun_20261123/` directories. Exact sessions/current pointers and acceptance
gates are in `NYL12_MECHANISM_ROADMAP.md`.

The five-call TS/native endpoint audit finished successfully: recovered SCC
inner trials were misclassified by the old parser. Fresh endpoint gradients
pass, and the two minima classify as E-S/acylenzyme. The TS candidate meets
its original force tolerance but remains unverified because of the Hessian
soft-mode uncertainty. Twelve targeted column evaluations are recommended,
not queued. The original extended-CG ASH274 NEB remains untouched. Azalea's
Mistral server remains stopped; no new mechanism-agent job is queued.

## CURRENT — live audit, 2026-11-23

Machine-readable audit: `provenance/current_jobs_20261123.json`. This section
supersedes every older running/queued snapshot below.

- **Running science:** one Violet Sander engine, extended-CG ASP272/ASH274
  NEB. Latest coherent audit: cycle 6336, relative PES peak 44.838 kcal/mol;
  not converged or a verified activation barrier. No science is active on Azalea.
- **Finished:** all 12 matched WT MD replicas (4/state × 100 ns) and queued
  occupancy/relay analyses; original45-model and new55-model structural pilots;
  63-window NEB-seeded US (30 ps/window), final MBAR and movie exports.
- **MD result:** mean Tier A/B occupancies are NH2/ASP/ASP 27.77%,
  NH2/ASP/ASH 5.18%, NH3/ASP/ASP 28.99%. Native-Thr retention is rare in all
  states; fixed-state occupancy does not assign protonation or validate a mechanism.
- **Still unresolved:** first-order TS verification (soft-mode/electronic
  diagnostics failed), qualified assisted proton relay, and NylC MCCE preparation
  (`REMARK` parser failures remain in current attempts). No new jobs queued.
- **Model servers:** Mistral/vLLM on Azalea was gracefully stopped at the PI's
  request after completed benchmarks; both GPU process lists are empty, with
  about95 GiB free each. Ollama remains running without GPU allocation. Weights,
  benchmark responses and deployment scripts are preserved.
- **Delivery/repo:** VMD syntax/frame fixes and 433MiB movie transfer fix are
  pushed in `51dbdb1`. The last implementation checkpoint passed130focused tests;
  this audit changed no scientific inputs and required no tests.

## Matched-mask campaign complete — 2026-11-22

All 52 new Boltz predictions and Amber minimizations completed; with three
qualified reused WT seeds, all 55 endpoints reached the force tolerance.
Reporting-only recovery completed with exit 0 in tmux
`nyl12_mask_publish_20261122`; see
`design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation/COMPLETED.json`
and `validation/reporting_repair/job_final/`. Original wrapper/report failures
are preserved. No calculations were repeated and no follow-up MD is queued.

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

# Current running and queued work

## Campaign-specific receipt snapshot — 2026-11-22

The Nyl12 structural-validation campaign receipt reports `running`, with 21/45
models normalized and 15/45 relaxed. This is campaign progress only; this update
did not audit live processes across the workstations. WT seed 01 minimization
converged, and its one-seed four-site geometry preflight is recorded in
`design/nyl12_ligandmpnn_pilot/boltz_validation/verification/checkpoint_receipt.json`.

PI authorized Tulip prediction and Geranium Amber relaxation for this campaign
after their availability check and confirmed both have the same architecture as
Violet. Existing runtimes sufficed; no dependency installation was needed. Both
Boltz prediction hosts use the same checkpoint SHA256, with runtime-build details
in the campaign's `dispatch/runtime_versions.json`. This authorization is scoped
to the structural campaign; no current host-wide liveness is asserted here.

## Latest NEB-US MBAR refresh — completed 2026-11-22

The requested MBAR refresh completed with exit 0 using 63/63 windows and six
accepted 5 ps blocks per window (30 ps/window; 1,890 ps aggregate; 37,800 raw
and 7,871 thinned frames). Sampled overlap is connected; target ESS is 7.2928.
All 50 connected 5 ps block bootstraps completed. The apparent maximum is
47.0475 kcal/mol relative to the same fixed E-S density bin as the prior
46.9414 kcal/mol result (+0.1062 kcal/mol). The fixed-peak-bin bootstrap 95%
interval is 45.4070–47.8696 kcal/mol and is not comprehensive PMF uncertainty.
Windows 028, 030 and 032 are flagged for drift.

This PMF remains unconverged and is not a validated activation barrier. The
refresh changed no sampling jobs and authorized no further extension. Output is
under `qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/umbrella_sampling/analysis/pmf_refresh_20261122/`.
The stable `analysis/pmf_latest/` pointer now targets its `compact_export/`;
the previous export is preserved. Frozen source/input hashes, same-bin numeric
regression and the independent 2D marginal passed publication checks. Updated
plot: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/umbrella_sampling/analysis/pmf_latest/primary_q_acyl_pmf.png`.

## Closeout — historical neutral Thr/Asp272/Asp274 direct-route NEB-seeded US — 2026-11-22

Campaign root: `qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/irc_0K/new_neb/trace_seed_25img_step005_max4000_20261113/`.

All 63 windows completed six accepted 5 ps blocks (30 ps/window). Violet's
eight workers exited 0 on November 22 between 10:43 and 10:47 UTC; azalea's 24
workers exited 0 on November 21 between 11:55 and 20:49 UTC. The latest
coverage-first receipts are under
`umbrella_sampling/dispatch/coverage_first_20261120/{violet,azalea}/`.
Only the local analysis PID 4100242 and sync PID 4098079 were stopped with
SIGTERM; both were absent afterward and the receipt is
`umbrella_sampling/dispatch/stop_20261122/violet_monitor_stop_receipt.json`.
No campaign Sander remained active on either host. Sampling completed; sampling
convergence did not.

The refreshed apparent PMF maximum is 47.0475 kcal/mol versus 46.9414 in the
same E-S density bin, with connected overlap but target ESS 7.29 and drift in
windows 028/030/032. It is not a converged PMF or validated barrier. Decision:
stop additional investment in this forced direct route. Preserve raw data,
restarts, sources and this negative finding; no extension is queued. This
closeout applies only to the historical neutral-Thr, ASP272/ASP274 direct-NEB-
seeded US at the long `trace_seed_25img_step005_max4000_20261113/` campaign
root. Other jobs, design work and ASH274 NEBs are unaffected.

The verified movie is available at
`umbrella_sampling/vmd_dynamic/view_umbrella.vmd`: 30 frames/window,
1,890 frames and 19,967 stripped solute atoms. All selected frames reproduce
the original four reaction distances within 0.00005 Å. Bundle checksums pass;
the laptop pull includes the final bundle and excludes duplicate host parts.

## Historical host-work snapshot — 2026-11-21 09:22 PM EST

The following inventory is from the last full live-process audit on violet and
azalea. Treat it as dated history, not current liveness. Machine-readable
historical receipt: `provenance/current_jobs.json`.

| Work | Status at 2026-11-21 audit | Queued work at that time |
| --- | --- | --- |
| Matched WT MD, violet GPUs | 10/12 accepted 100 ns replicas. NH2/ASP/ASP r04 at 33.7 ns; NH3/ASP/ASP r03 at 40.5 ns. | These are the final two replicas; native remaining-time estimates 5.9/5.2 h at inspection, approximately 03:15/02:30 EST November 22. Estimates may change. |
| Post-MD analysis, violet | Live controller waiting for all 12 accepted replicas. | Tier A/B, native Thr and signed-phi comparisons across four replicas/state, plus new-four relay screening and matched comparison. |
| NEB-seeded US | Eight violet Sander workers running 5 ps production block 05. No native jobs remain active on azalea. | Eight other violet windows queued (025,039,050,006,015,021,028,047), then unfinished blocks through 30 ps/window. |
| ASP272/ASH274 NEB | Five live Sander engines on violet: original plus four spacing/optimizer trials. | Continue their existing budgets; no extra paths queued. None has current scientific convergence acceptance. |
| Background support | US synchronization/MBAR watcher, NEB plot/movie watchers, MD queue/completion controller. | Automatic updates use completed accepted data; monitors are not additional independent scientific campaigns. |

US accepted coverage as of 18:37 EST: **1,770 ps**, all 63 sampled, with 8 windows
at 20 ps, 8 at 25 ps and 47 at 30 ps. Last separately frozen/plotted primary PMF
used 1,575 ps and had apparent maximum **46.94 kcal/mol**. It remains preliminary;
newer accepted sampling is not silently substituted into that frozen figure.

Latest spacing-trial records (PES maxima, not verified barriers):

| Trial | Cycle | Relative peak, kcal/mol | Maximum gradient / target |
| --- | ---: | ---: | ---: |
| Uniform 24, L-BFGS | 2069 | 33.922 | 16.0 |
| Uniform 48, L-BFGS | 1093 | 51.900 | 39.3 |
| Adaptive 24, L-BFGS | 2076 | 45.645 | 78.5 |
| Uniform 24, extended CG | 3507 | 42.687 | 2.90 |

Trial chemical-gap acceptance is unavailable because coherent-export cycles
mismatch the numerical records. The original band's exported profile is stale;
its native engine is live. Do not claim convergence from these records.

Completed or stopped:

- TS directional follow-up completed exit 0: 26/26 native probes, two
  optimizer-converged diagnostic descents. Both descent electronic gates failed;
  **no verified first-order TS or formal IRC**. No further TS calculation queued.
- NylC MCCE comparison controllers for 5XYG, 7YU0 and 7YU1 all stopped with
  `ValueError: Unsupported PDB record 'REMARK' on line 1`. No active rerun or
  accepted new pKas; preparation/parser repair is pending.
- Mechanism-agent model benchmarks and semantic review completed. Neither
  model is qualified for autonomous scientific dispatch. Ollama/vLLM servers
  remain up on azalea; no scientific task is being dispatched by them.
- Email notifications still await external ORNL relay authorization.
- No Gaussian, xTB or new parameterization calculation was found running on
  violet or azalea. No jobs have been assigned by us to geranium or tulip.

Original inputs, outputs and failed attempts are preserved. Resource-use
restrictions remain in AGENTS.md; this audit did not launch or migrate work.
