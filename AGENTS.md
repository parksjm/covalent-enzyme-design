# Working agreements

- Use `uv`; ask before adding new production dependencies.
- Clarify the chemistry and intended products before changing implementation.
- Use typed function signatures, Google-style public API docstrings, and Python
  `logging` for operational messages. Follow PEP 8 and use Ruff.
- Keep modules small, preferably under 300 lines. Preserve entry points and
  isolate enzyme-specific configuration.
- Fail explicitly on missing prerequisites or invalid scientific inputs.
  Do not substitute methods, protonation states, or parameters silently.
- Run authorized scientific jobs in tmux with scripts, input/source hashes,
  environments, logs, and separate process/scientific acceptance records.
- Preserve all original artifacts and unrelated user changes. Delete backups
  only after verification succeeds, using safe deletion.
- Use `pytest`; run `uv run pytest` and relevant formatting/lint checks after
  edits when feasible. Ensure at least 80% coverage of new code and report
  when coverage cannot be measured.
- Keep responses compact: changes, checks run, and a short rationale.
- Record every agent/tool action, including delegated work and failures, in the
  active append-only action journal. Record intent before execution and outcome
  afterward; include permission mode, reason, affected-file hashes and rollback
  guidance. If logging fails, stop dependent mutations and repair the logging
  path. See `docs/action_logging_and_rollback.md` for the active journal and
  checkpoint procedure. Keep the user informed of meaningful changes.
- Use descriptive artifact/cache directory names. Keep identity hashes in
  metadata; do not use hashes as the normal user-facing directory names.

Read the current resume note before calculations. TI1 protonation/connectivity,
measured product identities, and parameter quality are scientific requirements.

For this Nyl12 campaign, start with `docs/NYL12_THROUGH_NOV30.md` for active
priorities and compact status/evidence links; consult the long historical
roadmap/handoff only for the relevant unresolved question. From 1 December,
prefer deterministic checkpoint summaries and milestone/decision reviews to
routine repeated LLM monitoring, preserving scientific checks and provenance.
The initial repository import is a development snapshot; it does not certify
the existing topologies for MD or the current coordinates as optimized.

## Research efficiency

- Delegate bounded scripting, extraction, plotting and monitoring to lower-cost
  models. Give explicit inputs, outputs and checks; avoid full-history forks.
  Keep chemistry decisions and final interpretation with the lead agent.
- Reuse validated tools and parameter releases. Read relevant checkpoint
  sections, parse large files locally, and return concise tables or excerpts.
  Preserve full evidence on disk rather than repeatedly loading it into chat.
- Keep JSON/JSONL authoritative and CSV for tables. Reduce prompt payloads
  before considering a new serialization format or dependency.
- Prefer focused tests for changed behavior; run the full suite for shared
  infrastructure changes or release checkpoints. Do not repeat passing checks
  without new evidence. Inspect documentation/plot edits without adding tests
  that merely mirror them. Report any unmeasured code coverage.
- Define the hypothesis, protonation, atom mapping, method and acceptance gates
  before calculations. Screen existing models before launching new campaigns.
  Keep comparisons bounded and change one variable at a time where practical.
- Distinguish job completion, force convergence, chemical retention, geometric
  acceptance and sampling convergence. MD preparation need not meet the force
  tolerance required for a stationary-point calculation.
- Distinguish PES profiles from PMFs; an unconverged NEB maximum is not a
  validated barrier. Record all restraints, frozen atoms and proton inventories.
  Separate substrate hydrolysis from enzyme maturation/autocleavage.
- QM/MM PMFs and verified PES transition states are not prerequisites for PI
  experimental nomination based on sequence, assembly, structural and dynamic
  evidence, with expression/maturation/activity checks. Chemistry may precede
  testing if resources permit or follow improvement to explain mechanism.
  Computed barrier-improvement claims require matched, supported WT/design
  calculations; do not tune post hoc protocols to reproduce assay outcomes.
- RFD3 remodeling must couple mature alpha/beta repeats A/E, B/F, C/G and D/H
  with symmetry-equivalent motion masks and the native symmetry sampler.
  Approximate crystal D2 NCS does not qualify exact operators/native export or
  fixed single-site TS compatibility. Independently check final backbone
  residuals, fixed context and ligands. Native sampling projects denoised
  coordinates; raw Gaussian noise is independent. Historical chain-A-only,
  default-sampler parents omitted symmetry: their lineages cannot advance
  pending independent symmetry qualification. Immutable jobs may finish as
  diagnostics; original MD-backbone baseline lineages remain eligible.
- With consistent scissile-amide mapping, opposite catalytic Herschlag-phi
  signs identify flipped local poses. Other distances assess catalytic utility.
- Before template-guided prediction, verify actual native feature tensors,
  source-author/query mapping and omitted-position masks. Coordinate deletion
  and alignment records alone do not prove effective template visibility.
  Test all matched conditions and preserve source coordinates; see the
  27 November native-feature correction in `docs/NYL12_THROUGH_NOV30.md`.
- For the current Nyl12 E–S design campaign, the PI-approved Boltz pose is
  **positive phi** under the unchanged projected-angle definition. Require an
  explicit `target_phi_sign=+1` and the approved reference manifest; the earlier
  assistant negative-only target was incorrect. Preserve numerical signs and
  historical analyses. Follow the PI's default protein-based active-site camera
  in `docs/PRESENTATION_AND_FIGURE_STYLE.md` for figures and movies.
- Use independent MD seeds and matched comparison protocols. Update topology
  charges and counterions for protonation changes; never change states silently.
- Preserve the hard Nyl12 sequence minimum at canonical 113, 156, 186, 221–235
  and 272–276 in all four copies, plus structural/interface protections. See
  `nyl12_sequence_protection.py` and the latest audited 101-position receipt.
  K156 is protected by the PI's latest review of NylC Lys189 mutagenesis.
  L235 remains a conservative N-terminal-context protection without a
  demonstrated residue-specific catalytic role. Earlier K156-variable pilot
  outputs are diagnostics, not eligible candidates. Historical
  config claims about275–276 are protection policy, not mechanistic proof.
  Sequence identity and RFD3 coordinate mobility are separate. Current designs
  violating restored protections may finish as diagnostics but cannot advance.
- For local MSA searches, read `docs/LOCAL_MSA_WORKFLOW.md` first. Violet's
  working Pixi ColabFold is under enzyme-substrate-modeling and databases are
  in `/data/boltz_db`; do not infer failure from a broken global wrapper or an
  empty unrelated database directory. Native output uses the query header.
- Use scripted monitors and meaningful checkpoints, not frequent manual polls.
  Check CPU/GPU contention before stopping work; signal only verified processes.
- Default science hosts are violet and azalea. On 2026-11-22 the PI explicitly
  reauthorized geranium and tulip for the bounded 45-model LigandMPNN-variant
  validation campaign: violet/tulip predict with Boltz-2; geranium relaxes with
  Amber. This supersedes the earlier two-host restriction for that campaign.
  Aster and other workstations still require renewed PI authorization.
- Before any Amber GPU launch, read `docs/AMBER_WORKSTATION_RUNTIME.md`.
  Azalea's Blackwell GPUs require the already validated CUDA-12.8 executable
  `/data/jpc/code/pmemd26_blackwell_build_20261005/pmemd26_src/build_sysgcc_cuda128_v2/src/pmemd/src/pmemd.cuda_SPFP`.
  The generic `/data/jpc/code/pmemd26/bin/pmemd.cuda` works on Violet's Ada GPUs
  but fails on Azalea with "no kernel image is available". Existence, `ldd`
  success and an idle GPU do not establish CUDA architecture compatibility.
  Recover the exact successful host-specific command, hash and environment;
  never rediscover this by launching another generic-binary campaign.
- Geranium and tulip are shared: use both simultaneously
  only for work expected to finish within a few hours. Before each launch,
  verify that nobody else is
  using the proposed host (sessions, CPU/GPU processes and queued work), not
  just that a GPU is idle. Ask the PI when availability or duration is unclear.
- Keep pathways in their existing logical directories, preserve immutable
  inputs, and record continuations clearly. Archive obsolete data in `clutter`.
- Preserve references, negative results and decision rationale for publication.
  Numerical derivative or installation-wide checks need a specific diagnostic
  reason; prioritize validation of our custom code and scientific inputs.
