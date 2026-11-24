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
- With consistent scissile-amide mapping, opposite catalytic Herschlag-phi
  signs identify flipped local poses. Other distances assess catalytic utility.
- Use independent MD seeds and matched comparison protocols. Update topology
  charges and counterions for protonation changes; never change states silently.
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
