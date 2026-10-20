---
name: scientific-tmux-jobs
description: Launch, inspect, and resume reproducible computational-science jobs in tmux with input hashes, explicit environments, and retained logs. Use for local or workstation calculations where the user requires persistent scripted execution.
---

# Scientific jobs in tmux

Use the project's existing job launcher when available. Read
[the current pipeline adapter](references/pipeline.md) for its location and
interface. Keep scientific setup and acceptance checks in versioned project
code; this skill describes how to execute and preserve them.

## Before launch

- Read the current resume note and inspect existing job status and tmux sessions.
  An absent session does not establish success or justify rerunning a calculation.
- Identify the exact executable, working directory, input files, output directory,
  and environment overrides. Preserve a virtual environment's executable path:
  resolving its Python symlink can silently select a different environment.
- Check software access and required scratch directories on the execution host.
  Record missing prerequisites and continue independent work. Do not substitute
  another Hamiltonian, charge method, or scratch location without an explicit
  scientific decision.
- For an expensive parameter calculation, check the validated content-addressed
  cache first. Match chemistry, atom mapping, charge, protonation, method,
  conformer, and tool versions. A corrupt entry is an error.

## Execute and record

Launch the exact argument vector in a uniquely named detached tmux session.
Use a new output directory for each attempt. Preserve:

- The runnable script and command arguments.
- Input and relevant source-code SHA-256 hashes.
- Explicit, nonsecret environment settings and software versions.
- Start/end timestamps, exit status, and combined stdout/stderr.

Avoid recording credentials or dumping the whole environment. Do not alter
unrelated sessions. Launching a local command and submitting work to a remote
host remain subject to the user's task scope.

Hashes detect changes but cannot recover overwritten files. Preserve an
immutable source snapshot or committed revision plus patch for uncommitted
calculation code, and retain the actual input artifacts referenced by hashes.

## Inspect and resume

Check the status file, calculation log, and expected outputs. A zero exit status
or a restart file establishes neither convergence nor valid chemistry. Apply
the calculation-specific convergence, geometry, topology, and charge checks;
record their results separately from process completion.

For interrupted jobs, retain the failed attempt and record the exact checkpoint
and source attempt used for continuation. If a process vanishes without a final
status, mark it unresolved and inspect logs before deciding how to continue.
Do not overwrite a previous attempt or retry a deterministic failure unchanged.

Update the project resume note with artifact paths, completed checks, scientific
limitations, active sessions, and the next action. Keep changing run state in
project provenance rather than in this skill.
