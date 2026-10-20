# Current covalent-pipeline adapter

As of 2026-10-20 the development checkout is
`/data/jpc/covalent-enzyme-design` (private `parksjm/covalent-enzyme-design`).
Read `docs/standalone_handoff_20261020.md` first. Existing scientific artifacts
remain under `/data/jpc/enzyme-substrate-modeling`; historical relative output
paths in the original resume note refer to that older workspace. Pass explicit
paths while workstation-specific defaults are being removed.

Entry point: `scripts/run_nyl12_job.py`.
Implementation: `src/nyl12_covalent/jobs.py`.

```bash
python scripts/run_nyl12_job.py launch \
  --session unique_calculation_name \
  --output-dir /absolute/new/job-directory \
  --cwd /absolute/calculation-directory \
  --input /absolute/input-file \
  --env OMP_NUM_THREADS=4 \
  -- /absolute/environment/bin/python /absolute/calculation-script.py
```

Replace every example value with the inspected task's values. Select the
launcher Python environment explicitly as well. Multiple `--input` and `--env`
options are supported. The launcher creates `run.sh`, `command.json`,
`status.json`, and `job.log`. Input hashes cover only explicitly supplied files;
include relevant scripts and imported modules or an immutable source manifest.
Environment overrides are recorded verbatim, so do not pass secrets.

```bash
python scripts/run_nyl12_job.py status /absolute/job-directory
tmux list-sessions
```

Inspect calculation-specific files alongside the launcher log. Scientific
failures may be reported in an engine's own output rather than `job.log`.

The current `qmmm.py` functions are initial-run helpers, not checkpoint-resume
interfaces: they read the original `*.inpcrd` and reuse output names. Do not
rerun one to continue an interrupted optimization. First establish that a
restart matches its topology and scientific specification, then construct an
explicit continuation command using that restart as `-c` in a new attempt
directory. Record the original positional-restraint reference separately;
changing it changes the optimization problem.

The current Nyl12 handoff is
`docs/nyl12_covalent_resume_20261019.md`.
On migration, update this adapter to the new entry point without copying the
launcher implementation into the skill.
