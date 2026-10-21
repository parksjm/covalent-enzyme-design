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
