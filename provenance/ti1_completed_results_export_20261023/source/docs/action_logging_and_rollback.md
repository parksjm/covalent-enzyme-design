# Agent action journal and rollback

The user requested persistent logging of every agent action on 2026-10-21.
This applies to the primary agent and delegated agents, including reads,
commands, edits, job launches, permission escalations, tests, Git operations,
failures and retries. Record an intent before the action and an outcome after
it; journal writes themselves are represented by their own record rather than
recursively logging the logger. Record nonsecret command arguments in the
action/result text or retain an exact command manifest and reference it.

Use `scripts/record_agent_action.py`. Records include UTC time, agent, phase,
reason, permission mode, rollback guidance and optional before/after file
hashes. New schema-version-1 records have unique UUIDs and strict validation.
Concurrent appends use a file lock. Do not log credentials or entire environments.
If the journal cannot be written, stop dependent mutations, report the failure,
and fix the journal path first. If a tool cannot execute at all, its tool error
remains in the conversation transcript; append an honest failure record once
the logging path is available. Never invent a successful action or rewrite
history to imply that an action was logged earlier than it was.

The session bootstrap is preserved in
`provenance/agent_actions_20261021.jsonl`. Early records have heterogeneous
legacy fields; validation tolerates those explicitly and applies the strict
schema to new records. A documented bootstrap gap occurred when a logger
validation error failed to stop a shell that launched the second read-only QM
audit. Its original tool transcript, launch script and job provenance remain;
the gap is disclosed by a subsequent journal note. Shell batches now use
`set -e` so failed journal writes stop dependent commands.

At the Git checkpoint, freeze the versioned journal as a snapshot and copy it
to **`output_agent_actions/20261021.jsonl`**, which becomes the live journal for
subsequent actions, including commit/push outcomes. This avoids an endless
cycle of commits that record their own commits. A later checkpoint can version
a new snapshot; never overwrite the earlier snapshot. The live journal is
retained on violet and excluded from ordinary code commits by `output*/`.
The journal is a workflow audit trail, **not an operating-system security
boundary** or a guarantee against modification by an account with write access.
Permission controls and the retained tool transcript remain separate evidence.

Example (use the active journal path):

```bash
python3 scripts/record_agent_action.py \
  --log output_agent_actions/20261021.jsonl --agent root --phase intent \
  --action 'Run uv run pytest' --permission-mode sandbox \
  --reason 'Verify changes before commit' --rollback 'Tests use temporary fixtures'
```

For rollback:

- Code/document changes: inspect the commit diff and use a new `git revert`
  commit when authorized. Preserve unrelated work; do not use a destructive
  reset as a substitute for review.
- Scientific jobs: retain immutable inputs, source snapshots, exact command
  manifests, logs and acceptance records. A failed attempt remains evidence;
  a new attempt uses a new directory. Stop only a specifically identified job
  when requested, preserving its terminal/process records.
- Cache renames: consult `output_parameter_cache/.migrations/` and
  `provenance/cache_directory_migration_20261021.json`. Each entry records its
  original and descriptive name. Verify manifest/payload hashes, confirm the
  original path is free, then explicitly reverse that one rename if needed and
  append a reversal record. Do not change the cached manifest or delete the
  migration history.
- Laptop exports: original structures remain in their original run directories.
  The download manifest maps every copied file to its source and digest. A
  correction produces a new export; local rsync never uses `--delete`.

Hashes detect changes; recovery requires the retained artifacts or Git history.
