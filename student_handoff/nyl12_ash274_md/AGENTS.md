# Student-agent scope: Nyl12 ASH274 MD screening

Read START_HERE.md, MEMORY.md and STUDY_PLAN.md before acting. This directory
scopes the student MD arm; historic parent QM/MM plans are context, not orders
to launch additional calculations. Follow the current human supervisor's
instructions. Do not alter running PI jobs or contact collaborators automatically.

- Establish matched WT controls before making a mutation/design claim.
- ASH274 occupancy superiority is a hypothesis, not a settled finding.
- Keep 300 K, fixed microstates, exact released PA6 L4 RESP/GAFF2, ff19SB/TIP3P,
  neutral mature Thr termini and the validated MD protocol unless explicitly
  changing one declared variable. Update topology/ions/charge for every state.
- Require accepted, hash-matched prebuilt topology/restart pairs. Never infer
  chemistry from a filename, silently substitute parameters or fall back to CPU.
- Use the explicit precursor/chain-local/Amber residue maps; never transfer raw
  atom indices across mutant topologies. Protect all symmetry copies of the
  catalytic set and maturation junction. Fixed length in the first design arm.
- Reproduce site-AE historical metrics, then implement/report four-site-aware
  occupancy separately. Keep phi signed; document the exact scissile atom map.
- Analyze independent replicas and all attempted runs, not only successful
  frames. Binding-pose occupancy is neither a pKa nor a catalytic rate.
- Preserve all source artifacts. Run authorized calculations in tmux with scripts,
  explicit software/env/commands, input/source SHA values and durable logs.
  Keep process completion, preparation acceptance and scientific acceptance apart.
- Log intent and outcome for every agent/tool action in the student's append-only
  action journal. Retain rollback guidance and failed attempts. Do not overwrite
  delivered manifests or PI records; create a declared new dataset/version.
- Use uv for Python dependency management; ask before adding new production
  dependencies. Existing scientific binaries need their own official installation.
- Use typed small modules, logging, Google-style public docstrings, Ruff and
  meaningful pytest checks. Delegate bounded routine work to approved cheaper
  models. Chinese LLMs are excluded by the PI's project constraint.
- No production LigandMPNN campaign until WT/microstate and analysis gates are
  reviewed. E–S exploratory design is distinct from validated TS-guided redesign.
- No new QM/MM pathways, force-field refits, constant-pH MD or public data release
  is implied by this handoff. Such work needs a separate explicit scientific scope.
