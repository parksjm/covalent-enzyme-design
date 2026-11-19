# Nyl12–PA6 protonation and design-screening handoff

Date: 2026-11-19. This is a **student MD subproject**, complementary to the
ongoing Nyl12 mechanism/QM/MM project. It does not replace that project's queue.

## Read in this order

1. `MEMORY.md`: the short scientific checkpoint.
2. `STUDY_PLAN.md`: scope, controls, milestones and decision gates.
3. `STUDY_SUMMARY.md`: manuscript-style background, methods and findings.
4. `SOFTWARE_AND_REPRODUCIBILITY.md` and `tools/README.md`: environment/run setup.
5. `evidence/`: audited prior counts, mutant inventory, sequences and hashes.
6. `export_manifest.json` and `datasets/manifest.json` in the extracted archive:
   delivered file hashes and the three explicitly prepared WT microstates.
   The repository `data_manifest.json` records the original workstation sources.

`AGENTS.md` in this directory scopes work by the student's agents. The parent
repository's `docs/MASTER_PLAN_NYL12_COVALENT_REDESIGN.md`, mechanism roadmap,
standalone handoff and manuscript/source ledger preserve the wider history.
The archive includes that code/documentation snapshot under `project_code/`.
Earlier dated plans are context, not instructions to launch every old idea.
No previous ChatGPT conversation is required to understand this subproject.

## Deliverables and access

The private code/documentation repository is
https://github.com/parksjm/covalent-enzyme-design . The authoritative substrate
parameter project is https://github.com/parksjm/polyamide-amber-ff .
**The code repo alone does not contain the prepared solvated systems or bulk
trajectory data.** Give the student repository access for updates and choose
**one** curated archive: `core_data.tar.gz`, or `inspection.tar.gz`, which
includes everything in core plus selected QM/MM movies. Either archive has
a self-contained source snapshot. The full original 100 ns trajectories remain
separately addressable in the source manifest and are not included.

The core data include the pinned PA6 release, prepared WT systems/acceptance
receipts, source preparation recipes, MD inputs, prior geometric tables/SI
figures, all eight historical relaxed mutant complexes and provenance.
Those mutant complexes retain their legacy parameter lineage and are pose
references; rebuild and validate their new topologies before MD. The inspection data add selected portable VMD bundles
and QM/MM summaries. Do not copy an entire workstation output tree.
No executable Amber/Gaussian installation or account credential is part of
this handoff; obtain software through its official distribution.

See `DATA_ACCESS.md` for the exact archive transfer command and optional
large-data boundaries. `STUDENT_START_PROMPT.md` can start the student's agent
session. Historical occupancy/persistence summaries are in
`evidence/legacy_md_results/`; preserve their original definitions and FF caveats.

## First concrete tasks

- Verify the archive and all delivered SHA-256 values.
- Confirm WT sequence, mature four-alpha/four-beta assembly and residue map.
- Reproduce the six existing WT occupancy counts from cached evidence.
- Validate one delivered WT system with `tools/run_prepared_md.py plan` on the
  target computer. Plan mode writes validation/planning records but launches no MD.
- After the student/PI chooses hardware and independent seeds, run matched
  WT controls before rebuilding and running the eight mutant systems.

The portable runner consumes **already prepared and accepted** systems. It
never invents missing coordinates, changes protonation or converts a legacy
mutant topology to a new force field. Mutant system construction is an explicit
upstream work item, followed by the same validation gates as WT.

## Current evidence, stated accurately

ASH274 gives a higher pooled Tier A fraction, but the two replicas are
16.77% and 0%; this is not a reproducible improvement yet. The stricter joint
native-Thr screen does not show the same advantage. The old mutant panel used
AM1-BCC/GAFF substrate parameters; new ASH274 runs must be compared with
matched controls, not attributed solely to protonation by comparison with
that legacy panel.

No MD, modeling or LigandMPNN jobs were launched when this handoff was built.
The existing QM/MM jobs on violet continue independently.
