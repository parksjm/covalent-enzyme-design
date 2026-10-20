# Covalent enzyme design

Research pipeline under development for covalent enzyme–substrate models,
chemical-state validation, parameter caching, and QM/MM geometry refinement.
The current case study is Nyl12: PA6 L4 cleavage yielding L2 and poly(X,AMA)
cleavage yielding MAXA.

## Current status

Covalent AF3 TI1 predictions and Amber topologies exist in external project
storage. The repository contains the scripts and provenance needed to trace
their construction. **No converged TI1 QM/MM minimum is established yet.**
The rigid-pocket relay-water searches failed their contact checks. Neutral
nonreacting beta termini still use provisional charges; these systems are
not ready for classical production MD or mechanistic barrier comparisons.

The import namespace remains `nyl12_covalent` to preserve existing entry points.
This initial extraction is a development snapshot, not a portable release:
some builders still have workstation-specific defaults and Nyl12 residue
indices. Read the resume note and pass explicit input/output/tool paths.
Active calculations and large artifacts remain in the original project.

Start with [the standalone handoff](docs/standalone_handoff_20261020.md).

## Environment and checks

Use Python 3.12 and `uv`. The project manifest records the Python dependency
versions already used in the source environment.

```bash
uv run pytest
uv run ruff check src scripts tests
```

Verification on import: **51 tests pass** and Ruff passes. Overall line coverage
is **38%** because the imported chemistry/topology/prediction code still needs
tests. The new workflow modules and wrappers each exceed 80%; this is not a
claim of complete pipeline validation. See `provenance/coverage_20261020.json`.

The Python dependencies do not install external engines. Geometry/topology
work additionally requires the configured AmberTools26 tools and parameter
data. Prediction requires a separately installed, authorized AF3 environment
and model data. Gaussian/RESP is a planned charge-fitting route.
Keep external executable paths and scratch locations explicit.

## Entry points

Run scripts with `--help` for their arguments:

- `scripts/prepare_nyl12_af3_ti1.py`: custom ligand CCD and covalent AF3 input.
- `scripts/canonicalize_nyl12_af3_ti1.py`: inspectable models and atom maps.
- `scripts/build_nyl12_ti1_qmmm.py`: TI1 topology/protonation construction.
- `scripts/search_nyl12_relay_water.py`: rigid-pocket water candidate search.
- `scripts/diagnose_nyl12_relay_water_contacts.py`: explicit contact diagnostics.
- `scripts/audit_nyl12_ti1_qm_region.py`: proposed QM-boundary selection audit.
- `scripts/run_nyl12_job.py`: persistent tmux jobs with execution records.

The older acyl-state entry point remains for provenance. TI1 is the current
priority. Existing QM/MM runners do not implement checkpoint continuation;
never rerun into an existing result directory to resume a calculation.

## Provenance and storage

`provenance/source_inventory.json` identifies the exact imported source files,
their hashes, original Git tracking status, and original workspace revision.
Untracked and modified source files are not reproduced by that revision alone.
The initial commit preserves the reviewed extraction. Subsequent development
belongs here once the handoff explicitly switches the active workflow.

Keep model ensembles, trajectories, parameter caches, scratch, manuscript
drafts, external software, and model weights outside this repository. Reference
them by explicit paths and hashes in calculation manifests. Never interpret a
zero process exit status as scientific convergence.

Workflow skills are versioned under `tools/skills/`; they call the maintained
pipeline scripts. Changing scientific results and job state belong in project
notes, not in skill instructions.

## Next milestones

1. Repair and validate TI1 geometry, including relay-water placement.
2. Validate neutral-terminal charge models and state-specific parameters.
3. Replace workstation defaults with explicit configuration.
4. Expand portable chemistry/topology tests and document an integration example.
5. Add reaction-path and RFdiffusion3/LigandMPNN workflows after validated models.
