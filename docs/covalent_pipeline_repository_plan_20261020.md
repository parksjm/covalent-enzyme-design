# Standalone covalent enzyme design pipeline

Status: private repository created on 2026-10-20 at
`https://github.com/parksjm/covalent-enzyme-design`; preparing the initial source
snapshot. No scientific output files are being moved.

## Recommendation and boundary

Create `parksjm/covalent-enzyme-design`, initially private. The user selected
their personal account, and both the GitHub connector and CLI identify the
authenticated account as `parksjm`. The original organization-owner suggestion
was superseded by that selection. Local checkout:
`/data/jpc/covalent-enzyme-design`.

The first release should describe its actual scope: covalent intermediate
construction, parameter provenance, and QM/MM preparation/refinement. Reaction
paths and RFdiffusion3/LigandMPNN design are planned extensions, not implemented
capabilities. Keep Nyl12 as an explicit case study with PA6 L2 and poly(X,AMA)
MAXA product gates; generic workflows must not silently inherit these choices.

Include:

- `src/nyl12_covalent/`, initially preserving the import namespace and entry points.
- Covalent builder, AF3 preparation/canonicalization, job runner, relay-water
  search/diagnostic entry points and their tests.
- Versioned scientific workflow documentation, skill sources under
  `tools/skills/`, and small redistributable example inputs.
- Explicit software/environment configuration and a `uv` project definition.
  Resolve dependencies from the already used software; adding production
  dependencies still requires the user's confirmation.
- Artifact manifests that identify source inputs, checksums, software versions,
  command lines, seeds, charge methods, and acceptance results.

Keep outside Git:

- Large model ensembles, topologies, trajectories, job outputs, parameter caches,
  Gaussian scratch, local Python environments, and external software/model weights.
- Manuscript drafts and attached papers. Reference appropriate sources and
  access-controlled artifacts without copying them into the software repository.

## Concrete extraction work

1. Snapshot an explicit allowlist from this dirty workspace; retain source paths,
   source revision, hashes, and the fact that much of this code is currently
   untracked. Do not claim the current Git revision alone reproduces these files.
2. Move workstation-specific AF3, AmberTools, input structure/MSA, and output paths
   into explicit configuration. Keep scientific defaults separate from host paths.
3. Preserve the working import namespace and CLI arguments during extraction;
   add installed entry points without breaking existing scripted procedures.
4. Keep `jobs.py` and `parameter_cache.py` reusable. Isolate Nyl12 residue/chain
   mappings, substrate definitions, and protonation choices in the case-study
   layer. Current inspection and QM-region code contain Nyl12-specific indices.
5. Run the portable unit checks and a small documented integration example in
   the extracted checkout before switching active calculations to it. Real
   engine checks must report external prerequisites and scientific acceptance,
   not merely process exit status.
6. Review the exact staged files, origin metadata, exclusions, and remote
   owner/visibility before the initial push. Preserve current output paths and
   in-progress calculations; make migration a recorded change, not a silent move.

## Current evidence and readiness

The audited Python dependencies include NumPy, RDKit, and ParmEd; external tools
include AmberTools26 and, for prediction, AF3. Gaussian is optional for the
planned RESP route. Gaussian 16 C.02 is now accessible on violet and passed a
scripted tmux preflight; use `/data/scratch/gaussian_scratch/jpc`. The long-lived
Codex process needs `sg g16` to acquire the newly granted group. The monolithic
environment should not be copied wholesale.

On 2026-10-20, 14 focused tests passed for job execution, canonical inspection,
and parameter caching. This does not establish chemistry/topology coverage or
validate a QM/MM minimum. Both AF3 TI1 models still need geometry refinement.
The current rigid-pocket relay-water search has no passing candidates.

## Skills

- `scientific-tmux-jobs`: source created in `tools/skills/`; captures persistent
  execution and provenance. It calls project code instead of duplicating it.
- Covalent-intermediate modeling skill: create after the TI1 refinement workflow
  has passed its geometry and protonation checks. Record decision rules, state
  definitions, and validation steps; keep model-specific results in project docs.
- Consider additional skills after repeated use demonstrates a benefit:
  validated parameter/RESP cache workflows, canonical model inspection, and
  manuscript methods/provenance assembly. Prefer supporting references in the
  two main skills until a workflow warrants independent discovery.
