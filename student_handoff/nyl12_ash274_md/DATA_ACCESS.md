# Delivery and optional additional data

The public `covalent-enzyme-design` GitHub repository contains code and
documents. The separate parameter repository retains its prior visibility. The prepared solvated systems, restart files and inspection movies
need the separately generated archive. Choose **one**:

- `core_data.tar.gz`: complete student brief, source snapshot, PA6 parameter
  release, three prepared WT microstates, historical mutant starting poses,
  prior analysis tables and original input/provenance records.
- `inspection.tar.gz`: all core contents plus selected QM/MM VMD bundles,
  structures, plots and summaries. It is a superset of core.

Both archives and their individual `.sha256` files are published on violet at:

`/data/jpc/covalent-enzyme-design/output_student_handoff/nyl12_ash274_md_20261119/`

After cloning the code repository, run its fetcher with an authorized
SSH alias/account. This command transfers one archive, verifies it and extracts
it; it does not submit workstation jobs:

```bash
python student_handoff/nyl12_ash274_md/tools/fetch_bundle.py \
  --source YOUR_AUTHORIZED_VIOLET_ALIAS \
  --remote-dir /data/jpc/covalent-enzyme-design/output_student_handoff/nyl12_ash274_md_20261119 \
  --package inspection --destination ./nyl12_student_data --extract
```

If copying the fetcher without cloning, also copy its adjacent
`student_bundle_support.py`. Alternatively transfer the chosen archive and
checksum by an approved institutional file-transfer method. Repository access
and workstation SSH/data access are separate; neither is granted by these docs.
No invitations, account changes or messages were sent while building this handoff.

The archive's `export_manifest.json` gives exact original workstation sources,
hashes and portable destinations. `CONTENTS.sha256` covers delivered files.
Historical absolute paths in provenance are evidence; use portable dataset
manifests/configs for execution on the student's workstation.

The approximately 75.55 GB of six WT raw production trajectories are excluded.
Original paths are catalogued as `raw_reference_only` in the repository
`data_manifest.json`. Full-frame analysis from
`analysis/raw_comparison_manifest.json` requires obtaining those optional
trajectories and placing them at its declared relative paths. Cached completed
analyses are included and can be inspected without raw trajectories.

The full QM/MM workspace remains on violet under
`/data/jpc/covalent-enzyme-design/qmmm/`. This inspection archive deliberately
contains selected immutable published movie generations and compact results,
not every old path, failed calculation or live update. The supplied source
snapshot includes the broader manuscript, reference/source ledger, mechanism
roadmap and run scripts. The PI can authorize targeted access to further
QM/MM directories or legacy raw trajectories when needed; preserve associated
input manifests, topology/atom maps and hashes when copying them.

`n12c.cif` is the PI's reference from a coauthored manuscript in revision;
its final citation is pending. Its tracked copy became publicly visible with
the PI's explicit repository-visibility request on 2026-11-19. This release
does not publish or grant access to all raw workstation data. Software engines,
credentials, Gaussian scratch/checkpoints and bulk MD trajectories are not
distributed here.
