# Local MSA workflow on Violet — verified 27 November 2026

Read this before proposing a new installation or a public-server request.
The working setup was recovered from the enzyme-substrate-modeling project.
The earlier conclusion that local search was unavailable checked the wrong
environment/database locations and is superseded by a successful native run.

## Validated tools and data

- ColabFold executables:
  `/data/jpc/enzyme-substrate-modeling/.pixi/envs/colabfold/bin/`.
  Its `colabfold_search --help` works; MMseqs reports 18.8cc5c.
- Indexed, GPU-ready databases: `/data/boltz_db/uniref30_2302_db` and
  `/data/boltz_db/colabfold_envdb_202108_db`.
- LightningBoltz checkout: `/data/jpc/code/lightning-boltz`, with environment
  `/data/jpc/code/lboltz-env`. Its local-search option is separate from the
  public MSA-server option. A broken global wrapper is not evidence that the
  Pixi installation is broken.

The successful search-only command uses the full executable path:

```bash
PATH=/data/jpc/enzyme-substrate-modeling/.pixi/envs/colabfold/bin:$PATH \
CUDA_VISIBLE_DEVICES=GPU-ec7e1af2-0974-2893-d34b-411eee82bbd6 \
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=2 \
/data/jpc/enzyme-substrate-modeling/.pixi/envs/colabfold/bin/colabfold_search \
  --gpu 1 --prefilter-mode 1 --threads 2 --use-env 1 --use-templates 0 \
  QUERY.fasta /data/boltz_db NEW_OUTPUT_DIRECTORY
```

`--gpu 1` enables GPU search; CUDA_VISIBLE_DEVICES selects the device. Verify
availability before launch. The retained runner waits for an idle GPU0 and
checks memory/disk prerequisites. Do not interrupt existing jobs. Native
ColabFold names this output from the FASTA query header, not necessarily
`0.a3m`. No packages or databases were installed/downloaded for this run.

## Precursor-derived paired mature-chain inputs

The WT search query is one exact 328-residue sequence: alpha 232 plus beta 96.
This search representation does not change the mature eight-chain structure
used in prediction. Split homolog rows at the query-aligned boundary232|233,
not at a raw text offset. Preserve insertions/gaps and same-precursor source
rows; assign explicit shared integer pairing keys. Boundary insertions attach
to beta consistently with native deletion counts. Exclude duplicated segments
symmetrically so native per-chain deduplication cannot create orphan pairs.

Do not pair the older independently generated alpha/beta MSAs by matching
synthetic `homolog_NNNNN` labels. Cache one search per exact unique query and
database/settings version; equivalent copies and prediction seeds reuse it.

The successful source had 2,796 records, including the query. After symmetric
native-compatible filtering it supplied 1,387 homolog pairs; 1,408 rows were
excluded. Both WT controls passed the actual native dataset check: all eight
protein chains carry coherent same-precursor pairing, record.templates is
empty, and every template feature equals the stock dummy-template output.
Categorical dummy features are not necessarily numerically zero.

## Retained execution and failure evidence

Campaign: `design/nyl12_relay_ts_ensemble_20261126/broader_pilot/`
`wt_paired_precursor_controls_20261127/`.

`precursor_search/` retains the exact FASTA, database metadata, runnable
scripts, source-code snapshots and logs. Native search completed with exit 0;
the first wrapper failed after that because it expected `0.a3m`. The separate
`accept_completed_search.py` validated the existing query-named result without
rerunning the search. The original failed wrapper receipt remains unchanged;
`source_acceptance.json` and `source_acceptance_job01/` establish recovery.
The future search entry point and regression tests use native query naming.

Native prediction is stock Boltz, with fresh output directories and no crystal
templates, substrate coordinates or contact restraints. The two matched seeds
are 20261126/20261127; four zwitterionic PA6 L4 molecules and eight mature protein
chains remain. Canonical all 48 geometry audits, identical Amber relaxation,
positive-phi gates and curated comparison are queued in persistent tmux jobs.
Input acceptance is not evidence that pocket geometry or activity improved.

Use public ColabFold only as an explicitly recorded alternative if the local
route becomes infeasible; preserve failure evidence rather than silently
changing the source. Local search is currently proven operational.

## Canonical export and relaxation checkpoint

Both native predictions finished with exit 0. The full-precision CIF geometry
audit is retained separately from PDB serialization. A 0.65-degree dihedral
change for one distant, rejected pairing arose from three-decimal PDB rounding;
it was not an atom-map discrepancy. The corrected audit checks exact native
versus full-precision CIF geometry and the actual serialized-PDB → sequential
Amber mapping → restored-PDB roundtrip, without loosening mapping tolerances.
All 48 raw-to-PDB metric, tier and assignment deltas remain explicit. No
selected tier or assignment changed in either WT case. Failed receipts and
first-attempt code/derived artifacts are retained in
`canonical_precision_attempt01_archive/`; recovery source hashes are bound in
`canonical_precision_recovery_binding.json`. Native inference was not rerun.

Persistent recovery sessions are `nyl12_paired_WT_collect_precision01` and
`nyl12_paired_WT_amber_precision01`. The unrelaxed seed 20261126 has one Tier A
and three Tier B positive-phi sites; seed 20261127 has no passing sites. These
are preliminary geometry observations, not relaxation or activity acceptance.
