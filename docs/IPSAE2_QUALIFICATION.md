# ipSAE2 installation and Boltz scoring qualification

Updated 28 November 2026. The PI authorized an isolated installation and
qualification on WT plus one design before integration. This record supplements
the [initial source assessment](ipsae2_code_review_20261127.md).

## Installed runtime

- Official source: <https://github.com/DunbrackLab/IPSAE2>.
- Package 0.1.11; exact commit `fe560f4c2461eeb91aebded61ac8c1e83b264f60`.
- Environment: `environments/ipsae2_20261128/.venv`, Python 3.12.11.
- Recreate from the committed `pyproject.toml` and `uv.lock` with
  `uv sync --project environments/ipsae2_20261128 --frozen`.
- Executable: `environments/ipsae2_20261128/.venv/bin/ipsae`.
- The project and Boltz dependency environments were not modified. Scoring is
  CPU-only with numerical-library threads limited to one and GPU visibility empty.

The environment lock is portable; runtime acceptance hashes and absolute paths
are host-specific. Recreating the environment requires verifying its executable
and source hashes and issuing a new technical-qualification receipt, rather than
reusing another installation's runtime receipt.

## Accepted checks

Evidence resides in `analysis/ipsae2_qualification_20261128/`:

| Check | Outcome |
|---|---|
| Upstream tests and bundled example | 53 tests passed; example values matched |
| Real WT/design inputs | `WT_seed_20261126`, `relay_prot101_06_seed_20261126` |
| Native matrix shape | 1444 × 1444 for each: 4 × (232 + 96) + 4 × 33 tokens |
| Source identity | Native token asymmetry IDs, CIF atom order and token pLDDT verified |
| Identical-copy negative controls | Swapped subunit/ligand blocks rejected |
| PAE axes | Native row-alignment/column-scored orientation retained |
| Independent numerical formula | 448 directional values per model; maximum error <5 × 10⁻⁷ |
| Repeated chain-pair scoring | Byte-identical TSV outputs |
| Native files before/after | Coordinate, PAE, confidence, pLDDT and structure hashes unchanged |

The numerical comparisons cover original chain pairs, whole-protein grouping
and four mature-protomer groups. They validate implementation and input handling,
not a calibrated Nyl12 activity predictor. A subsequent finite-input hardening
check also rejects nonfinite CIF confidence values; both qualified actual files
already contained finite values.

Biological pairing must be determined per model. Protein-only correspondence to
the native assembly gives WT A–E/B–F/C–G/D–H, whereas design 06 uses
A–H/B–G/C–F/D–E in native label IDs. Fixed native letter pairs would therefore
misgroup the design. These native IDs precede canonicalization and must not be
confused with the canonical chain identities used for Amber and geometry.

## Integration and interpretation

`scripts/score_boltz_ipsae2.py` runs optional postprocessing separately from
the frozen Boltz/Amber controller. Accepted global WT/design qualification and
fresh per-model mapping/hash receipts are prerequisites. Retain exact commands,
logs, raw TSVs, parsed JSON and explicit failure records. A scoring failure must
not stop or alter prediction, Amber relaxation, QM/MM or US jobs.

Score unmodified native predictions with matching native PAE. Do not supply
canonicalized or relaxed coordinates against the original PAE. Retain original
protein interfaces, correctly paired mature protomers and each of the four
ligands separately. Use explicit protein-aligned `ipsae12` or `ipsae21` for
ligand comparisons; grouped scalar `ipsae` takes a maximum that can select the
reverse direction. Keep both directions, native ipTM, n0, d0 and contact fields.
An absent PAE-selected contact is not proof of absent geometric contact.

Whole-protein, protein-aligned scores in the qualification seed are:

| Native ligand | WT | Design 06 |
|---|---:|---:|
| I | 0.026397 | 0.017145 |
| J | 0.016310 | 0.011036 |
| K | 0.014456 | 0.009317 |
| L | 0.018853 | 0.013247 |

All 33 heavy atoms of each PA6 L4 enter this direction with d0 = 2 Å; no
ligand PAE cutoff removes uncertain distal atoms. These low scores describe
whole-ligand prediction uncertainty. They neither establish an incorrect
reactive-amide geometry nor show that one sequence is less active. Identical
ligand letters across predictions do not themselves match biological sites;
site comparisons require the explicit structural mapping. The four copies are
not independent prediction replicates, and this table is one seed only.

Keep ipSAE2 supplementary to positive-phi Tier A/B, clashes, assembly and
chemical checks. No universal 0.8 or PA6-specific cutoff is adopted. High pose
confidence is not binding affinity, transition-state stabilization, a reduced
barrier or improved hydrolysis. Whole flexible-ligand scores and amide-specific
geometry answer different questions.

## Production execution and checks

**Completed:** all six matched models scored successfully; the postprocessing
job exited 0. All 24 whole-protein/ligand rows are retained in
`analysis/ipsae2_qualification_20261128/matched_ligand_confidence.csv`;
`production_acceptance.json` confirms unchanged native-artifact hashes.
The full production config, per-case mappings, physical-pairing evidence, score
JSON/TSVs and process records are archived under `evidence/production_receipts.tar.gz`
with an independent member-hash manifest. Original server files remain in place.

The six-case matched campaign is configured in
`design/nyl12_relay_ts_ensemble_20261126/broader_pilot/original_backbone_matched_20261127/scoring/ipsae2/sidecar_config.json`.
Session `nyl12_ipsae2_score_20261128` processes only completed native outputs,
retains per-case failures and exits when all six cases are terminal. Its
`jobs/scoring01/status.json` records process completion separately from the
per-model scientific/input acceptance and scoring status. A one-thread, no-GPU
environment and CPU niceness 10 limit interference with ongoing calculations.

The reusable preparation entry point is `scripts/prepare_boltz_ipsae2.py`;
it requires six completed native cases, explicit runtime/qualification paths,
the native Boltz Python and the existing protein-only pairing adapter. It
produces fresh hash-bound per-model manifests and mappings, then the scorer
can run with `scripts/score_boltz_ipsae2.py --config PATH --watch`. Preparation
is deliberately separate from frozen prediction/Amber inputs.

All 38 focused project tests passed, covering mapping, preparation, scoring
and terminal-state handling. Combined coverage is 89%, with every module at
least 80%; Ruff lint and formatting passed. These tests supplement the 53
upstream tests and real-model numerical qualification.

## Retained evidence

- `environment/`: installation metadata, exact commands, upstream tests/example.
- `runtime.json`: executable and installed-source hashes.
- `*_native_mapping.json`: strict native-order receipts; exact compressed copies
  are versioned under `evidence/`, bound by `evidence_manifest.json`.
- `physical_pairing.json`: independently verified alpha/beta pairing.
- `run_qualification.py`: independent numerical checks and immutable-source audit.
- `numerical_attempt01/qualification_acceptance.json`: accepted WT/design gate.
- Per-model `numerical_attempt01/`: raw score TSVs, process logs and parsed JSON.
- Matched campaign `scoring/ipsae2/`: production config, receipts and scores.

The locked installation, mapping verifier, scorer and configured postprocessing
entry point are reusable. A different predictor representation, chain layout or
ligand tokenization requires explicit input qualification; it must not silently
reuse this Nyl12-specific chain/token assumption.
