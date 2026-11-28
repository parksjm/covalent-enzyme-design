# Validation checkpoint — 27 November 2026

The release check `uv run --offline pytest -q --tb=short` completed with
3,223 passed, 17 skipped and three failures in 516.58 s. All three failures
were test fixtures, not weakened production acceptance gates:

- One hydration report test read a live source manifest whose recorded
  status-file hash was stale. It now uses an isolated producer-shaped fixture
  and explicitly verifies rejection after source mutation. The live manifest
  was not resealed or changed. Before any future hydration follow-on relying
  on that historical report, regenerate/qualify its derived analysis from
  accepted immutable sources; the original hash gate still rejects stale data.
- Two laptop-pull tests mocked every optional SSH probe as successful but
  returned an empty immutable publication path. Production correctly rejected
  this. The mocks now declare absent optional bundles explicitly; no pull or
  checksum guard was relaxed.

After these test-only changes, both affected complete files passed:
`uv run --offline pytest -q tests/test_hydration_followon.py tests/test_nyl12_inspection_download.py`
→ 38 passed. Restored protection/broad/ensemble/score-resume tests separately
passed 45 cases. Ruff checks and formatting passed for changed tests/policy.
No second full-suite run was claimed or required to rerun already passing
production tests after these isolated fixture repairs.

Additional executed, source-bound validations include:

- Azalea 100 ps continuation: 22 focused tests, helper 93% coverage. All 57
  windows entered new native Sander production from accepted 60 ps coordinate,
  velocity, time and box states. One-time analysis discard is unchanged.
- Replacement 101-position LigandMPNN run: 7 namespace tests, 100% coverage;
  native four-TS/E–S parsers and GPU prerequisite passed. Generation and
  scoring use the same frozen policy package.
- Existing score-only recovery helper: 9 tests, 81% coverage; historical
  K156-variable launch wrapper: 9 tests, 88% coverage. Their generated
  candidates remain ineligible under the restored policy; completed outputs
  and cancellation evidence are retained.
- Supplemental movie publication: 48 focused tests; source coordinate/box
  parity for 1,400 frames, manifest/time/map checks and VMD cell verification.

Raw trajectories, licensed native executables, downloadable model weights
and workstation caches are not source-code dependencies to upload blindly.
The repository records scripts, parameter/model versions, curated small input
structures, manifests, reports and figures; full data/source snapshots remain
at their explicit recorded workstation paths and in inspection exports.
The compressed action-journal checkpoint records this release boundary; later
Git commit/push outcomes continue in the live append-only journal.
