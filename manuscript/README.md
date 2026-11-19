# Covalent enzyme design manuscript working draft

Started 2026-10-20 at the user's request. This is a development manuscript,
not a report of validated catalytic improvement. Preserve this distinction
when converting planned methods to past tense.

- [Introduction](INTRODUCTION_DRAFT.md)
- [Mechanism execution roadmap and paper scope](../docs/NYL12_MECHANISM_ROADMAP.md)
- [Methods: completed work and planned extensions](METHODS_DRAFT.md)
- The force-field parameterization section in the Methods draft covers the
  planned ff19SB/GAFF2/RESP workflow, state definitions and validation gates.
  It is included with the bibliography in the zero-argument laptop pull.
- [References and verification status](REFERENCES.md)
- [Local source and claim ledger](SOURCE_LEDGER.md)
- Machine-readable source hashes: `../provenance/manuscript_sources_20261020.json`.

The original Nyl12 manuscript files remain untouched. This draft adapts their
product/register framing and methodological limitations, and adds the current
covalent-state workflow. Historical Boltz/MD results are not new results of
this repository. Experimental activity claims require the original assay
table, product definition, normalization and replicate structure before reuse.

Update protocol: retain dated source/decision records; tie each numerical claim
to an artifact and its hash; record failures as well as successes; state whether
a method is implemented, executed, converged, scientifically accepted, or only
planned. Store evolving calculation state in new dated provenance snapshots.
Do not overwrite past snapshots to make them agree with later job state.

Student protonation/occupancy study summary and matched comparison plan:
`student_handoff/nyl12_ash274_md/STUDY_SUMMARY.md` and `STUDY_PLAN.md` in the
repository root. These remain exploratory, separate from verified kinetic claims.
