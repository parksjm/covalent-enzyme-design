# Research repository checkpoint: 2026-11-19

This checkpoint saves the custom mechanism workflow, manuscript drafts,
controlling roadmap and bounded parallel NEB diagnostics developed since the
2026-10-26 commit. The private repository is
`parksjm/covalent-enzyme-design`, branch `codex/initial-import`.

The current execution plan is `NYL12_MECHANISM_ROADMAP.md`; the handoff holds
resume details. The manuscript Methods and Results distinguish numerical
path diagnostics from chemical/force convergence and verified barriers.
The four parallel spacing tests are documented in their existing campaign
README. L-BFGS remains primary; CG is a short diagnostic.

Versioned content includes source/tests, small configurations, references,
selected canonical campaign scripts, and the exporter preparation/calibration
records required by the current code. The supplied `cif/n12c.cif` remains an
unpublished coauthor-provided reference; citation details await publication.
One Model 23 PDB is retained explicitly as a geometry-test fixture.

Raw trajectories, restart/topology collections, Gaussian outputs, generated
inspection histories, bulk models and MCCE/analysis output stay on the
workstations. Their locations and identities remain in run records and the
action journal. This repository checkpoint does not claim those data are
archived remotely or that the working pipeline is a turnkey public release.
The authoritative polyamide parameter project remains separate.

`provenance/checkpoints/20261119/` retains an immutable action-journal snapshot,
a reviewed source-file manifest/allowlist and verification results. The live
append-only journal remains `output_agent_actions/20261021.jsonl`; commit/push
outcomes are appended there after the snapshot, avoiding recursive commits.

A clean machine still needs the declared scientific software and source data.
Some historical orientation analysis imports the older
`enzyme-substrate-modeling` workspace. That limitation is recorded rather
than silently vendoring another project during this checkpoint.
