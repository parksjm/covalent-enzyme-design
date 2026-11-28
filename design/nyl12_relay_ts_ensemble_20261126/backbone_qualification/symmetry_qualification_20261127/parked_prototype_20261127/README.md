# Parked RFD3 symmetry prototype

The PI chose the existing MD-derived TS-like ensemble as the current design baseline on 27 November 2026. No RFD3 GPU qualification or replacement LigandMPNN generation was launched.

The source snapshot, original worktree files and reversible development patch preserve the unfinished implementation. They are not production-qualified. The two actual native CPU preflights failed before feature generation on E233 versus mature beta-local numbering. Fix that mapping and requalify all native inputs before any future execution. Native projection tests alone did not qualify either scientific parent.

For reproducing the exact-projection rejection calculation, use the archived `source_snapshot/src` on PYTHONPATH; the geometry module is retained there. The JSON evidence embeds the analysis source and the source-structure hashes. Its exact assembly projection was rejected because it creates clashes and large tail displacements.

The original export implementation was restored after 52 focused prototype/export tests passed. The prototype snapshot records the abandoned changes; no installed Foundry files or running scientific jobs were modified.
