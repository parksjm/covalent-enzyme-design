# Nyl12–PA6 research progress briefing

Prepared 23 November 2026. Open **Nyl12_PA6_progress_20261123.pptx** for the
editable 40-slide deck, or the matching PDF for convenient viewing.

The deck groups the work by scientific question: experimental baseline,
parameterization, structural sampling and protonation, reaction paths and
free energies, structural design, and the next decisions. The Eyring table
includes 300 K, 65 °C and 75 °C. Four atomic-view panels are explicitly marked
VMD placeholders; they are not fabricated molecular renderings. Their source
bundles are recorded in slide notes.

## Scientific interpretation

- **35.81 kcal/mol is a NEB potential-energy maximum, ΔE, not ΔG‡.**
  The associated path family has already supplied a refined TS candidate
  and the NEB-seeded US campaign. Their different reactant references,
  relaxation and sampling prohibit treating their numerical differences as
  a pure entropy correction.
- The 47.05-kcal/mol PMF maximum is an apparent projected height. Sampling
  convergence and a validated activation free energy remain unestablished.
- Crystal calibration retains every native-Thr threshold. All ten original
  crystal sites pass; known-active control rejection still limits activity
  prediction. Derived symmetry/reference copies are not independent experiments.
- New N-terminal-relay releases are running at deck preparation. Their
  finite initial gradients are accepted; no final chemical outcome is claimed.
- MCCE native completion is separate from pKa validity. Current NylC outputs
  remain unqualified; no protonation states were assigned.
- The current modular PA6 release is MD-ready within its declared
  zwitterionic L4 scope. The dated October parameterization deck was consulted;
  current architecture/profile/MD figures supersede its historical hold status.
  Compact-conformer electrostatics and extension to other states remain limited.

## Sources and rebuilding

`slides.json` is the editable scientific narrative. `asset_sources.json` maps
each copied figure to its original path and SHA-256 digest. Full references
and caveats are also embedded in the PowerPoint speaker notes.
`evidence_inventory.json` records the independent source audit.

From the repository checkout, run `bash docs/presentations/nyl12_progress_20261123/rebuild.sh`.
The builder uses an isolated `uv run --no-project --with python-pptx==1.0.2`
report environment; it does not add project production dependencies. PDF
rendering requires the existing LibreOffice executable. The build JSON records
the source, renderer, PowerPoint and asset hashes.

The editable deck, PDF, JSON sources and figures are included in
`scripts/pull_nyl12_inspection.sh`, retaining the workstation directory layout
on the laptop. Verify the downloaded bundle with `sha256sum -c MANIFEST.sha256`
or `shasum -a 256 -c MANIFEST.sha256`.

Renderer verification: eight focused tests, 99% coverage; all layout variants,
source failures, full source notes, image hashes and slide bounds are checked.
The full deck was rendered through LibreOffice and its 40-page contact sheet
visually inspected. Those checks do not establish scientific validity of a
profile or validate a new mechanism.

Calibration verification used the existing Biopython analysis environment:13 passed,97% new-code coverage. The default project environment selectively skips three CIF-parser tests when this optional analysis prerequisite is absent; scalar and control checks still run. No production dependency was added.

Final repository checkpoint: `uv run pytest -q` — **2,060 passed,17 skipped**. Skips comprise13 retained-preparation-fixture prerequisites, three optional CIF-parser tests, and the optional PowerPoint test module in the default environment. All parser/report checks were also run successfully in their explicit tool environments. Relay preparation/export/monitor checks:18 passed,91% coverage.
