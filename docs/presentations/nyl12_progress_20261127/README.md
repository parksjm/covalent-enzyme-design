# Nyl12 PA6 progress briefing —27November2026

Current editable presentation: `Nyl12_PA6_progress_20261127.pptx`; PDF preview alongside.

30slides tell the experimental motivation, finite-L4/register limits, validated force field, WT/protonation studies, QM/MM checks, competing acylation pathways, current provisional PMF, structural sequence pilot and November deadline plan. Notes carry definitions, full source paths and limits.

`slides.json` is authoritative content; `build_deck.py` builds editable text/tables around figures. `plot_evidence.py` plots snapshotted arrays without recalculating the PMF. `evidence/` freezes55ps PMF and WT site records. `assets/structures/provenance.json` and renderer scripts preserve the new molecular figures.

Run `bash docs/presentations/nyl12_progress_20261127/rebuild.sh` from the repository. Uses existing offline `uv` environments, pinned python-pptx1.0.2, project Matplotlib, LibreOffice and pdftoppm; adds no dependencies. Prior outputs are preserved in central `clutter`. Date is generated in America/New_York. A later rebuild retains27November scientific evidence unless content is explicitly updated.

`source_manifest.json`, `.build.json` and `MANIFEST.sha256` record provenance; `quality_checks.json` describes checks and limitations. Prior25November deck remains unchanged at its prior logical path.

Scientific correction: the PI-approved orientation has positive phi under the current ordered projected-torsion definition. Earlier negative-only target advice was wrong and withdrawn. Numerical signs and molecular source images were preserved. Active design target is explicit +1 plus A/B geometry after relaxation.
