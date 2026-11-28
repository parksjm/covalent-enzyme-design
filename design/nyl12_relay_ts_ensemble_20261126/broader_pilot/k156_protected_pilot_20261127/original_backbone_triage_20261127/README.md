# Original MD-backbone sequence triage

All 24 candidates in the parent pilot completed generation and matched scoring. The current comparison uses only its 12 original MD-backbone candidates, without regenerating sequences.

`rerank.py` reuses the two original TS-parent cross-scores and retains E–S scores. RFD3-template scores are excluded from this derived ranking. Run with `PYTHONPATH=src uv run --offline python -B <absolute path to rerank.py> <absolute pilot path>` from the project root, using a new output directory because outputs are immutable. Raw score files and their original rankings remain unchanged.

`relay_prot101_06` and `relay_prot101_07` are the leading likelihood diagnostics from the two original generation parents. This is a provisional pair for sequence-diversity and structural review, not a catalytic ranking or experimental nomination. Retain all 101 protections and four-copy sequence ties. Next use a bounded, matched WT/design Boltz comparison, canonicalization, clash screening, Amber relaxation and all-four-site analysis; reserve MD for a few finalists. No predictions or scientific jobs were launched by this reanalysis.
