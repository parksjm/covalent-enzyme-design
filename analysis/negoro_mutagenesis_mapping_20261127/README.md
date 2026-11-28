# Gap-aware NylC–Nyl12 protection correspondence

The JSON is authoritative; CSV repeats its 355 NylC-reference rows. Both count
nongap residues independently in the existing alignment. The complete Nyl12
query was verified against all 328 WT identities in the current protection
inventory. The reference was verified against the preserved 355-residue public
3AXG FASTA. Source SHA-256 hashes are recorded in JSON. Known critical-site
mappings were checked explicitly. No policy or figure was changed.

**Important:** this reference is NylCA/3AXG, not necessarily the same variant as
NylCp2 used in a given paper. For example, reference residues G122/Y130 should
not be confused with parental p2 D122/H130. Establish each paper's construct
numbering and identity before transferring evidence. The existing alignment is
a simple global sequence alignment (match 2, mismatch −1, gap −2), not a new
structure-validated residue correspondence. A query gap is not evidence that a
nearby residue is the functional equivalent. Local low-identity loops need
additional scrutiny. There are 39 NylC positions aligned to a query gap.

| NylC reference | Nyl12 | Current protected status |
|---|---|---|
| D36 | T21 | Variable |
| I75 | Query gap | Unassigned |
| L106 | F72 | Variable |
| S111 | I77 | Interface |
| G122 | S89 | Variable |
| Y130 | N97 | Interface |
| F134 | F101 | Interface |
| A137 | M104 | Interface |
| S141 | S108 | Variable |
| Y146 | F113 | Core + interface |
| P157 | P124 | Variable |
| A179 | H146 | Variable |
| K189 | K156 | Core |
| D191 | L158 | Variable |
| N219 | N186 | Core |
| M225 | L192 | Variable |
| N235 | L202 | Variable |
| E263 | E229 | Core |
| N266 | N232 | Core; terminal modeled alpha residue |
| T267 | T233 | Core; first modeled beta residue |
| H295 | S261 | Interface |
| Q299 | H265 | Interface |
| D306 | D272 | Core |
| G307 | G273 | Core |
| D308 | D274 | Core |

Physical alpha chains A–D cover canonical 1–232; beta chains E–H cover canonical
233–328 (local 1–96). These mappings describe sequence homology, not proof that
the catalytic geometry or biochemical effects are identical.

## Scripted exact reproduction

From the repository root, using the existing uv environment (standard library
only; no new dependency):

```bash
uv run python analysis/negoro_mutagenesis_mapping_20261127/build_mapping.py \
  --output /tmp/nylc_mapping_reproduction_20261127_final
uv run ruff check analysis/negoro_mutagenesis_mapping_20261127/build_mapping.py
```

The output must be a new directory. The preserved alignment is reused exactly;
no alignment is recomputed. Its original settings were match 2, mismatch −1,
gap −2. The CLI verifies full equality with the preserved authoritative JSON,
including every row and source hash, plus nine source-specific anchors. It
writes the same JSON/CSV and a `reproduction.json` containing Python version,
script/source/comparison hashes and validation results. The initial scripted
verification succeeded on 27 November 2026 in the project's uv environment.
Source input paths and hashes are explicit in the authoritative JSON and source.
