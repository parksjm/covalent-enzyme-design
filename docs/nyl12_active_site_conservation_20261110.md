# Nyl12 active-site conservation, 2026-11-10

The existing Boltz alignments strongly conserve Asp272 and Asp274. This supports
their functional or structural importance; it does not determine their pKa,
protonation, proton-transfer direction, or catalytic role in Nyl12.

This is an analysis of **existing local Boltz MSAs**, not a new BLAST search.
No sequence was submitted to an external search service. The source CSVs have
sequence and pairing-key columns but lack accession/taxonomy annotations; these
are broad-family sequence candidates rather than confirmed nylon hydrolases.

| Nyl12 position | Public NylC 3AXG correspondence | Query residue / non-gap calls | Conservation | Gap fraction |
|---|---|---:|---:|---:|
| Phe113 | Tyr146 | 718 / 1268 | 56.62% Phe; 95.90% Phe/Tyr | 0% |
| Lys156 | Lys189 | 876 / 1226 | 71.45% | 3.31% |
| Asn186 | Asn219 | 1263 / 1268 | 99.61% | 0% |
| Thr233 | Thr267 | 2113 / 2113 | 100% | 0.61% |
| Asp272 | Asp306 | 2121 / 2126 | 99.76% | 0% |
| Asp274 | Asp308 | 2125 / 2126 | 99.95% | 0% |

Both Asp positions are Asp in 2120/2126 retained beta sequences. Their query
motif is `YATILDGDVLYTV`; the corresponding public NylC motif includes
`FHTDMDGDTLFAV`. Numbering correspondence was independently checked by global
alignment of the modeled Nyl12 alpha+beta query against the
[public 3AXG precursor FASTA](https://www.rcsb.org/fasta/entry/3AXG/display).
The full aligned sequences are retained in the analysis artifacts.

The alpha MSA contains 3137 raw rows and 3124 unique aligned sequences; 1268
nonquery sequences pass the filter. The beta MSA contains 4913 raw rows and
4859 unique aligned sequences; 2126 pass. Query lengths are 232 and 96 residues.
Lowercase insertion characters are removed; deletion gaps retain their columns.
Exact aligned sequences are deduplicated, all exact query matches are excluded,
and at least 80% of query match columns must be occupied. Coverage counts
non-gap calls. Frequencies use non-gap denominators; separate standard-residue
denominators and ambiguous-X counts are also retained in JSON. The two chain
MSAs are analyzed independently. No phylogenetic or identity-cluster weighting
was performed, so these counts are not an effective independent sample size.

The [2023 NylC study](https://doi.org/10.1111/febs.16755) proposes an
Asp308–Asp306–Thr267 network for substrate hydrolysis and studies its role in
autocleavage. The numbering alignment and conservation make this network a
specific hypothesis to test in Nyl12; they do not establish that the entire
proposed NylC mechanism transfers to Nyl12. Phe113 replacing NylC Tyr146 is an
explicit chemical difference. Conservation can also reflect maturation or
folding requirements.

Reproduce with `uv run python scripts/analyze_nyl12_site_conservation.py`,
supplying the two original MSA CSVs and retained public FASTA as recorded in
`command_argv` in
`analysis/nyl12_active_site_conservation_20261110/conservation.json`; choose a
new output directory. The original temporary FASTA argument can be replaced
with the retained `public_nylc_3axg.fasta`. The output records source/script
SHA256 hashes, Python version, filters, raw/unique/filtered counts, numbering
mapping and public reference origin. Existing input CSVs remain unchanged.

Verification: five focused tests pass, including insertion handling, gap/X
denominators, deduplication, coverage, numbering alignment, source retention and
overwrite refusal; script statement coverage is 94%. Ruff format/check pass.

## QM-region additions assessed from actual coordinates

The scripted read-only audit `scripts/audit_nyl12_qm_neighbors.py` records atom
membership and distances for E–S, the restrained attack09 TI1-like seed and the
converged unbiased acylenzyme. Results and input/source hashes are in
`analysis/nyl12_active_site_conservation_20261110/{es,ti1,acyl}_qm_neighbors.json`.
No running QM region or protonation state was changed.

| Candidate | E–S | TI1-like | Acylenzyme | Interpretation |
|---|---:|---:|---:|---|
| Lys156 NZ–Asn186 OD1 (Å) | 3.237 | 3.120 | 3.143 | Persistent second-shell contact; useful QM-region sensitivity test |
| Gly273 backbone N–Thr233 OG1 (Å) | 2.985 | 2.926 | 3.275 | Nearby backbone donor, currently MM-only |

The TI1-like Gly273 H–OG1 distance is1.919Å. Its proximity is a reason to test
quantum treatment of this peptide donor; no extra Gly conservation statistic
was calculated here. Lys156 has71.45% Lys among non-gap calls in this MSA and
maps to NylC Lys189, discussed in the NylC catalytic/substrate-binding network.
The Lys–Asn contact can affect polarization; it does not establish proton transfer
through Asn or prove that Lys is a general base in Nyl12. Keep the existing Lys
protonation for a first sensitivity test. Transferring protonated Lys into QM
changes the QM formal charge and boundary charge treatment; audit both explicitly.

First proposed expansion: Lys156 side chain plus an appropriately bounded
peptide segment containing the Gly273 donor and adjacent Asp272/Asp274 backbone.
Preserve intact amide groups; do not cut through a peptide C–N bond or add an
isolated NH donor. Select exact caps/boundaries only during implementation.
The already included Asp sidechains do not quantum-mechanically represent this
backbone donor. Leu271 carbonyl (3.514Å from leaving N in attack09) is a secondary
candidate if it forms a persistent relevant hydrogen bond or belongs to the
chosen intact peptide fragment.

Key atoms of Phe113, Asn186, Thr233, Thr234, Asp272 and Asp274 are already QM.
The closest Arg159 guanidinium contact to the five reactive heavy atoms is about
19.97Å in attack09; all four sequence-equivalent Arg copies were checked. There
is no current proximity-based reason to add Arg159 to QM, despite its reported
substrate-binding/product-release role. Do not equate a residue's functional
importance with a requirement for QM treatment.

No current one-water relay is established. A water becomes a priority QM addition
when a physically supported relay or deacylation geometry is selected. MM water
cannot exchange protons. For matched comparison, preserve identical QM membership
across both endpoints and every image of a new path; do not switch it mid-NEB.
First compare selected existing E–S, acylenzyme and path geometries with the
expanded region; no guarantee that expansion fixes the barrier is implied.

Sources: [NylC2023](https://doi.org/10.1111/febs.16755),
[Nyl12 structural study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13015463/).
