# ipSAE2 bounded implementation review — 27 November 2026

## Revision and scope

Primary source: https://github.com/DunbrackLab/IPSAE2, pinned commit
`fe560f4c2461eeb91aebded61ac8c1e83b264f60`, committed
2026-11-13T23:38:21-05:00; package version 0.1.11, alpha research release.
Reference checkout: `/tmp/nyl12_ipsae2_review_20261127`.
Read-only source review, a raw-CIF inventory and a historical native PDB/PAE
reader check; no dependency
installation, scoring benchmark, prediction launch, or campaign edit.

## Recommendation

Use as an additional assembly/pose-confidence descriptor after validating native
mapping on matched WT/design outputs. Do not treat scores as measured binding,
catalytic competence, or broadly calibrated ligand success probabilities.
The README establishes six-decimal agreement with an older implementation for
six contacting AURKA–TPX2–ATP–Mg pairs in each AF3/Boltz example, plus batch
execution on 1,600 seeded AURKA models per predictor. This establishes numerical
regression and throughput, not calibration on Nyl12 or PA6 oligomers.

### Pipeline decision and bounded qualification

**Recommend adding ipSAE2 as a supplementary confidence descriptor after
technical qualification; do not add it as a catalytic pass/fail gate.** It
can address assembly and substrate-placement confidence while our geometry
and physical validation address different questions. Keep the current matched
Boltz/Amber and QM/MM/US jobs unchanged.

The proposed insertion point is immediately after native Boltz prediction,
before any canonicalization, trimming or Amber relaxation. Retain the raw
coordinate/PAE/confidence triplet and native token/atom order. Join scores to
later canonical site identities through the explicit mapping; a relaxed model
does not acquire a new predictor confidence from the original PAE.

For Nyl12, retain protein-chain interface scores, biologically paired alpha/beta
group scores and protein-aligned PA6 scores for all four ligands. The group
membership must use verified native label IDs. Report both directions and
explicitly select the protein-aligned direction for ligand comparisons;
the grouped scalar maximum can select the reverse direction. Preserve missing
contacts and undefined scores explicitly, rather than replacing them by zero.
Summarize all four sites and seed variation without treating copies as
independent replicates. A whole-complex maximum is insufficient.

For this study, a further limitation follows from scoring the entire flexible
PA6 L4 ligand: poorly localized distal segments could lower confidence even
when the scissile amide is plausibly placed. Conversely, a confident pose can
have the wrong catalytic orientation. This is a study-specific inference,
not a demonstrated ipSAE2 result for PA6. Retain whole-ligand scores alongside
our amide-specific geometry and do not crop/reorder PAE inputs to improve a
score or label a custom local diagnostic as the published ipSAE metric.

Next qualification: a pinned, separate uv environment; upstream examples and
tests; then one raw WT and one design, with native-order verification and
direction/group checks, before processing all six matched models. Confirm
actual PAE availability rather than inferring it from CLI flags. Existing
historical WT matrices are available, so no prediction rerun is justified
merely to obtain a test input. Do not upgrade the active Boltz environment.
No new production dependency has been added; installation remains subject to
the project's confirmation rule.

The supplied primary preprint is Dunbrack (2026),
[DOI 10.64898/2026.11.02.763449](https://doi.org/10.64898/2026.11.02.763449).
The publisher full text was inaccessible during this review; bibliography and
the author abstract were cross-checked against public indexing, while the
implementation was inspected directly. The abstract reports a held-out
Boltz-2 threshold of 0.8 for a 5 Å ligand-RMSD classification task. Such a
criterion does not establish catalytic geometry, correct cleavage register,
PA6 affinity, barrier reduction or improved hydrolysis, and will not become a
Nyl12 hard filter. Full benchmark details and transferability need review
before any quantitative calibration claim.

## Required native artifacts and exact mapping

Supply each unmodified native prediction CIF, its corresponding native
`pae_*_model_N.npz`, and matching `confidence_*_model_N.json` explicitly.
Confidence is optional for calculating ipSAE but required to retain native
ipTM comparisons. NPZ loader accepts `pae`, `token_pair_pae`, or
`predicted_aligned_error`. Native confidence alone cannot reconstruct PAE.
Do not score canonicalized/reordered/truncated PDB coordinates against native
PAE without an independently verified mapping.

`readers.py:101–128` uses `_atom_site.label_asym_id`; author IDs are only a
fallback when label IDs are absent. Heavy atoms remain in coordinate encounter
order; H/D and water are excluded. Canonical amino acids use one token per
residue; ordinary ligands can use one token per heavy atom. `infer_tokens`
checks matrix dimension and any supported source chain metadata, rejects
multiple distinct compatible mappings, and checks compatible atom-to-token
metadata. Boltz writer NPZ contains only `pae`, so count compatibility alone
must not be promoted to proof of source order. Retain native processed token
and structure artifacts and explicitly verify chain blocks and atom order.
Absent frame metadata, ligand reference frames are inferred geometrically;
record this policy, especially when interpreting ligand-aligned directions.

Inspected historical raw reference:
`design/nyl12_relay_ts_ensemble_20261126/broader_pilot/wt_paired_precursor_controls_20261127/canonical_precision_attempt01_archive/native_models/WT_paired_seed_20261127/original_native_prediction.cif`.

| Native label IDs | Type | Residues per chain | Heavy atoms per chain | Expected PAE tokens per chain |
|---|---|---:|---:|---:|
| A–D | protein | 232 | 1,736 | 232 |
| E–H | protein | 96 | 744 | 96 |
| I–L | ligand LIG1 | 1 | 33 | 33 |

For this reference, expected native matrix is **1444 × 1444**:
4 × (232 + 96) + 4 × 33. These counts are not a validation of every design.
Queued matched-campaign eligibility remains pending actual matching native
outputs, dimensions and source token evidence; this audit found no prediction
triplets in the locally enumerated original_backbone_matched_20261127 tree.
Other hosts were not inspected.

## Groups and interpretation

After checking chain identity and precursor pairing, use disjoint virtual
protein groups `--group A+E --group B+F --group C+G --group D+H` to assess
mature protomer assembly and each protomer versus individual I/J/K/L ligands.
Retain original chain-level rows. Assess whole protein versus each ligand in a
separate grouping run with `--group A+B+C+D+E+F+G+H`; whole protein versus
`--group I+J+K+L` is a further pooled descriptor, capable of hiding bad sites.
Avoid comparing its scalar directly to ordinary protein–ligand rows because
grouped scalar selection differs. Preserve all four ligand-site results and
both directional columns rather than ranking on one pooled maximum.

`ipsae12` aligns on chn1 and scores chn2; `ipsae21` reverses this.
Ordinary same-type scalar is the larger valid direction; ordinary mixed-type
scalar aligns on higher-ranked protein > nucleic acid > ligand > ion.
Grouped scalar always takes the larger valid direction, even for mixed types.
The algorithm maximizes over eligible alignment frames, then averages the TM
kernel over selected scored units. Protein/nucleic-acid units default to
PAE ≤10 Å; ligand-only scored entities include all heavy-atom tokens, without
a PAE cutoff, with d0=2 Å. Protein d0 is max(1 Å,
1.24*cuberoot(n_protein−15)−1.8), where n_protein counts qualifying protein
residues in that frame, not full chain length. Mixed protein groups apply the
PAE cutoff to every scored group unit but only protein residues set d0.
Contact ≤5 Å filters which pairs are reported, not which selected units enter
the score. Use `--all-pairs` to retain noncontacting assembly comparisons.
Report n0, d0, alignment tokens, contacts and both directions with the scores.

## Installed Boltz flag: important local implementation caveat

Installed source root:
`/data/jpc/enzyme-substrate-modeling/.pixi/envs/boltz/lib/python3.10/site-packages/boltz`.
`main.py:1081` declares **`--write_full_pae`**; Click boolean flag and Python
signature default are false (:1306), despite help saying true. It is placed
in `predict_args` (:1624). However, the installed Boltz2 `predict_step`
(`model/models/boltz2.py:1100–1101`) returns PAE whenever confidence prediction
and alpha_pae > 0; it does not read write_full_pae. Writer
`data/write/writer.py:232–238` saves PAE whenever the prediction contains it.
Thus absence of this flag does **not** establish missing PAE for this installed
Boltz2. Check actual artifacts before proposing reruns or modifying a frozen
campaign. Setting it is an explicit output request; no dependence on it was
found in local Boltz2 sampling/scoring code. This is a source audit, not a
numerical replay proving identical trajectories under changed invocations.

## Requirements, provenance and remaining gates

Python ≥3.10; MIT license. Declared production requirements: NumPy ≥1.24,
pandas ≥2.0, SciPy ≥1.10, Gemmi ≥0.6.7, Matplotlib ≥3.7. No GPU, PyTorch,
Boltz or AF3 installation required for the scorer. Existing Boltz environment
has NumPy/pandas/SciPy/Gemmi but no Matplotlib; no package installed here.
A separate uv-managed environment and PI confirmation precede adding a new
production dependency. No wall-time/memory benchmark was run.

`CITATION.cff` is stale (0.1.5, 2026-10-28) and says a separate manuscript is
in preparation; it does not verify the supplied November preprint citation.
Primary paper text was not reviewed here. Before production, pin revision and
lock environment; validate one WT plus one eligible design end-to-end with
native token tables, expected dimensions, group membership and full provenance;
then assess matched seeds without introducing universal pass thresholds.

Source SHA256:
- readers.py: ecc400becfd10f350185c773277925d34808791dfab98c1c8bf84d1f0ad7456c
- scores.py: 0de5a233ed888e2135c0ad20a1c36104e8a6e11c18f1224faa37c7a4b8b323dd
- installed Boltz main.py: da55e701d35917bf569d28b7ad064f1da42e7acb20a1b4b20473a57821ae45c7

Verification: source inspection and read-only raw CIF parsing with the pinned
reader; dependency-presence check. No code changes, automated tests or coverage
measurement; no scientific acceptance claims.

## Additional historical WT check

Matching historical native PDB/PAE/confidence directory (no CIF present):
`design/nyl12_ligandmpnn_pilot/boltz_validation/predictions/WT/seed_01/boltz_results_WT/predictions/WT`.

- Native WT PDB/PAE shape: (1444, 1444); inferred policy modified=atom;ligand=atom.
- Chain token counts: {'A': 232, 'B': 232, 'C': 232, 'D': 232, 'E': 96, 'F': 96, 'G': 96, 'H': 96, 'I': 33, 'J': 33, 'K': 33, 'L': 33}
- numpy: 1.26.4
- pandas: 2.3.3
- scipy: 1.13.1
- gemmi: 0.6.5
- matplotlib: absent

Count agreement is not independent proof of identical-copy order. Exchanging identical protein copies or identical ligands can preserve counts and chain block lengths. Verify source token asymmetry IDs/native serialization before translating canonical identities. For grouped protein-versus-ligand reporting use the explicit protein-aligned direction; grouped scalar max may choose the reverse ligand-aligned direction. This historical PDB check does not certify a current raw-CIF mapping.
