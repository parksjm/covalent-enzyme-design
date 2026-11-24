# WT PA6 pose comparison

Run `pymol -r /absolute/path/to/view_poses.py`. Use PyMOL scenes `flipped`, `positive`, and `overlay`; the loader locates PDBs relative to itself. `comparison.pse` and rendered PNGs are provided. Protein carbons are cyan, PA6 carbons orange, nitrogen blue and oxygen red. Small spheres identify Thr233 Nα and PA6 N59; black/yellow dashed lines measure their separation. Beta chain local residue 1 is precursor Thr233; 40/42 are Asp272/274.

| Source | Site | Signed projected φ | Nα–leaving N |
|---|---|---:|---:|
| WT seed 2026112204 | AE, flipped | −105.481° | 4.58775 Å |
| WT seed 2026112202 | BF, positive | +96.435° | 3.61882 Å |

These are individual examples, not group averages. Relaxed WT negative φ: 4.65 Å mean, range 4.24–5.47 Å (16 sites); positive φ: 3.61 Å mean, range 3.51–3.67 Å (4 sites). Both depicted sites fail the held-out native-Thr gate. Sign describes orientation under consistent scissile-amide mapping; distance alone does not establish relay competence or proton-transfer direction. Negative poses may require reorganization or water-mediated relay. The flipped-only MD choice tests this hypothesis and is subject to starting-state selection bias.

The positive structure was rigidly aligned to the flipped structure using 16 N/CA/C/O backbone atoms from canonical Phe113, Asn186, Asp272 and Asp274. Thr233 and PA6 were excluded from fitting. The extracts retain source heavy atoms; the positive pose underwent only a rigid transform and PDB rounding to 0.001 Å, with no added hydrogen or optimization. Alignment RMSD, original and extracted file SHA-256 hashes, exact source paths and full precision measurements are in `manifest.json`. `build.py` reproduces the extraction inside the source checkout using the existing alignment and geometry functions.

Both aligned poses use the same anchor-defined camera: Thr233 Nα below Phe113 and Asn186; bond valence display is disabled. Camera changes do not alter atom coordinates.

The retained calibration and ensemble comparison tables disagree on the flipped example’s tier label; tier classification is omitted from this figure. Their distance and signed φ agree.
