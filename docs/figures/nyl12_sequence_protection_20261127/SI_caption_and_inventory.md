# SI figure: protected and variable Nyl12 sequence positions

**Figure Sx. Sequence identities protected during symmetry-tied Nyl12 design.**
(A) Complete modeled WT TS-like tetramer, with 95 other protected protein positions magenta,
six mechanistic core positions teal (F113, K156, N186, T233, D272 and D274),
and 227 variable positions very light gray and PA6 carbon orange (N blue; O red). Protection
status deliberately replaces the usual chain colors. This is a modeled TS-like
parent, not an experimental structure or a validated transition state. All 328
precursor positions are modeled in each of four copies. (B) Canonical precursor
mask, shown as four successive 82-residue rows. The restored K156 policy fixes
101/328 positions and permits 227/328 to vary: 404 protected and 908 variable
physical residues in the tetramer (1312 total). Sequence identity protection does
not constrain backbone mobility.

Magenta plus teal marks the complete 101-position protected union. The six
teal positions are a mechanistic subset, including homolog-informed K156;
this label does not claim direct experimental proof of catalysis in Nyl12.

The union contains 24 core policy positions and 77 additional recurring interface
positions. The core includes 17 PI hard protections (F113, N186, Q221–N232,
T233, D272 and D274), plus seven precautionary/context positions (M1, K156,
T234, L235, G273, V275 and L276). These categories state design policy, not a
claim that every residue is individually experimentally essential. K156 retention
is a precaution informed by homologous NylC Lys189 evidence; L235 retains mature
beta N-terminal context without demonstrated residue-specific catalytic function.
The full 221–232 alpha tail is retained by PI instruction.

Interface protection requires an inter-repeat heavy-atom contact at ≤4.0 Å in
at least three of four equivalent copies **in each of all three original TS-like
parents**. Intrinsic alpha/beta pairs A/E, B/F, C/G and D/H are excluded from
this contact census. These recurring contacts are a structural design heuristic,
not sequence-conservation evidence or proof of biological necessity. There are
no detected disulfide protections. Overlap is retained in the machine-readable
per-position reasons, so the union is not double-counted.

Source: `design/nyl12_relay_ts_ensemble_20261126/templates/medoid_01_delivery_full_010_frame_068/complex_heavy.pdb`.
Census uses this parent and original medoids 02 and 03, with exact paths and
SHA-256 values in `canonical_protection_restored_K156_20261127.json`. Rendering
uses a proper rigid protein-landmark camera (Phe113 upper left, Asn186 upper
right, Thr233 below), white background, valence 0 and 0.2 Å sticks. No coordinate
optimization is performed. `manifest.json`, `protected_positions.pml`, the PyMOL
session and `build_figure.py` preserve reproducibility. The crystal structure's
missing tail is not silently substituted for this complete modeled assembly.

## Exact protected WT identities (canonical precursor numbering)

M1, Q41, I46, Q47, G48, G49, S50, V51, G52, T53, F54, F55, T56, Q57, E58, K59, Q61, Q62, E64, Y66, I77, L78, L80, E81, S88, F91, K96, N97, E98, V99, Q100, F101, R102, M104, F113, Y115, M121, Y123, K126, K156, L157, N186, K195, K208, N212, V215, L218, L219, Q221, P222, K223, K224, P225, F226, W227, P228, E229, P230, K231, N232, T233, T234, L235, P246, R247, K250, Q251, R254, Q255, H257, H258, S261, Q262, H265, A268, T269, I270, L271, D272, G273, D274, V275, L276, R282, D287, G292, A293, D294, I295, E296, D298, L299, N300, A301, K302, I304, Y305, L306, M308, V309, E312

## Core policy positions (24)

M1, F113, K156, N186, Q221, P222, K223, K224, P225, F226, W227, P228, E229, P230, K231, N232, T233, T234, L235, D272, G273, D274, V275, L276

## Complete mutability inventory and candidate audit

`canonical_position_inventory.json` is authoritative; the CSV provides the same
328 rows for SI tables, including WT identity, chain-local mapping, equivalent
copies and all applicable protection reasons. `candidate_eligibility_audit.json`
records the read-only snapshot of all 24 preserved `relay_corr01_*` sequences.
All 24 mutate K156 (I, L or T), with no other violations of the 101-position mask;
therefore none is eligible under the restored policy. Existing generation and
scoring artifacts are preserved. Twenty-one completed candidate score receipts
were present at the audit snapshot; this does not certify scientific acceptance.
