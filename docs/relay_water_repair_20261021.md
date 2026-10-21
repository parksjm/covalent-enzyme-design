# Scripted TI1 relay-water repair — 2026-10-21

## Chemistry and scope

The substrate names are PA6 L4 (`pa6_l4`) and MAXA|MAXA (`maxa_maxa`).
Historical `pa6_l2` and `poly_xama_maxa` input identifiers are provenance
aliases, not the current substrate names. L2 and MAXA remain the respective
hydrolysis products; renaming does not change their molecular graphs.

For these TI1 seeds, catalytic Thr233 is PDB E1 / Amber residue 929. Its
N-terminal group is NH3+, OG1 has no proton and is covalently bonded to the
tetrahedral substrate carbon; the substrate oxygen is an oxyanion. The relay
must connect **N–H…Ow–H…OG1**. A pair of O…N/O distances alone is insufficient.
This protocol places neutral H2O and does not move a proton between atoms.

The original water positions fail the retained steric gate (maximum contact
shortfall 0.6866 Å PA6 / 0.6613 Å MAXA|MAXA; allowed 0.10 Å). The prior
5832-point fixed-pocket oxygen search found zero eligible sites. Excluding
only Asn186 side-chain atoms and Thr929 N hydrogens exposed 9 / 5 oxygen
sites. This motivates bounded local rotations, not manual coordinate editing.

## Reproducible search and acceptance

The versioned entry points are `scripts/launch_nyl12_relay_repair.py` and
`scripts/repair_nyl12_relay_water.py`. The launcher retains input copies and
source snapshots, hashes, software versions, explicit arguments/environment,
and a runnable tmux job script. Results use new attempt directories. Original
models and previous failed attempts remain intact. The live action journal
records preparation, launch, inspection, exports, and repository operations.

1. Validate the original TI1 atom inventory and connectivity. Search the same
   water-O sphere-intersection grid: distances to N and OG1 each 2.6–3.4 Å in
   0.1 Å increments, with 72 azimuths. Keep every site passing the optimistic
   fixed-environment screen, not merely the best-ranked site.
2. Rotate Asn186 chi1/chi2 and Thr NH3 hydrogens in explicit 15° increments.
   These are delta rotations from the input structure. Attached hydrogens move
   with the Asn side chain; all backbone, substrate and Thr heavy atoms stay
   fixed. Require Asn ND2–oxyanion O 2.4–3.5 Å.
3. Check moved atoms against the full environment with only topology 1–2 and
   1–3 exclusions. Geometric contact limits (Å): H–H 1.40, H–C 2.10, H–N/O
   1.55 for polar H or 2.00 for nonpolar H, H–S 2.30, C–C 2.80, C–N/O 2.75,
   C–S 3.00, N/O–N/O 2.50, N/O–S 2.95, S–S 3.20. Maximum allowed shortfall
   is 0.10 Å. These are screening criteria, not a force field or energy model.
4. Require N–H…Ow angle ≥140° and H…Ow 1.4–2.5 Å. Construct rigid H2O
   (O–H 0.9572 Å, H–O–H 104.52°), aim H1 at OG1, and sample the H2 roll
   every 10°. Do not substitute an oxyanion-directed water for the OG1 bridge.
   Require all water atoms to pass the original full-environment steric gate.
5. Rank passing candidates geometrically by small Asn/heavy and NH3/H
   displacements, with a modest preference for 2.85 Å bridge distances.
   Retain sampling settings, stage failures, original and selected geometry,
   and the exact selected rotations. A finite-grid failure is not proof that
   no valid geometry exists in continuous coordinate space.
6. Independently verify fixed coordinates and non-water bond-length invariance,
   write an Amber restart, and rerun the unchanged TI1 preflight on that file.
   Export a PDB only after acceptance. Record d1, d2, d3, a_BD and
   herschlag_phi before and after. Process exit status and scientific acceptance
   are separate records; completion can legitimately yield no accepted seed.

Accepted outputs, if found, are **geometric starting seeds**, not optimized
TI1 minima, TS structures, parameter validation or evidence of catalytic
activity. Fixed pre-existing distortions are not repaired by this local search.
No failed candidate is silently promoted. A remaining backbone obstruction
requires a separately documented bounded-repair protocol before refinement.

## Inspection export

`scripts/prepare_nyl12_inspection_download.py` generates canonical model names
and copies the versioned `scripts/nyl12_inspection_view.pml` into the bundle as
`protein/view_models.pml`. The viewer includes `set valence, 0`,
`show sticks, resi 113+186 and chain A`, and `util.cnc`. All exported files have
source paths and hashes. No coordinates are adjusted interactively in PyMOL.

## AF3 alternative: executed feature preflight, prediction pending

The user proposed a minimal alpha–beta/covalent-substrate/water test. Scripts
`launch_nyl12_af3_water_test.py` and `check_af3_water_features.py` launched the
PA6 L4 test in tmux with source/input snapshots. Explicit HOH is retained as a
non-polymer ligand despite the pipeline water-chain filter. Actual featurization
retained both N/OG1-to-water restraints. Inspect the attempt job status/log and
prediction outputs before claiming inference success. AF3 lacks explicit water
H orientations and does not validate terminal protonation.
