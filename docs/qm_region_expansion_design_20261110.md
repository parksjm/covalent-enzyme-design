# Controlled QM expansion design — 2026-11-10

Design only. No calculation input, protonation state, coordinates or running job was changed.

Recommended one-based Amber mask:

```text
qmmask='@1621-1655,2268-2283,2727-2734,13905-13931,14549-14579,19914-19917,19925-19931,19937-19941'
qmcharge=-1
```

The proposed region has **133 real atoms + 9 hydrogen links = 142 centers**, versus
98 + 8 currently. Its capped electron count is 498, consistent with a proposed
closed-shell singlet. Formal charge changes from −2 to −1 because the added
Lys156 sidechain remains NZH3+; no proton is added, removed or transferred.
The selected MM partial-charge sum (−0.814849) is not the formal QM charge.

| Component | Exact Amber atom indices | Reason |
|---|---|---|
| Lys156 sidechain | 2268–2283 | CB through NZ and HZ1/HZ2/HZ3; leaves N/CA/carbonyl in MM |
| Leu271 carbonyl cap through Asp274 before its carbonyl | 14549–14579 | Adds intact peptide backbone around Gly273 while retaining both existing Asp sidechains |
| Gly273 peptide donor | N14563–H14564 | Biological Gly273 is Amber residue 969; Amber residue 273 is unrelated GLN |

The contiguous block is Leu271 C/O + all Asp272 + all Gly273 + Asp274 excluding
its final C/O. Biological Asp272/Gly273/Asp274 map to Amber residues 968/969/970.
Including Leu271 C/O is deliberate: it avoids cutting its peptide C–N bond into
Asp272. The backbone expansion adds 19 atoms beyond the 12 Asp-sidechain atoms
already selected; Lys adds 16, for 35 added real QM atoms total.

New boundary cuts are exclusively carbon–carbon:

| QM atom | MM atom |
|---|---|
| Lys156 CB2268 | Lys156 CA2266 |
| Leu271 C14549 | Leu271 CA14534 |
| Asp274 CA14572 | Asp274 C14580 |

The former Asp272 CA14553–CB14555 and Asp274 CA14572–CB14574 cuts disappear.
No selected boundary cuts a peptide carbonyl C–N bond. The **inherited
ILE112 N1619–CA1621 polar boundary remains**; it is an N–CA cut, and this controlled
expansion does not repair or validate it. Exact all-nine boundary records and
all 133 mapped atoms are in [the machine-readable proposal](../analysis/nyl12_active_site_conservation_20261110/proposed_qm_expansion.json).

All proposed QM atoms lie inside the existing fixed-reference 12 Å movable shell.
Keep that shell, the physical states, DFTB3/3OB parameter set and convergence
settings identical. Retain qmmm_int=1, adjust_q=2, qm_ewald=1, qm_pme=1,
qmcut=10 Å, MM cut=10 Å, and the current hydrogen-link treatment. No explicit
`polarcut` keyword was found in the inspected base input or Sander source
subtree; the inherited polar boundary and embedding/cutoff sensitivity still
need separate qualification. Enlarging QM is a sensitivity test, not proof that
polarization or electrostatic truncation errors are resolved.

Before interpreting results, verify the installed engine reports 133 real QM
atoms, nine links and charge −1, and audit unchanged proton inventory. Use valid
reconstructed XYZ-based coordinates while the native numbered DL-FIND NetCDF
writer offset defect remains unresolved. Compare reaction energy differences
within each QM partition and their geometries; raw absolute total energies
across different partitions are not directly comparable. No calculation was
launched to validate this design.

Structural verification: exact topology/mask selection, all crossing bonds,
closed-shell electron parity, Lys NZ three-proton inventory, intended residue
identities, and unchanged active-shell inclusion passed. Source topology,
reference and config hashes are retained in the JSON proposal.
