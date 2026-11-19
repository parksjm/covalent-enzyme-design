# Alternative proton routes for mature Nyl12 PA6 hydrolysis — 2026-11-10

Documentation assessment only: no calculation, protonation change, job control,
or NEB input change. Rankings below are priorities for discrimination using
existing evidence, not established mechanistic probabilities or barrier rankings.

## Ranking and chemical interpretation

1. **Direct OG1–HG1 donation to the leaving amine during acylation** is the
   best-supported local donor geometry. OG1 is near leaving N, whereas terminal
   N is farther away. The running unconverged band already samples an OG1 donor
   route after a terminal-N proton detour; its initialization does not establish
   that the detour is necessary. A compact four-center O–H–N–C transition geometry
   remains a hypothesis, with potentially unfavorable strain/proton orientation.
   The retained restrained TI1-like seed has HG1–leaving-N **4.093 Å** despite
   OG1–leaving-N **2.448 Å**: proximity of heavy atoms alone is insufficient.
2. **Asp272/Asp274-assisted activation or electrostatic stabilization coupled
   to OG1 donation** merits consideration before assuming direct Asp acid donation
   to leaving N. Both carboxylates contact terminal N much more closely than
   leaving N. This supports a possible role in N-terminal protonation/activation
   or polarization, but does not demonstrate an operative relay. A carboxylate
   cannot donate an acid proton it does not possess; specify its uptake source,
   proton inventory and geometry before proposing Asp-mediated donation. Direct
   Asp-to-leaving-N acid transfer has poor support in these snapshots.
3. **Water-assisted acylation, including an Asp/water network**, remains an open
   alternative with less structural support in the inspected seeds. The earlier
   attack09 audit found no water O within **5 Å** of leaving N. It found water
   O49306 (Amber residue11239) **3.141 Å** from Thr N. This does not establish a
   one-water bridge and does not exclude water recruitment, another conformation,
   multiple waters, or a different substrate pose. One snapshot cannot prove an
   absent water relay. Existing MM waters cannot exchange protons in this model.

**NH3+ clarification:** the proposed protonated N-terminal microstate still uses
**OG1–HG1 as the leaving-amine proton donor**; NH3+ need not donate that proton.
This is chemically distinct from an NH3+-to-leaving-N hypothesis. A concerted
OG1 attack/proton donation can be investigated without requiring NH3+ to accept
HG1. Whether that route is favorable, or an Asp/water-assisted activation step
is needed, remains untested. Keep a complete proton ledger and define matched
ES/acyl endpoints with unchanged atom inventory and charge within each microstate.
Compare barriers relative to each microstate's own ES; include protonation
populations when comparing rates. The neutral-NH2 calculation does not select
the preferred N-terminal protonation state.

## Retained distance evidence

Distances in Å from existing neighbor-audit JSON; nearest Asp carboxylate oxygen
is OD2 for every listed pair. Biological Thr233/Asp272/Asp274 correspond to
Amber residues929/968/970. ES is the frame1603 optimized structure; TI1-like is
the restrained corrected attack09 seed, not an accepted stable intermediate;
acyl is the reaction-unrestrained converged endpoint, without a Hessian test.

| Pair | ES | TI1-like | Acyl |
|---|---:|---:|---:|
| Thr233 OG1–leaving N | 3.040 | 2.448 | 3.970 |
| Thr233 terminal N–leaving N | 4.444 | 4.937 | 5.614 |
| Thr233 terminal N–OG1 | 2.692 | 2.612 | 2.656 |
| Asp272 OD2–Thr233 terminal N | 3.161 | 2.605 | 3.225 |
| Asp274 OD2–Thr233 terminal N | 3.698 | 2.743 | 3.021 |
| Asp272 OD2–OG1 | 5.443 | 4.323 | 5.385 |
| Asp274 OD2–OG1 | 4.920 | 4.651 | 5.245 |
| Asp272 OD2–leaving N | 6.242 | 5.988 | 7.207 |
| Asp274 OD2–leaving N | 7.125 | 6.915 | 7.022 |
| Reported nearby water11239 O–Thr233 terminal N | 3.433 | 3.141 | 3.199 |

The water row reports the retained nearby-water record, not a complete solvent
network analysis. JSON distances are Cartesian, without minimum-image wrapping;
no new distances, angles, occupancy statistics, free energies or pKa values were
computed here. Donor–H–acceptor angles and physically populated water networks
are needed before asserting a transfer-ready hydrogen-bond chain.

## Hydrolysis, maturation and matched comparisons

These proposals concern **mature-enzyme substrate hydrolysis**. Acylation yields
an ester acylenzyme plus amine; deacylation consumes water, passes through a
proposed TI2 and regenerates Thr–OH with the carboxyl product. A water reactant
for deacylation is distinct from a water relay during acylation. Neither is
established by merely obtaining an acylenzyme endpoint.

The retained NylC2023 discussion proposes an Asp308–Asp306–Thr267 network for
hydrolysis as well as investigating autocleavage. Nyl12 numbering correspondence
and Asp conservation support testing a related network, but maturation defects
alone cannot identify catalytic acid/base function in the mature enzyme.
Conservation also reflects folding and maturation constraints.

The separate Lys156/Gly273 QM expansion is a controlled sensitivity design,
not evidence for an additional proton shuttle: Lys NZ–Asn OD1 is3.237/3.120/3.143Å
and Gly273 backbone N–OG1 is2.985/2.926/3.275Å (ES/TI1-like/acyl). Its proposed
133-real-atom region includes both existing Asp sidechains and intact nearby
backbone; changing region does not establish protonation or resolve all boundary
errors. Keep the running NEB unchanged. Any future water test needs an identified,
physically supported water included consistently in QM at both endpoints and all
images. Treat substrate orientation, N-terminal microstate and QM partition as
separate controlled comparisons before combining them.

## Sources and scope limits

- [Backward acylation record](acylenzyme_backward_path_20261110.md): HG1 distance,
  water audit, endpoint identity, band limitations and retained primary-literature
  interpretations (NylC DOI10.1111/febs.16755; four-center warning
  DOI10.1021/acs.jcim.1c00666).
- [Existing conservation and structural assessment](nyl12_active_site_conservation_20261110.md).
- Existing geometry and source SHA256 records:
  [ES](../analysis/nyl12_active_site_conservation_20261110/es_qm_neighbors.json),
  [TI1-like](../analysis/nyl12_active_site_conservation_20261110/ti1_qm_neighbors.json),
  [acyl](../analysis/nyl12_active_site_conservation_20261110/acyl_qm_neighbors.json).
- [QM expansion design](qm_region_expansion_design_20261110.md).

No new literature search or scientific calculation was performed. The DFTB3 band
is unconverged; a band maximum is not an activation barrier, and no inspected
snapshot proves route viability or exclusion. The four-center protease precedent
was high energy, and the PET ester analogy does not establish PA6 amide chemistry.
