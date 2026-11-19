# Negoro 2023 NylC mechanism audit

Audited 2026-11-12. Source: Negoro et al., *The FEBS Journal* **290**,
3400–3421, [doi:10.1111/febs.16755](https://doi.org/10.1111/febs.16755).
The full main-text extraction (1,301 lines, including methods and references)
was read in bounded chunks from
`analysis/negoro_2023_protonation/Negoro_2023_FEBS.txt`.
Supplementary information and rendered figure artwork were not independently
audited here; captions and references to SI below come from the main text.

## What the paper establishes about protonation

The paper **proposes** protonated Asp306 as a general acid in precursor
autocleavage; it does not directly determine its protonation or a microscopic
pKa. On p. 3414, step 2, the authors argue that Asp306-COOH donates a proton
to the scissile amide nitrogen during conversion of the tetrahedral intermediate
to an ester. Their neutral-pH protonation assertion is part of this mechanistic
proposal. The argument combines proximity, complete loss of maturation in D306A,
and hydroxylamine cleavage of the Cys267 precursor, consistent with a thioester.
These observations do not uniquely distinguish protonation microstates.

The corresponding proposed maturation step 3 uses **Asp306-COO−** as a base to
activate solvent water for ester hydrolysis (p. 3414). Therefore, the paper's
mechanism changes Asp306's protonation across the reaction. It does not support
holding both Asp residues neutral, nor does it directly establish which oxygen
of either Asp is protonated in the mature substrate-bound enzyme.

The 500 ns classical MD concerns the **uncleaved precursor**, reconstructed from
the Y130/A266/A267 crystal structure. The final 200 ns provided water/geometry
analysis (pp. 3412–3413, 3418). The main methods identify AMBER20, ff19SB and
TIP3P but do not specify explicit Asp306/Asp308 or Lys189 protonation choices,
pKa calculations, constant-pH sampling, or reactive proton-transfer calculations.
Classical hydrogen-bond observations are evidence for accessibility, not a
demonstration of the proposed proton transfer or its barrier.

## Maturation and substrate hydrolysis are distinct mechanisms

For **maturation**, the authors favor Lys189 as the base that activates Thr267-OH
through water-mediated or water-inaccessible networks involving Asn219 and
Gly265. No Thr267-Oγ/Asp306 or Asp308 interaction was observed during their
500 ns precursor MD. They nevertheless explicitly retain the possibility that
Asp306/Asp308 participate in Thr activation during slower motion beyond their
simulation (p. 3414). Their absence-of-contact result is not a definitive
exclusion of an Asp-mediated route.

For **mature substrate hydrolysis**, the proposed mechanism assumes transfer of
BapA chemistry to NylC: Asp308 activates Thr267-Oγ through Asp306 and the newly
available Thr267 α-amino group (p. 3415, Fig. 10). This is a homology-based
proposal, with a different N-terminal inventory from the uncleaved precursor.
It is not the Lys-only maturation mechanism repeated without qualification.

D306A/D308A completely block precursor processing. The resulting loss of
hydrolysis activity cannot by itself isolate their contribution to the mature
catalytic reaction: mature enzyme was not produced. Conversely, Lys/Tyr mutant
hydrolysis activities were corrected for their partial cleavage, making those
results more informative about mature substrate turnover (Table 1).

The assays use Ahx-cyclic oligomers, with TNBS monitoring newly formed amino
groups at 30 °C (p. 3417). TLC detects Ahx linear dimer/monomer products. These
data do not establish the project's specific PA6 L4→L2+L2 reaction register,
poly(X,AMA) product identity, or a quantitative barrier for either substrate.

## Lys189, Tyr146 and Asn219

| NylC site | Evidence and interpretation |
| --- | --- |
| Lys189 | K189A leaves approximately one-third eventual maturation but only 0.6% cleavage-corrected H-GYAQ hydrolysis activity (0.032 U/mg). Strong functional evidence for mature catalysis; its exact role is not resolved. It is not conserved in BapA. |
| Tyr146 backbone NH | Main-chain NH binds precursor water Wat94 (2.86 Å); together with Asn219-NH2 it is proposed to stabilize substrate oxyanion in mature hydrolysis (p. 3415). This donor can remain after Tyr→Phe substitution. |
| Tyr146 side-chain OH | Y146A retains up to 92% maturation but 3.6% cleavage-corrected H-GYAQ hydrolysis activity. NylCA Y146A kcat/Km falls to about 3%; both kcat and apparent affinity change (Table 2). Side-chain involvement is supported, but proton donation is suggested by analogy with NylB, without a substrate-bound NylC reaction structure or direct Tyr proton-transfer evidence. |
| Asn219 | N219A maturation plateaus near 31–33%; water-network involvement is supported. Its side-chain NH2 is proposed as an oxyanion donor during mature hydrolysis; the main text does not separately quantify mature N219A hydrolysis kinetics. |

Y146A removes an aromatic ring as well as the OH. It cannot establish that a
natural Phe homolog lacks activity, or that its side chain can donate a proton.
The Lys189-NH2 base formulation likewise requires an explicitly chosen neutral
Lys microstate in any test; the free-Lys pKa quoted in the paper does not
determine the active-site pKa.

## Sequence mapping to Nyl12

Mapping was recomputed by independently counting nongap residues in both rows
of the existing alignment, rather than assuming an offset:
`analysis/nyl12_active_site_conservation_20261110/nyl12_nylc_numbering_alignment.fasta`.
Its Nyl12 alpha query is also present in `conservation.json`; this confirms F113
and K156 directly in the actual modeled sequence. The NylC row is 3AXG/NylCA
(355 residues), and Nyl12 is the modeled concatenated alpha+beta sequence
(232+96 residues). Thus these are homolog alignment assignments, not proof
that their active-site geometries are identical.

| NylC | Nyl12 | Consequence |
| --- | --- | --- |
| Phe134 | Phe101 | Conserved aromatic residue |
| Tyr146 | **Phe113** | Backbone NH remains; Tyr phenolic OH is absent |
| Lys189 | Lys156 | Conserved Lys identity; protonation remains a hypothesis |
| Asn219 | Asn186 | Conserved Asn identity |
| Asn266 | Asn232 | End of modeled alpha chain |
| Thr267 | Thr233 | Start of modeled beta chain |
| Asp306 | Asp272 | Asp retained |
| Gly307 | Gly273 | Gly retained |
| Asp308 | Asp274 | Asp retained |

The offset changes from 33 to 34 near the cleavage site. In particular, replacing
NylC Tyr146 with a hypothetical Nyl12 Tyr would silently introduce chemistry
that the actual sequence lacks. F113 side-chain packing and its backbone NH
should be evaluated separately from any Tyr phenolic mechanism.

## Mature Fig. 10 crystal: all Asp oxygen distances

Fig. 10 explicitly identifies mature NylCp2-G122Y130A36Q263 **5Y0M**, not the
precursor structures newly reported in 2023. The deposited PDB was retrieved
from [RCSB](https://files.rcsb.org/download/5Y0M.pdb), preserved without edits as
`analysis/negoro_2023_protonation/5Y0M_deposited.pdb`.
It reports 1.03 Å resolution and two independent deposited chains A/B; the
biological tetramer uses twofold-related copies. No OD alternative locations
occur at Asp306 or Asp308 in these chains.

Distances below are unrelaxed Euclidean measurements in Å. The minimum over
all four OD1/OD2 pairs is invariant to oxygen-name swapping. Applying the
deposited rigid symmetry operation to an entire site preserves every distance,
so the two additional tetramer sites reproduce A and B, respectively.

| Site | 306 OD1–308 OD1 | OD1–OD2 | OD2–OD1 | OD2–OD2 | Minimum |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 4.0376 | 4.6759 | 4.4276 | 3.9163 | **3.9163** |
| B | 4.0287 | 4.6616 | 4.4190 | 3.9110 | **3.9110** |

The mature Fig. 10 deposited reference therefore does **not** contain a short
Asp306–Asp308 oxygen contact. A drawn relay or the precursor discussion's
suggestion of strong Asp–Asp hydrogen bonding should not be reported as an
observed short O–O bond in this structure. These distances neither establish
Asp charge states nor exclude rearrangement in a bound/reacting state.

Reproduction: `python3 analysis/negoro_2023_protonation/audit_crystal_and_mapping.py`.
The authoritative JSON stores every atom coordinate, pair distance, residue
mapping, source URL, timestamp and source/script SHA256:
`analysis/negoro_2023_protonation/5Y0M_asp_pair_and_mapping_audit.json`.
No molecular dynamics, protonation changes or energy calculations were run.

## Complete Table 1 structure panel

The follow-up panel uses every Table 1 PDB identifier: **5XYG, 3AXG, 5XYO,
5XYT, 5Y0L, 5Y0M, 7YU0, 7YU1 and 7YU2**. This includes the mature wild types
and thermostability variants, the Ala266/Ala267 and Cys267 uncleaved precursors,
and the hydroxylamine-cleaved Cys267 variant. The exact variant labels, all
resolved sites, every oxygen pair, deposited coordinates and source SHA256 are
retained in
`analysis/negoro_2023_protonation/table1_job/panel/table1_asp_pair_summary.json`;
the companion CSV provides compact site-level comparisons.

The extraction ran in tmux with a retained script/source snapshot, input hashes,
command, environment and process log under
`analysis/negoro_2023_protonation/table1_job/`. Process completion was exit 0.
Across measured sites, the smallest oxygen-label-invariant Asp306/Asp308
distance is **3.5794 Å (5XYO chain B)** and the largest is **4.1736 Å (7YU2
chain A)**. None of the measured site minima is below 3.5 Å. This broader panel
therefore also provides no observed short Asp–Asp oxygen contact, while still
not determining protonation.

| PDB | Variant/state | Measured deposited sites | Range of site-minimum O–O distances (Å) |
| --- | --- | ---: | ---: |
| 5XYG | NylCp2 wild type, mature | 2 | 4.0985–4.1535 |
| 3AXG | NylCA wild type, mature | 15 | 3.8092–4.0669 |
| 5XYO | NylCp2 D122G, mature | 2 | 3.5794–3.8284 |
| 5XYT | NylCp2 H130Y, mature | 2 | 3.8701–4.0311 |
| 5Y0L | NylCp2 D122G/H130Y, mature | 2 | 3.8618–3.9588 |
| 5Y0M | NylCp2 D122G/H130Y/D36A/E263Q, mature | 2 | 3.9110–3.9163 |
| 7YU0 | NylCp2 H130Y/N266A/T267A, precursor | 2 | 4.1220–4.1417 |
| 7YU1 | NylCp2 D122G/H130Y/T267C, precursor | 2 | 3.9752–4.0344 |
| 7YU2 | NylCp2 D122G/H130Y/T267C, hydroxylamine-cleaved | 2 | 4.0995–4.1736 |

Coverage: **31 deposited sites** were measured. An independent ATOM-record
inventory found a complete OD1/OD2 pair for both target residues in every
deposited 3AXG protein chain. Its 15 measured chains therefore do not indicate
a missing Asp oxygen in an otherwise deposited chain. No additional sixteenth
site is assumed or reconstructed. Symmetry-related copies are not counted as
independent observations. The updated extractor records deposited chains and
target-atom inventories; the original panel and source snapshots are preserved.

## Implication for the project

Treat D272/D274 protonation and Lys156 participation as explicit competing
microstate/mechanism hypotheses. The paper motivates examining them but does
not select a validated mature Nyl12 microstate. Preserve the project's measured
product and TI1 connectivity gates; audit the retained Thr233 N-terminal relay,
Asn186/F113 backbone donors, and local water geometry before selecting further
calculations. A short contact in a project model must be compared with its own
proton inventory and trajectory, not labeled crystal agreement on the basis of
the mature Fig. 10 drawing.
