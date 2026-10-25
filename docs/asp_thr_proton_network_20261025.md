# Asp–Thr proton network: observed geometry and next hypotheses

User inspection identified the Asp272 (Amber residue 968) contact with
catalytic Thr233 (Amber 929). An independent full-system geometry audit
confirms that Asp272 OD1 accepts contacts from both the terminal amine H1
and the hydroxyl proton HG1.

| Contact, H-bond donor–H–acceptor | Original ES H···A / angle | Failed forward endpoint |
|---|---:|---:|
| N–H1···Asp272 OD1 | 2.153 Å / 170.1° | 1.932 Å / 164.9° |
| OG1–HG1···Asp272 OD1 | 1.984 Å / 157.3° | 2.146 Å / 134.3° |
| OG1–HG1···terminal N | 2.570 Å / 96.3° | 1.841 Å / 131.8° |

This suggests a competing, well-oriented Asp proton acceptor. It does not
prove that the amine cannot rotate or accept a proton. An N–H donor contact
does not itself consume the amine lone pair. The contact could orient the
amine and stabilize a protonated state. N lone-pair orientation is not directly
shown by the VMD bond representation. Continuing the imposed direct-N scan
could force transfer, but that would not establish the preferred mechanism.

No water currently hydrogen-bonds to these sites: with triclinic periodic
images, nearest water O distances in the original ES are 8.70/6.44 Å to the
two Asp272 oxygens and 10.45/9.87 Å to Thr N/OG1. All 27690 WAT residues
were checked; this is not a stripped-display selection artifact. These
coordinates do not establish equilibrium hydration or exclude slower entry.
There is no geometrically justified new relay-water site from this audit.

Recommended comparison, after fixing minimizer convergence: direct transfer
to the terminal N versus an Asp272-assisted route, allowing the other QM
atoms to respond. The latter would require its own explicit proton-state
assignment, not automatic acceptance as the previously prescribed NH3+/oxyanion
TI1. Keep both Asp residues; do not remove their interactions merely to make
the prescribed route proceed. Do not add a water solely to lower an energy.

Related experimental evidence: the NylC study reports an Asp308–Asp306–Thr267
triad and loss of autocleavage for the two Asp-to-Ala substitutions. These
are the paper's numbers; sequence correspondence to Nyl12 must be verified
before transferring residue-specific mechanistic assignments.
[Primary study, PMID 36799721](https://pubmed.ncbi.nlm.nih.gov/36799721/).

The tighter ES job stopped after 223 cycles with REPEATED LINMIN FAILURE:
RMS 0.016484 versus target 0.0005, GMAX 0.85962. Exit zero did not establish
convergence. Its contacts remain close to the original ES values; no new
barrier or TI1 stability result is available. No continuation was launched.

Audit script, hashes and geometry JSON:
`output_nyl12_pa6_dftb_pilot_20261024/asp_thr_review_20261025/`.
VMD shared viewer now overlays only non-QM substrate carbons in orange,
keeps QM Licorice and retains protein NewCartoon as a default-hidden rep.
The original ligand rep is retained underneath to preserve heteroatom bonds.
