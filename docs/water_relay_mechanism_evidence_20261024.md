# Water-relay hypotheses and limits of the packing result — 2026-10-24

The user correctly distinguishes a remembered relay mechanism in other enzymes
from an established relay in NylC. They now question whether front placement is
required. Primary literature was broadened to N-terminal nucleophile hydrolases
and threonine proteases. These are competing mechanistic models, not proof that
all family members use the same general base or relay geometry.

## Primary sources

- Wei D, Fang L, Tang M, Zhan CG. Fundamental reaction pathway for peptide
  metabolism by proteasome: insights from first-principles QM/MM free energy
  calculations. J Phys Chem B (2013)117,13418–13434.
  DOI10.1021/jp405337v; https://pubmed.ncbi.nlm.nih.gov/24111489/ .
  Particularly relevant positive example: the reported preferred substrate-
  hydrolysis pathway activates Thr1 via water-assisted proton transfer during
  acylation. It distinguishes this shuttle from the water attacking the
  acyl-enzyme in deacylation. These are calculations on a proteasome peptide
  substrate, not NylC or nylon. The abstract and primary indexed text were read;
  PMC full-page retrieval intermittently returned CAPTCHA.
- Crystal structures of penicillin acylase enzyme-substrate complexes:
  structural insights into the catalytic mechanism (2001).
  https://www.sciencedirect.com/science/article/pii/S0022283601950430 .
  Primary abstract describes enhancement of the N-terminal Ser nucleophile
  by its alpha-amino group, potentially mediated by a bridging water. Proposal
  based on enzyme/substrate structures, not proof of a required relay in NylC.
- Chilov GG, Sidorova AV, Svedas VK. Quantum chemical studies of the catalytic
  mechanism of N-terminal nucleophile hydrolase. Biochemistry(Mosc)(2007)
  72,495–500. DOI10.1134/S0006297907050057.
  https://pubmed.ncbi.nlm.nih.gov/17573703/ .
  Penicillin-acylase quantum model supports activation by its own terminal
  amine without a required bridging water. Counterexample to universality;
  the specific solvation environment still affects the calculated reaction.
- A unified mechanism for proteolysis and autocatalytic activation in the
  20S proteasome (2016). DOI10.1038/ncomms10900.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4792962/ .
  Describes historical Thr-OH-to-terminal-amine water-shuttle proposals and
  re-examines the general base using structural/mutational evidence involving
  Lys33/Asp17. Do not present older water proposals as settled consensus.
- Negoro et al. X-ray crystallographic and mutational analysis of the NylC
  precursor: catalytic mechanism of autocleavage and substrate hydrolysis of
  nylon hydrolase. FEBS J(2023)290,3400–3421. DOI10.1111/febs.16755.
  https://febs.onlinelibrary.wiley.com/doi/10.1111/febs.16755 .
  Substrate-hydrolysis section proposes Asp308–Asp306–Thr267, including the
  terminal amino group, and an oxyanion hole. Its water-accessible/water-
  inaccessible MD analysis concerns PRECURSOR AUTOCLEAVAGE. Do not conflate
  those water-relay figures with established mature-nylon hydrolysis.

## What our geometry says

The failed front search rejects the simultaneous assumptions of the current
fixed AF3-derived protein/substrate conformation, the specified front hemisphere,
and a single water contacting both terminalN and OG1 within2.6–3.2A under the
stated clash thresholds. It is not a reaction-path or free-energy calculation.
It cannot choose between direct transfer, another water position, another
protein/substrate conformation, or a carboxylate-assisted network.

The source pocket is unrefined; it may be an incorrect reactive conformation.
Conversely, a relay may operate during activation/TS1 without the same water
remaining as a closed bridge in the TI1 minimum. Water used later to hydrolyze
an acyl-enzyme has a separate chemical role and need not be this same molecule.
Gas-phase TI1 collapse did not establish the enzyme's mechanism or disprove
TI1; it may reflect missing stabilization, boundaries, or an incorrect state.

## Next comparisons; not yet adopted chemistry

1. Check candidate water locations without a front-side restriction; same
   full-pocket steric screens apply. A second geometric search is authorized
   by the resumed placement investigation and the user's reconsideration.
2. Map NylC literature residues to the actual Nyl12 sequence before assigning
   Asp272 a catalytic correspondence. Its proximity makes its electrostatics
   worth examining rather than treating it only as a steric obstacle.
3. Compare direct terminal-amine transfer, physically accessible water-assisted
   transfer, and a mapped carboxylate-assisted network in consistent models.
   Review protonation, total charge and boundaries before adding an Asp to QM.
4. Check native crystal geometry and substrate pose. Do not deform the protein
   solely to accommodate the assumed ring or resume fitting a rejected TI1.

No new mechanism, residue mapping, protonation state, water-free force field,
Hessian or RESP protocol is accepted by this literature review. No new Gaussian
job was launched. Computational geometry search results are documented separately
in docs/front_water_feasibility_20261024.md.
