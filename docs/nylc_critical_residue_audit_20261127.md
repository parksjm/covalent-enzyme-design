# NylC critical-residue evidence audit

This audit distinguishes mature-substrate hydrolysis from precursor autocleavage,
expression, oligomer assembly, and thermostability. The principal sources are
Negoro et al. 2023 [DOI:10.1111/febs.16755](https://doi.org/10.1111/febs.16755),
2012 [DOI:10.1074/jbc.M111.321992](https://doi.org/10.1074/jbc.M111.321992),
and 2018 [DOI:10.1038/s41598-018-27860-w](https://doi.org/10.1038/s41598-018-27860-w).
In the supplied 2023 PDF, Table 1 is p. 3 (journal p. 3402), the K189A
discussion/Fig. 8 are pp. 11–12 (journal pp. 3410–3411), and Fig. 10 is p. 17
(journal p. 3416).
See the [machine-readable evidence table](../analysis/negoro_mutagenesis_mapping_20261127/nylc_critical_residue_evidence.json)
and its CSV rendering for construct backgrounds, numbers, mappings, and citations.

## Catalytic and maturation evidence

| NylC → Nyl12 | Primary evidence | Interpretation / confidence |
|---|---|---|
| N266 → N232 | N266A: 42–57% cleavage; N266G: 40–57%. | Strong cleavage-site maturation sensitivity; no independent mature-substrate role. |
| T267 → T233 | T267C: <1% cleavage; T267S: 9.7% initially, >99% after 48 h, cleavage-corrected activity 0.72 U mg⁻¹. The reported T267A result is N266A/T267A with H130Y, not an isolated point mutant. | Strong maturation role; substrate nucleophile role plausible but inferred because impaired maturation confounds low activity. |
| D306 → D272 | D306A in NylCp2 and H-GYAQ: <1% cleavage, <0.1 U mg⁻¹ raw activity. | Strong maturation requirement; mature-substrate triad role remains unisolated by these mutants. |
| D308 → D274 | D308A in H-GYAQ: <1% cleavage, <0.1 U mg⁻¹ raw activity. | Same maturation confound; substrate-triad role is proposed. |
| G307 → G273 | Backbone NH proposed to help stabilize oxyanion; no side-chain mutagenesis. | Proposal only; low confidence in residue-specific necessity. |
| K189 → K156 | H-GYAQ K189A: 6–33% cleavage (67–94% remains uncleaved); cleavage-normalized activity 0.032 U mg⁻¹ (~0.6% parental); trace Ahx-dimer after 24 h. | Strong hydrolysis-associated loss beyond maturation; exact chemical role remains proposed. |
| Y146 → F113 | H-GYAQ Y146A: cleavage-normalized activity 0.18 U mg⁻¹ (~3.6% parental); NylCA Y146A kcat/Km ~3% of parent. | Supports Tyr contribution, but Ala removes ring and OH; no Y146F result. Nyl12 F113 retains backbone NH, lacks Tyr OH. |
| N219 → N186 | N219A cleavage 31–33%; mature hydrolysis was not independently quantified. | Strong maturation-network effect; proposed mature oxyanion donor, not proven by this mutant. |
| F134 → F101 | F134A cleavage 44–59%; activity not tested. | Active-site-adjacent maturation effect, not established catalytic necessity. |

Table 1 of the 2023 paper reports the main mutation values; its footnote d
corrects partially processed-enzyme activities by cleavage ratio. The assay
uses Ahx cyclic oligomer. Lys189 and Tyr146 effects therefore bear on this
substrate system, not a demonstrated Nyl12 substrate register or product
profile. Lys189's general-base/network role and the Asp308–Asp306–Thr267
substrate triad are proposed by structure/MD/mechanistic analogy; they are not
direct observations of proton transfer.

## Other functionally altered residues are not catalytic essentials

L137A (NylC→Nyl12 M104) disrupts NylCp2 assembly and processing: <1% cleavage,
<0.1 U mg⁻¹, and Tm 41 °C versus 52 °C WT. D122G partly rescues processing
(33–34%). However, NylCA naturally carries A137 and is >99% cleaved, showing
context dependence. This is an assembly/maturation warning, not evidence that
L137 is chemically essential. D122G/H130Y/D36A/E263Q are interface/stability
engineering substitutions (GYAQ Tm 88 °C); they are active but are not
demonstrated essential catalytic residues. D122 and H130 substitutions also
change activity and substrate affinity, so these are not stability-only
phenotypes. See 2012 Tables 1–2 and 2018 Tables 1–3.

The 2018 processability screen additionally reports expression, aggregation,
or processing phenotypes at positions 75, 106, 141, 157, 179, 191, 235, 295,
and 299. They do not identify additional catalytic-essential residues. Mapped
positions and outcome classes are retained in the machine-readable table.

## Mapping and design consequence

Mappings use independent residue counts through the existing alignment of
NylC_A/3AXG (355 residues) to the current Nyl12 concatenated alpha/beta query
(328 residues); no constant offset is applied. N232 is the last alpha-chain
residue, and T233 is the first beta-chain residue. These are homolog mappings,
not proof of identical geometry. The alignment and hashes are identified in
the JSON receipt; source numbering for NylCp2 constructs is presumed to match
the common 355-residue NylC numbering and is qualified there. Some p2
experimental alleles differ from the 3AXG reference allele at the same site:
p2 D122/H130/L137/V225 align by position to 3AXG G122/Y130/A137/M225 and to
Nyl12 S89/N97/M104/L192, respectively. The machine table records the actual
3AXG reference identity separately to avoid implying identity at those sites.

Among the presently hard-protected positions, the direct NylC evidence is
strongest for K156's hydrolysis-associated role and for N186/N232/T233/D272/D274
in maturation or the proposed catalytic network. K156 protection is a
precautionary NylC-based analogy; there is no direct Nyl12 K156-mutant proof.
F113 preserves the backbone NH but not Tyr's OH. M104/F101 are protected in
the current mask for recurring inter-repeat interface reasons, with NylC
mutagenesis supporting context-dependent structural effects rather than
chemical essentiality. This audit does not justify protecting every residue
that changed activity, stability, or expression in another NylC background.
