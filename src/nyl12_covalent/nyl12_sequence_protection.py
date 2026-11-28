"""Minimum Nyl12 sequence protections in the 328-residue precursor numbering.

The PI explicitly protects 221–232 and hydrolysis/maturation determinants.
This union retains historical beta233–235 and Asp272-through-Leu276
protections. Lys156 is protected following the PI's latest review: evidence
for the homologous NylC Lys189 motivates precautionary retention, without
claiming direct proof of Nyl12 essentiality. Leu235 remains
conservatively fixed as mature beta N-terminal context for this pilot; a
residue-specific catalytic role has not been established.
Legacy config assertions about 275–276 are protection policy,
not experimental proof that those residues directly catalyze autocleavage.
Additional structural/interface protections remain additive. Coordinate
mobility is a separate RFD3 selection; omission from a template never permits
sequence mutation at these positions.
"""

MINIMUM_PROTECTED_PRECURSOR = frozenset(
    {113, 156, 186, *range(221, 236), *range(272, 277)}
)
MINIMUM_PROTECTED_ALPHA = frozenset(
    position for position in MINIMUM_PROTECTED_PRECURSOR if position <= 232
)
MINIMUM_PROTECTED_BETA_LOCAL = frozenset(
    position - 232 for position in MINIMUM_PROTECTED_PRECURSOR if position > 232
)
