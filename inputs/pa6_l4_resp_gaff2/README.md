# PA6 L4 zwitterion RESP/GAFF2 inputs

The primary state is the explicit 79-atom PA6 L4 zwitterion: N56 is NH3+,
C53/O48/O51 is COO-, and the total charge is zero. This is the declared model
for a solvent-exposed soluble oligomer near neutral pH; it is not a claim that
termini buried in bulk or crystalline nylon have the same protonation state.

The DFT geometry stage freezes the 25 listed heavy-atom backbone and amide
torsions while optimizing all other coordinates with omegaB97X-D/def2-TZVP and
SMD water. These constraints preserve each conformer's identity and prevent an
artificial intramolecular ion pair or proton transfer during isolated-cluster
optimization. The later HF/6-31G(d) Merz-Kollman calculations use the accepted
DFT coordinates without geometry optimization. Frequencies are not required
for this charge-sampling stage; bonded-parameter benchmarks are separate.
