## ACTIVE SMD/water-bridge diagnostic — 2026-10-23

User requested a decision and tmux launch before leaving in5min.
New root: output_nyl12_smd_water_bridge_20261023/.
ACTIVE tmux: nyl12_bridge63_smd_20261023; pilot svp_pilot_01/.
63atoms, neutral singlet, wB97XD/def2SVP SCRF=(SMD,Solvent=Water),
16cores32GB,20cycles maximum, same frozen21heavy plus chi1.
No reactive-bond/proton/water restraints during optimization; no Hessian/RESP.

Prior water was not a validated bridge. New deterministic search moves ONLY
water O/H/H, retains every nonwater coordinate, and screens graph/clashes.
1800 orientations evaluated; selected N-Owater=OG1-Owater=2.70A,
N-H3-Owater=151.38deg, Owater-H-OG1=180deg, OH0.96A/HOH104.5deg.
Initial63atom geometry passes. Final acceptance adds both water heavy-atom
contacts2.5--3.3A and best NH-water / waterH-OG1 angles>=140deg, alongside
existing NH3, oxyanion, graph/stereo and convergence checks. A retained
state would be a constrained aqueous-continuum cluster diagnostic, not a
native enzyme minimum or accepted FF reference. Water continuum is not a
measured enzyme dielectric. Do not autoextend or fit parameters.

Ruff passes and real seed/source geometry checks pass. New bridge-gate unit
coverage has not been measured under the user's5min deadline. Earlier12tests
apply to the preceding gas-phase pipeline, not these added water criteria.
Source snapshots, input hashes and action journal retained. Job manifest pins
launched assessor; separate assessment_source_hash_at_checkpoint.json discloses
preformat digest in earlier staging source_manifest. No method silently changed.

Literature: Meadows etal2026 JCTC22,4012--4025 DOI10.1021/acs.jctc.6c00162
explicitly uses implicit water to stabilize glycine zwitterion for GAFF/RESP
parameter variants: https://pmc.ncbi.nlm.nih.gov/articles/PMC13130867/ .
Enzyme SMD dielectric4 precedent: https://doi.org/10.3389/fchem.2018.00249
(singlepoint enzyme energies; not proof SMD Hessians transfer to Amber).
Recommendation: solvated geometry is defensible; ESP charge protocol and
bonded fitting remain separate decisions. Do not adopt fully solvated ESP
charges or raw continuum Hessian as automatic ff19SB/GAFF2 parameters.
No replacement RESP or force-field fit launched. Diffuse/larger-basis and
actual protein-environment checks remain pending if TI1 survives.

