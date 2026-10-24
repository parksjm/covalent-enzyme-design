# Front-side relay-water feasibility — 2026-10-24

The user requested a front-side water without clashes to Asp272 or other
residues, followed by parameterization/testing if feasible. No acceptable
placement was found with the source protein/substrate heavy atoms fixed.
No new Gaussian calculation or parameter fit was launched.

The source is the complete hydrogenated alpha/beta model in
output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation/.
The complete source substrate is retained in clash checks. The precise intact
63-atom cluster supplies cap hydrogen coordinates. The original water is
excluded and all three terminal ammonium H are excluded from the oxygen
screen, allowing optimistic reorientation. Other protein/substrate atoms
remain in place. Source coordinate frame agreement is asserted explicitly.

Front is the positive projection onto the user oxygen direction, after both
are projected perpendicular to the N-to-OG1 axis. This is independent of camera
orientation. Distance grid: 0.025 A; azimuth grid: 0.5 degrees. This samples a
hemisphere, not a proof about every continuous coordinate or other conformers.

| Screen | N/O-water distance range | Oxygen samples | Accepted |
|---|---|---:|---:|
| Full protein/substrate, all atoms except movable NH3 H | 2.6–3.2 A | 225,000 | 0 |
| Same geometry, ignore all H (diagnostic only) | 2.6–3.2 A | 225,000 | 0 |
| Expanded range (diagnostic only) | 2.6–3.3 A | 302,760 | 0 |

Bondi radii (A): H1.20, C1.70, N1.55, O1.52, S1.80. A contact fails when
its vdW overlap exceeds0.40 A. Explicit donor-H--waterO hydrogen-bond
exceptions require angle>=140 degrees, donor--O2.6–3.3 A and H--O>=1.5 A.
Intended N--Owater and OG1--Owater have minimum2.6 A. These are declared
steric screens, not force-field energies or scientific minimum acceptance.
The safe local environment cutoff follows the triangle inequality; distant
atoms cannot be omitted if they could overlap any searched oxygen position.

Best all-atom clearance is still0.684 A BELOW the allowed threshold, with
substrate H23 at1.636 A from waterO. The best heavy-only point still misses
by0.217 A, including contacts to substrate C70/C62 and Asp272 CA/OD1.
Expanded-range best margin is−0.594 A. Thus merely rotating ammonium/water
hydrogens cannot produce an accepted bridge in this fixed pocket.

A front bridge may become possible after limited protein/substrate heavy-atom
relaxation or in another validated conformation. This search does not prove
that an enzyme-compatible bridge is impossible. Any subsequent relaxation
must respect the user's requirement to validate TI1 parameters before using
them for protein refinement; do not quietly use the rejected legacy TI1 force
field to move the pocket. Native-structure comparison and an explicitly
validated relaxation/embedding protocol are the next scientific decisions.
SMD remains a proposed charge-stabilization diagnostic; it cannot repair a
steric overlap by itself. No RESP/Hessian/force-field transfer is accepted.

Artifacts: output_nyl12_front_water_20261024/. Source-pinned tmux jobs retain
search inputs, exact commands, process status and separate scientific reports.
inspection/view_front_search.pml shows full fixed pocket and two REJECTED
oxygen probes: user position and best sampled clearance. Neither is a water
seed. Asp Amber272 is chainE residue40; terminal Thr233 is chainE residue1.
The ordinary rsync package includes the viewer and diagnostic reports under
qm/front_water_feasibility/. No new optimization trajectory is claimed.

## Opposite-side diagnostic and completed checks

The user's reconsideration of front placement prompted the same oxygen screen
on the opposite hemisphere:0/225000 pass with all atoms and0/225000 with
heavy atoms only. Best all-atom margin−0.67204A (Asn186HD22), heavy-only
−0.41121A (Asn186ND2). Even allowing every potential additional N/O heavy
Hbond contact down to2.60A leaves the heavy-only necessary upper bound below
zero (−0.34121A). Merely changing hemisphere does not provide a validated
bridge in this fixed model. Opposite search is also a finite diagnostic, not
a barrier calculation or evidence that a water-mediated pathway is impossible.

Front search/diagnostic tests:8pass,97%/97.06%linecoverage, Ruffpass. Separate
SMD-assessor tests:8pass,93.94%linecoverage, including lostbridge and reverse
proton transfer. Final opposite tests are retained in their own job records.
See docs/water_relay_mechanism_evidence_20261024.md for the broader literature
and distinction between mature substrate hydrolysis and precursor autoprocessing.
