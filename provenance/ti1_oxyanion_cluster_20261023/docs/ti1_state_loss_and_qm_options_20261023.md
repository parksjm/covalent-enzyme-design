# TI1 state loss and QM cost decision — 2026-10-23

## Stopped calculation and measured mechanism

At the user's request, the water-free39 worker was identity-checked and sent
SIGTERM. Launcher ended at 16:20:27 EDT (20:20:27 UTC), exit1. All inputs,
log, checkpoint and scratch were preserved. No replacement QM job was started.
The user's image interpretation is supported by the 19-frame stopped log:

| Coordinate (A) | First frame | Last complete frame |
|---|---:|---:|
| Thr OG1–PA6 C54 | 1.3353 | 2.1393 |
| PA6 C54–O49 | 1.3427 | 1.2201 |
| Thr N–transferring H3 | 1.0100 | 1.4878 |
| Thr OG1–transferring H3 | 1.9931 | 0.9161 |

H3 has bookkeeping origin THR:HG1 in the atom map. The geometry supports
proton transfer from the ammonium to OG1, cleavage of OG1–C54 and recovery
of carbonyl character: a reactant-like state, not a TI1 parameter reference.
There is no converged endpoint or reaction barrier claim. Data/script/process
records: output_nyl12_waterfree_reference_20261023/state_loss_review_01/.

## The model needs the local stabilizing environment

The earlier recommendation that water was unnecessary for the primary local
reference was too categorical for this charged, reactive state. The standard
modified-residue recipe does not guarantee that an isolated reactive adduct
has the enzyme-stabilized TI1 minimum. Backbone freezes are not substitutes
for missing electrostatics or hydrogen bonds. A solvent/basis/functional change
that happens to retain the bond is not by itself physical validation.

A read-only scan of the original AF3-derived protein identifies two candidate
oxyanion donors omitted from the cluster:

| Source donor | H–O49 (A) | Donor–O49 (A) | D–H–O49 angle (deg) |
|---|---:|---:|---:|
| Chain A Phe113 backbone N–H | 1.7133 | 2.6966 | 163.36 |
| Chain A Asn186 ND2–HD22 | 2.1559 | 2.9776 | 137.25 |

These are source-model contacts, not experimentally validated donor assignments.
The next proposed cluster retains the short substrate and native Thr NH3+,
restores relay water, and explicitly includes these donor groups at consistent
protein-derived positions. Small peptide/side-chain amide fragments can supply
them without retaining the entire PA6 substrate: approximately63atoms if using
N-methylacetamide for the backbone donor and acetamide for the Asn side chain.
Exact cap construction and mapping must be checked before any run. Retain
protein boundary geometry; reaction-center bonds and transferring protons must
be free when assessing TI1 stability. Protein electrostatic embedding should
be explicitly defined for enzyme-state validation. Water alone does not replace
oxyanion-hole donors. No automatic restraints are added to force a TI1 minimum.

The next acceptance checks must track every intended ammonium N–H and competing
OG1–H distance, in addition to oxyanion–H, bond windows and stereochemistry.
Generic graph screening did flag later frames, but an O49–H-only transfer
coordinate missed the specific reverse transfer mechanism. This is a protocol
improvement proposal, not a silently changed acceptance rule for old data.

## Lower-cost QM recommendation; not yet adopted

First correct the physical model. Keep omegaB97X-D while evaluating a smaller
basis so functional and basis effects are distinguishable. def2-SVP is a
reasonable model-development option; diffuse augmentation should be tested
for the oxyanion (for example6-31+G(d) as a practical comparison). Plain
6-31G(d) lacks diffuse functions and is not the preferred charged-site default.
No numerical accuracy penalty or speedup is established without paired tests.
For the old39atom composition, spherical def2-SVP contractions would imply
339 functions versus the634 reported for def2-TZVP; this is an analytic size
estimate, not a runtime benchmark. Geometries and especially local force
constants need larger/diffuse-basis validation before adopting a cheaper fit
reference. Higher-level single-point energies alone do not validate a Hessian.

B3LYP-D3(BJ)/DZVP has high-quality force-field precedent in OpenFF Sage2.0.
DZVP is not the same basis as def2-SVP. That precedent does not make B3LYP
mandatory for Amber compatibility or establish accuracy for this charged TI1.
There is no need to switch both basis and functional or rerun every completed
reference now. Preserve the existing data as controls. HF/6-31G(d) RESP
remains a separate, previously selected charge protocol. No extra dispersion
may be added to omegaB97X-D. No accepted parameters or MD have been generated.

References checked:

- Weigend and Ahlrichs, PCCP2005,7,3297–3305. DOI10.1039/B508541A.
  https://pubs.rsc.org/en/Content/ArticleLanding/2005/CP/B508541A
- Chai and Head-Gordon, PCCP2008,10,6615–6620. DOI10.1039/B810189B.
  https://pubs.rsc.org/en/content/articlelanding/2008/cp/b810189b/unauth
- Boothroyd et al., JCTC2023,19,3251–3275. DOI10.1021/acs.jctc.3c00039.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10269353/ (QM method, section2.2.4).

## Implemented donor-cluster pilot — 2026-10-23

The user adopted the two simple amide donor fragments and explicit water,
then selected omegaB97X-D/def2-SVP for the next short stability test. The
smaller basis is now used for this development pilot; it is not an accepted
parameter-fitting basis. Existing TZVP results remain intact as controls.

The cluster has 63 atoms (C15H37N5O6), charge 0, singlet. Retain all 42
starting-core/water coordinates exactly; do not use the reactant-like endpoint.
NMA retains ILE112 CA/HA/C/O and PHE113 N/H/CA/HA. Replace the four cut
CA--N/CB or CA--C/CB directions with C--H caps at 1.09 A, giving CH3CONHCH3.
Acetamide retains ASN186 CB/HB2/HB3/CG/OD1/ND2/HD21/HD22; replace CB--CA
with C--H at 1.09 A, giving CH3CONH2. Preserved source donor geometry gives
the contacts in the table above. Protein frame matches the intact seed with
zero displacement; water matches the separately TIP3P-preconditioned source.
Coordinate-derived stereochemistry is Thr CA S, CB R, attacked C54 R.

All 9 added donor heavy atoms are frozen at those protein positions in the
first pilot. The original 12 scaffold heavy-atom freezes and AF3-seed chi1
remain; all 37 hydrogens, water O and reactive OG1/C54/O49/N59 are free.
There is no restraint on the reactive bond or transferring protons. This is
a held-donor cluster approximation without remaining protein electrostatics,
not an unconstrained gas-phase minimum or complete enzyme environment.
The 20-cycle pilot is an optimization only; no Hessian or RESP is launched.
A larger/diffuse-basis check remains planned if the state is retained.

Explicit new prospective proton gates require all intended ammonium N--H
lengths in 0.85--1.25 A and competing OG1--H lengths >=1.20 A; these
conservative geometry screens do not establish a barrier. Existing graph,
oxyanion H, stereochemistry, reaction-bond and convergence checks also apply.
A synthetic reverse-transfer corruption is rejected. Twelve focused tests
pass, with >91% line coverage of each of the four new preparation/assessment
modules. Original failed preparation attempts and a coverage-instrumentation
warning are retained with corrected successful runs; no failed result is hidden.

Parameter fitting remains blocked. The donor fragments and water represent
an environment; do not merge their charges into a modified Thr residue.
Future RESP fragment charge constraints/partition and bonded-force-constant
transfer need separate validation. No production dependencies were added.

Files: output_nyl12_oxyanion_cluster_20261023/{seed_01,svp_pilot_01}/.
Live tmux: nyl12_oxyanion63_svp_20261023; exact input, policy, run script,
job status and source manifests are retained. Laptop copies are under
qm/ti1_oxyanion_cluster/. No native TI1 or force field has been accepted.
