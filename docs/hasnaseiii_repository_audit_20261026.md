# hASNaseIII repository audit

## Provenance

- Repository: <https://github.com/emedio/hASNaseIII>
- Audited commit: `8d1c076bc72393d4d7d794a8cf1d72bc92b3ff42`
- Commit date: 2023-07-14
- Audit date: 2026-10-26
- Related paper: *Elucidation of the Active Form and Reaction Mechanism in
  Human Asparaginase Type III Using Multiscale Simulations*, J. Chem. Inf.
  Model. 2023, 63, 5676–5688.

The repository has no root README, executable workflow scripts, or license
file. Do not copy its inputs into this project. The mechanistic setup can still
inform an independently implemented, cited workflow.

## Transferable contents

`QMMM_inputs/template_in` couples Amber QM/MM to Gaussian 16 through
`qm_theory='EXTERN'`. It specifies B3LYP with `IOp(3/124=30)` and
6-31+G(d), 48 Gaussian threads, a 15 A QM cutoff, 1 fs Langevin dynamics at
310 K, and two million steps per window. These settings document the authors'
production method; they are not automatically appropriate for Nyl12.

`QMMM_inputs/CVs` defines six bond-distance collective variables:

1. catalytic Thr HG1 to the Thr N-terminal N;
2. catalytic Thr HG1 to Thr OG1;
3. Thr OG1 to the substrate carbonyl carbon;
4. substrate carbonyl O to carbonyl C;
5. substrate leaving-group N to Thr HG1; and
6. substrate leaving-group N to carbonyl C.

`QMMM_inputs/STRING` configures a string/PMF calculation. `group_file` lists
96 windows, seeded in mechanistic blocks from reactant, first-proton-transfer,
nucleophilic, and second-proton-transfer structures. `guess` provides six
initial nodes in the six-dimensional collective-variable space. The repository
also provides seven transition-state PDBs and movies for acyl formation and
hydrolysis.

## Relevance to Nyl12

The work concerns substrate hydrolysis by the mature enzyme, so its reaction
coordinates are more relevant than coordinates for autocatalytic precursor
cleavage. It supports three changes to the Nyl12 reaction-path strategy:

- use a coupled, multidimensional coordinate description instead of forcing
  proton transfer and nucleophilic attack through one distance at a time;
- track carbonyl rehybridization and scissile C–N cleavage explicitly, in
  addition to Thr proton transfer and OG1–C bond formation; and
- treat pathway free energies with a string/umbrella method after stable
  reactant and intermediate basins have been established.

It does **not** establish that Nyl12 follows the same direct Thr-to-N-terminal
proton transfer, that Asp272 is absent from the Nyl12 coordinate, or that its
B3LYP setup is the best Nyl12 Hamiltonian. The Nyl12 coordinate set must retain
Asp272 proton-shuttle alternatives until calculations and structural evidence
exclude them.

## Recommended use

Complete the current unrestrained Model 23 DFTB3/MM E-S optimization first.
Then define an independently implemented Nyl12 multidimensional path with, at
minimum, OG1–HG1, HG1–acceptor, OG1–substrate C, substrate C–O, and substrate
C–N coordinates. Add Asp272 O–H and H–N coordinates for the relay hypothesis.
Use DFTB3 for tractable path exploration, validate representative stationary
points and proton-transfer energetics with range-separated-hybrid DFT, and use
DL-FIND only after a defensible approximate saddle geometry is available.

## Supporting Information geometry benchmark — 2026-10-27

The official ACS Supporting Information was retrieved from ACS Figshare record
24041393 (file `ci3c00900_si_001.pdf`, DOI
`10.1021/acs.jcim.3c00900.s001`, CC BY-NC 4.0). It describes restrained
heating followed by gradual restraint release and three independent 1 us
unrestrained Michaelis-complex simulations. It does not provide Cartesian
coordinates for a Michaelis-complex snapshot. Its data-availability statement
instead points to the GitHub repository, which deposits transition-state PDBs.

The corresponding CB–OG1–substrate-carbon angle was measured from the authors'
deposited structures, using TH2 residue 478 atoms CB/OG1 and substrate ZN619
atom CG:

| Structure | CB–OG1–CG angle | OG1–CG distance | Interpretation |
|---|---:|---:|---|
| TS1 | 135.211 degrees | 3.365 A | Thr proton transfer precedes close attack |
| TS2 | 116.949 degrees | 2.470 A | Nucleophilic-attack transition state |

The TS2 value supports using a near-tetrahedral CB–OG1–C orientation as a
mechanistic diagnostic. It does not justify a universal hard cutoff: TS1 is
outside 105–120 degrees, and the Nyl12 Model 23 start and degraded DFTB3/MM
endpoint both fall inside that interval. Combine this angle with attack
distance, Burgi–Dunitz angle, projected phi, Thr chi1/psi, and crystal-aligned
Thr displacements.
