# TI1 parameter generation protocol — 2026-10-22

## User-selected reference and scope

The user authorized generation and validation now, then specified a truncated
ACE/Thr/NME anion with a retained tetrahedral modification and constrained
backbone/side-chain torsions. This is the first **bonded-parameter surrogate**.
It is distinct from the real N-terminal Thr NH3+/oxyanion TI1: do not transfer
its charges to that native state or declare the native force field ready.

Use the unmodified user-selected `cif/n12c.cif`, SHA-256
`e318fddc7142be1efabd6ebbd80ad6d45bf98f6d83f1ea7c1e60add9646d7b8c`.
The explicit selection is author chain A, Thr233. The mmCIF has one model and
multiple chains; other chains include sequence variants and are not pooled.
Measured psi = 134.110276428 degrees, chi1 = -58.964417878 degrees.
No preceding chain-A carbonyl exists, so native phi is undefined. ACE-cap phi
values -120 and -60 degrees are declared model choices for cap sensitivity.
Both R and S attacked-carbon configurations are explicit comparison states;
neither is silently declared the biologically preferred face.

## Local graph

Use ACE-THR(-O-C*(O−)(CH3)(NHCH3))-NME: Ac-Thr-NME plus an attacked
N-methylacetamide surrogate, with the Thr hydroxyl proton omitted to an
external reservoir. Formula C10H20N3O4−, 37 atoms, charge -1, 134 electrons,
singlet; Thr CA S and CB R. The C*–N bond remains intact. Methyl termini
represent the local PA6 amide environment; they are not a full L4 substrate.
Relay water is absent from this local bonded reference; the actual enzyme
QM subsystem will still include explicit relay water with TIP3P MM naming/LJ.

The native full-L4 graph is separately scripted: 103 atoms including explicit
water, NH3+ reacting Thr, substrate oxyanion, charge 0. It is reserved for
native-state charge/transfer validation. Real PA6 termini remain NH3+/COO−.
The current provisional protein topology remains barred from refinement.

## Calculation stages

1. Deterministic ETKDG seed construction, then set declared phi/psi/chi1.
   This constructs starting coordinates; it is not an energy optimization.
2. Gas-phase wB97XD/def2TZVP constrained optimization in Gaussian 16, with
   UltraFine integration, Tight optimization/SCF, CalcFC and XQC. No extra
   dispersion correction. Freeze only the three declared dihedrals; the
   covalent TI1 reaction-center bonds and all hydrogens remain free.
3. Check normal completion, ordered atoms, L-Thr and C* stereochemistry,
   covalent connectivity/proton inventory and frozen-torsion agreement before
   starting the fresh DFT frequency/Hessian calculation.
4. Obtain the exact Cartesian DFT Hessian in a separate checkpoint, followed
   by a separate HF/6-31G(d) MK ESP calculation/checkpoint. Surrogate RESP,
   when fitted, is for testing this surrogate; native charges require the
   actual NH3+/oxyanion coupled state and explicit cap-transfer constraints.

Constrained references are not free minima. Raw full-molecule frequencies
include constrained directions; their signs do not alone establish a minimum
on the mobile-coordinate manifold. Preserve gradients, constraints and raw
Hessians for conditional-curvature analysis before fitting/validation.

## Screens declared before fitting

Optimization screening requires the intended charge/multiplicity, all input
stereocenters/connectivity, no detected proton transfer and frozen torsions
within 0.10 degree. Reactive bond windows (angstroms) are conservative chemical
screens, not fitted equilibrium targets: OG1–C* 1.30–1.75, C*–O− 1.25–1.65,
C*–N 1.25–1.70, C*–C 1.35–1.75. Retain every failure; do not relax these gates
after observing results. Record O− contacts to cap/backbone H explicitly;
fixed dihedrals do not guarantee absence of artificial gas-phase H bonding.
The scripted checker explicitly rejects O−–H distances below 1.20 angstrom
as a conservative possible-proton-transfer screen, in addition to checking
all original bonds and newly short nonbonded contacts. This is a screening
criterion, not a fitted interaction distance.

Candidate bonded terms will be derived from QM data with documented Amber
harmonic conventions. Existing AmberTools modified-Seminario routines are
available; source/units must be checked before using them. No new dependency
is needed or authorized. Independently compare candidate MM structures and
conditional curvature/scans with held-out QM data, including the new ether
junction, oxyanion and intact amide. Prospective local geometry tolerances:
bond errors <=0.03 angstrom; angle RMS <=3 degrees and maximum <=5 degrees.
Compare cap orientations: local bond targets within 0.03 angstrom and force
constants within 20%; larger sensitivity requires a larger reference model.
Held-out scan energies should have RMSE <=1 kcal/mol within the tested local
basin. These are project acceptance criteria, not literature-certified bounds.

Native-state RESP must separately enforce exact total charge (residual <1e-6),
document standard peptide/cap constraints, use shared multi-conformer charges,
report fit/held-out ESP errors and retain the true reactive proton inventory.
Do not force arbitrary integer partial charges on the two covalently joined
residues. Complete transfer/geometry/torsion/curvature checks and exact-topology
digests before declaring parameters eligible for AF3-model MM relaxation.

All preparation, engine calculations and real-input checks use scripted tmux
jobs, source snapshots/hashes, explicit resources/environments, process logs
and separate scientific acceptance records. Current action journal remains
`output_agent_actions/20261021.jsonl`. AF3 fork changes remain the user's work.
