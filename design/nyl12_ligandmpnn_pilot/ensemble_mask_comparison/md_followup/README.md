# Matched dynamics of two designs and two controls

The PI requested productive GPU follow-up while the window033 Hessian runs.
This bounded campaign tests persistence of modeled binding and preorganization,
not chemical barriers or experimentally improved activity. Eight trajectories
are planned: two independent 100 ns velocity replicas for each of four sequences.
Each sequence has one replica on Violet and one on Azalea. Only those hosts are
assigned. Existing QM/MM and other calculations are preserved.

## Starting structures and pose policy

| Sequence | Relaxed model seed | Initial signed phi, AE/BF/CG/DH (degrees) | Initial tiers | Native gate passes |
|---|---|---|---|---|
| WT | 2026112204 | −105.5 / −100.6 / −101.8 / −101.5 | A/A/A/A | 0/4 |
| Broader mask_07 | 2026112204 | −100.0 / −118.7 / −101.1 / −113.8 | A/B/A/B | 2/4 |
| Focused mask_03 | 2026112203 | −113.0 / −100.4 / −116.2 / −98.8 | B/A/B/B | 2/4 |
| V111C | 2026112202 | −92.3 / −95.9 / −91.6 / −92.6 | A/B/A/A | 3/4 |

The PI requires **flipped-only starting models**: all four assigned scissile
amides have negative signed projected Herschlag phi under the existing mapping.
No production restraint holds this sign. A later pose change is a result, not a
reason to delete frames. Positive-phi historical data remain preserved.
The initially proposed mixed-pose WT seed03 was never launched; its prepared
copy and superseded selection are retained in `clutter/`.

One structure per sequence is deliberately selected from existing models.
Starting native-Thr conformation and catalytic geometry differ. These are
protocol-matched retention tests, not a matched-coordinate mutation experiment
or an unbiased prediction of equilibrium occupancy. WT remains a control despite
its native-gate failures. The previously active K59D filter failure prevents use
of this gate as an activity classifier. Report sites nested within trajectories,
with velocity replicas as the independent dynamic replicates.

## Inputs and protocol

Reuse the exact accepted, fully solvated minimization endpoints; no extra
prediction, parameter fit, re-solvation or minimization is performed. The input
adapter verifies native command/input hashes, process completion, force tolerance,
released ligand parameters, symmetry/protected identities, hydrogen inventory,
and equivalence to the previously hashed water/ion-stripped endpoint. Derived
portable receipts retain all upstream receipts and explicitly identify their
provenance. Visualization-export hashes refer to stripped artifacts and are never
misrepresented as full solvated topology/restart hashes.

All four substrates per tetramer are PA6 L4 zwitterions, using the pinned modular
RESP/GAFF2 release with ff19SB/TIP3P and 0.15 M salt plus sequence-specific
counterions. All four catalytic Thr N termini are NH2; Asp272 and Asp274 are ASP.
No protonation or force-field changes are introduced. This force field models
noncovalent E–S, not TI1 or reactive proton transfer.

The existing `pilot_md_protocol.py` and prepared-MD runner are frozen under
`source_snapshot/`: 100 ps heating, 100 ps density relaxation, 100 ps weak-restraint
relaxation, 200 ps unrestrained equilibration, then 10 + 90 ns production at
300 K. Temporary heavy-solute positional restraints are 5, 2, then 0.5
kcal/mol/angstrom² during the three preparation stages only. Equilibration and
production use `ntr=0,nmropt=0`; all atoms move. Production is 2 fs, SHAKE on
H bonds, Langevin thermostat and Monte Carlo NPT at 1 atm, PME with 10 angstrom
cutoff. Trajectories are saved every 10 ps. No reaction coordinate is biased.

## Execution and acceptance

`configs/` contains frozen dataset/replica/host/GPU assignments and independent
stage seeds. `systems/` contains self-contained inputs and portable receipts.
`md/SEQUENCE/replica_0N/dynamics/` holds native inputs, outputs and trajectories.
`queues/HOST_gpuN/` contains persistent parent tmux receipts; each native replica
also has its own tmux process receipt. A failed stage stops its queue; the other
GPU queues remain independent. Occupied GPUs are rejected rather than shared.
Azalea uses its previously validated Blackwell CUDA12.8 executable, never the
generic Ada binary. Exact paths, hashes and the known failure are in
`docs/AMBER_WORKSTATION_RUNTIME.md`. Its corrected configs have the explicit
`azalea_blackwell_` prefix; failed generic-engine plans remain in `clutter/`.

On the assigned host, using its verified scientific Python, the reproducible
entry point is `launch_host.py --host HOST`; add `--prepare-plans` only when that
host has no existing plans. Existing outputs are never overwritten. Native
completion and temperature checks do not establish sampling convergence or
catalytic competence. Inspect early temperature/density/finite-energy behavior
and later replicate agreement before interpreting geometry.
The entrypoint verifies the exact host-specific executable hash and config
before launching anything. A missing or incorrect profile fails explicitly.

The analysis deliverable is per-site signed phi, d1/d2/d3, Bürgi–Dunitz angle,
native-Thr conformation, and substrate-retention time/block distributions, with
all-frame and bound-frame denominators reported separately. Also assess assembly
integrity and Asp contacts. The required four-ligand streaming analysis adapter
is a follow-up implementation task; no existing single-ligand result is reused
as a four-site analysis, and analysis is not yet automatically queued.

See `selection.json` for exact source paths and `sequence_comparison.fasta` for
full precursor and processed-chain sequences. New reaction-path work still
depends on the independent window033 curvature result.
