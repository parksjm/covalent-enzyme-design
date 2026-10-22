# TI1 capped-reference progress snapshots

These four 37-atom ACE/Thr/NME anions are unconverged Gaussian fitting references, not enzyme structures or accepted minima.
Each contains all 20 explicit hydrogens. No relay water belongs to this small bonded-fitting surrogate.
For each case, cd into its directory and run `pymol view_progress.pml`. Use the movie controls to inspect all saved coordinate tables.
VMD can open `trajectory.xyz`. `latest_geometry.xyz` and `.sdf` contain the latest table, which may be an unevaluated trial.
SDF connectivity is the intended graph, not a proof that QM preserved it. Consult `geometry_checks.json` and `trajectory_metrics.csv`.
Inspect the tetrahedral C*, intact OG1-C* and C*-N bonds, oxyanion proton inventory, and cap/backbone H bonds.
Only phi, psi and chi1 are frozen. Close cap H bonds are reported and may signal gas-phase cap sensitivity.
The XYZ tables are optimization snapshots, not dynamics or guaranteed one-frame-per-step output.
