# Water-free TI1 optimization trajectory

Snapshot: 14 geometry frames, 39 atoms and 23 explicit H per frame.

On your laptop, change PyMOL's working directory to the downloaded folder containing this README (normally nyl12-inspection-20261021/qm/waterfree_ti1_optimization). Then run:

```pml
@view_optimization.pml
mplay
```

Use mstop or the frame slider to inspect individual states. VMD can load optimization_trajectory.xyz as a multi-frame XYZ file.

The movie contains complete Input orientation tables from the immutable log snapshot, including any trial/repeated geometries. These are not MD times or guaranteed accepted optimizer steps. It is not an optimized minimum. All hydrogens are explicit; water is intentionally absent from this39-atom model.

PDB CONECT records retain the initial chemical graph for display and do not prove bonds survive. Check trajectory_metrics.csv for the four reactive bond distances, minimum oxyanion-H distance, frozen-atom shifts, chi1 and graph screens.

The exporter can be rerun with a new --output directory while Gaussian runs. Each export reads the log once, keeps that snapshot and never overwrites previous exports.
