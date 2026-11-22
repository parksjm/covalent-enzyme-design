# Dynamic NEB-seeded umbrella-sampling movie

Open `view_umbrella.vmd` in VMD. The portable bundle requires `path.parm7` and
`path.dcd` in the same directory. All three are included by the laptop pull script.

The verified movie contains 30 snapshots from each of 63 windows: five evenly
spaced frames from each accepted 5 ps production block (1,890 frames total).
Frames are chronological within each window; windows follow the original NEB
path order. Window boundaries are discontinuities between independent biased
simulations, not physical reaction events. The JSON/CSV frame manifests identify
every source host, trajectory, window, block and original frame.

Bulk water and ions are stripped; any atom in the declared QM subsystem is
retained. The protein cartoon is available but initially hidden. The QM subsystem
is in licorice; only MM substrate carbons are orange. Reactive bonds use dynamic
cutoffs (1.4 Å for H–N/O and 1.85 Å for the reactive heavy-atom bonds). HG1 is a
small sphere. Oxyanion-hole distances and the Bürgi–Dunitz angle are labeled.

The full trajectories, topology, restarts and accepted-block records are preserved
on their original workstations. The US campaign is closed after 30 ps/window;
this visualization does not start or extend simulations. The latest PMF remains
preliminary and unconverged (see `../analysis/pmf_latest/`).

Export code, commands, logs, source hashes and per-host parts are retained under
`source_snapshot/` and `work/`. Intermediate host parts are excluded from the
laptop pull so only the final combined trajectory is downloaded.

Native verification: 1,890 frames × 19,967 atoms; DCD 453,002,956 bytes.
All four saved reaction distances agree for every frame within 0.00005 Å.
The movie assembly and bundle checksums passed. The VMD GUI was not available
for an interactive visual check on this host.
