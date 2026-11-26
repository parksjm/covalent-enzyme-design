# Concurrent acylenzyme hydration branches

| Run | Accepted starting time | New sampling cap | Native end time |
|---|---:|---:|---:|
| Main continuation (parent directory) | 51 ps | 50 ps | 101 ps |
| replica_02_from_26ps_extend04 | 26 ps | 25 ps | 51 ps |
| replica_03_from_41ps_extend07 | 41 ps | 25 ps | 66 ps |

All run on one Violet CPU each, in tmux, at 300 K with the same DFTB3/MM model,
proton inventory, 12 A movable region and retained complete amine fragment.
Each has a distinct Langevin stream and continues accepted velocities without
reheating. Old trajectory prefixes are shared aliases for chronological context;
they are not independent samples and must never be summed across replicas.
Branch statistics use only its newly generated segment names in replica_plan.json.

Each 5 ps block must pass native/SCC/finite-state/chemical/boundary checks, then
refresh water analysis and the stripped VMD movie. A water with O-C <=3.5 A and
O-C-O angle95–125deg triggers review; these are screening criteria only. Waters
remain MM and cannot perform hydrolysis in these runs. A candidate needs QM
promotion and reoptimization before a deacylation path can be computed.

Each branch has replica_plan.json, replica_execution.json, job/status.json,
analysis/report.json and vmd_progress/view_hydration.vmd. The latest VMD movie
initially contains only the borrowed prefix; the first new frames publish after
the next complete accepted 5 ps block. Native data and source/input hashes remain
on Violet; the laptop pull copies curated movie/status/analysis outputs.
