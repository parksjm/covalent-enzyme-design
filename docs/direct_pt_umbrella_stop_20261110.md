# Direct OG1→ThrN proton-transfer umbrellas: intentional stop

Stopped on violet at 2026-11-11T01:50:27.692609+00:00 under the user's authorization to end
simulations without a useful immediate scientific goal. This supersedes the
handoff instruction to continue this campaign to 30 ps/window.

## Scientific decision

The standalone coordinate is d(OG1,HG1)−d(ThrN,HG1), using DFTB3/3OB QM/MM.
The 7-ps PMF remains essentially uphill: final umbrella target about 15.49
kcal/mol above the first, without a transferred-proton minimum. Sparse matched
embedded DFT correction is also uphill, about 15.85 kcal/mol at the endpoint.
These endpoint costs are not accepted activation barriers.

The latest analyzed prefix is **12 ps/window**, with 48,000 raw and 28,268
decorrelated samples; overlap scalar 0.01281, minimum bidirectional adjacent
overlap 0.18255. Changes from 7 ps are 0.414 kcal/mol RMS and 0.562 maximum.
Sampling is connected but not established as converged. Further sampling could
improve standalone PT precision but cannot establish coupled nucleophilic
addition, leaving-group protonation and acylation or alternative-pose chemistry.
No immediate precision target justified continuing that cost under the user's
current priority. This is a resource/prioritization decision, not proof that
standalone transfer is impossible or TI1 does not exist.

## Process acceptance and recovery

Controller PID 2008323 was verified by parent, full command, output, cwd and
procfs start time; 20 child Sander processes were verified by parent, executable,
per-stage cwd, input path, full command and start time. SIGTERM stopped the
controller first, preventing later block launches, then each verified remaining
child. No SIGKILL was needed. Subsequent process inventory found no local
Sander or umbrella-extension process and the target tmux session
`nyl12_direct_umbrella_30ps_queue_20261109` was absent. No remote signal was
issued; azalea NEB and unrelated jobs were untouched.

All original artifacts remain under
`qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/07_direct_pt_umbrella_30ps_violet/`.
Original campaign/process status files remain historical evidence: stale
`running` values or wrapper signal exit status do not establish spontaneous
failure. Supplemental authoritative classification:
`provenance/direct_pt_umbrella_stop_20261110/intentional_stop.json`.
The same directory retains the stop script, exact pre-signal process inventory,
original status snapshot, signal outcomes, per-window progress and SHA256
values of all last completed restarts. All 20 hashes were reverified.

Interrupted-block trajectories and restarts require validation before use and
are not marked completed. Any authorized continuation must use a validated
retained restart and a new attempt directory. Do not rerun the original
launcher or stop script. Nothing was deleted or overwritten by the stop script.

## Final per-window progress

Every window completed 5-ps extension equilibration and four 5-ps production
blocks: **20 ps extension production + 2 ps pilot = 22 ps/window complete raw
data**. This exceeds the latest analyzed 12-ps prefix. No final 22-ps MBAR or
all-block chemistry audit was performed during stopping; no converged PMF is
claimed. Block 5 was interrupted. Steps below are last written NSTEP values,
not an exact count of all executed steps or certified usable checkpoints.
Each step is 0.0005 ps.

| Window | Complete extension production (ps) | Interrupted block | Last output step |
|---|---:|---:|---:|
| window_00 | 20 | 5 | 400 |
| window_01 | 20 | 5 | 300 |
| window_02 | 20 | 5 | 200 |
| window_03 | 20 | 5 | 300 |
| window_04 | 20 | 5 | 400 |
| window_05 | 20 | 5 | 300 |
| window_06 | 20 | 5 | 200 |
| window_07 | 20 | 5 | 300 |
| window_08 | 20 | 5 | 400 |
| window_09 | 20 | 5 | 100 |
| window_10 | 20 | 5 | 300 |
| window_11 | 20 | 5 | 400 |
| window_12 | 20 | 5 | 300 |
| window_13 | 20 | 5 | no NSTEP written |
| window_14 | 20 | 5 | 200 |
| window_15 | 20 | 5 | 100 |
| window_16 | 20 | 5 | 300 |
| window_17 | 20 | 5 | no NSTEP written |
| window_18 | 20 | 5 | 100 |
| window_19 | 20 | 5 | 400 |

Checks: exact process identities, termination, session absence, and all 20
retained completed restart hashes. No production code or dependencies changed;
no pytest run was needed for the operational stop. Journal:
`output_agent_actions/20261021.jsonl`.
