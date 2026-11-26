# Direct-US continuation: 30 to 60 ps per window

All 63 original windows continue from their accepted 30 ps restarts on Violet.
Forty single-CPU tmux workers dynamically claim one 5 ps block at a time,
prioritizing complete common sampling horizons. New blocks 07–12 are in each
existing window directory; original prepared.json and accepted data remain
unchanged. Only Langevin seeds change between the cloned production inputs.
Current 16-file 3ob parameter, engine and topology hashes are pinned in plan.json.

Native commands/statuses are in worker_jobs/worker_XX; private worker event logs
and claims are in workers/. Failed windows remain explicit in receipts/.
Root review and immutable source/plan manifests are retained here. The absence
of a worker session does not establish scientific success; inspect receipts.

A separate derived-only watcher writes direct_35ps_analysis.json through
60 ps in analysis/fixed_burnin_20261126/. It concatenates continuous production,
discards the first 5 ps once, keeps all subsequent frames and performs 200
within-window contiguous 5 ps/2.5 ps block-bootstrap MBAR refits plus discard
and temporal sensitivities. It plots comparisons only after the same-horizon
relay report exists. Plots stay flat in the campaign plots/ directory and
pmf_latest.png advances to the latest published direct profile. Native biases
remain their original q_acyl/q_H potentials; no chemical frames are censored.

All current US jobs remain running. No NEB was launched. Completion of 60 ps
is not convergence; assess effective support, temporal/orthogonal stability,
proton inventory, endpoint basins and targeted opposite-start controls.
