# Asp272 two-stage HG1 shuttle pilot

User explicitly authorized this construction hypothesis on 2026-10-25.
Active run: `output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2/`.
tmux: `nyl12_asp272_relay2_20261025`. Read status.json and per-window records
for progress; initial native optimizer evaluations were verified at launch.

Stage 1 uses 17 coupled targets: shorten OG1–C54 from its initial distance
to 1.55 Å while increasing d(OG1,HG1)−d(Asp272 OD1,HG1) from its measured
initial value to +0.8 Å. OD1 is selected because it accepts the observed
hydroxyl hydrogen bond. All other QM coordinates remain free.

Stage 2 starts only after a converged, geometrically plausible covalent
intermediate with HG1 on Asp272. It holds the attack target at 1.55 Å and
uses 17 targets for d(Asp272 OD1,HG1)−d(terminal N,HG1), ending at +0.8 Å.
The original HG1 identity and both original terminal N–H bonds are tracked.
The intended final state has HG1 on terminal N, an oxyanion, and restored
Asp272 carboxylate. These assignments require final coordinate inspection.

The 103251-atom solvated system, 98 original QM atoms plus eight links,
charge −2/singlet, DFTB3/3ob-3-1 Hamiltonian and original MM restraint
reference/constant are retained. No water is added and no classical TI1
parameters are introduced. Start from the original ES endpoint, not the
failed direct-N forward scan or the stalled strict-ES endpoint.

Use native Sander `ntmin=5` with DL-FIND Cartesian L-BFGS and analytic
gradients. The maximum gradient-component tolerance is 1e-4 Hartree/Bohr
(approximately 0.119 kcal/mol/Å), together with the native RMS-gradient,
step and energy-change checks. This replaces the inadequate global-RMS-only
acceptance and avoids the prior conjugate-gradient line search. All atoms
remain active. The native optimizer is being tested here; convergence is
not guaranteed. A failed optimizer stops the sequence without a fallback.

Ordinary harmonic NMR restraints use rk2=rk3=100 kcal/mol/Å² in Amber's
kΔ² convention. Each window must report native convergence, finite coordinates,
the unchanged box and both achieved coordinates within 0.15 Å of target.
Chemical gates check the covalent bond, C–O/C–N geometry and proton identity
before advancing between stages. No reverse scan, MD or TS search is queued.

Keep ordinary DISANG separate from the optional `&dlfind disang` hard-
constraint converter; that converter is not used. NMR reads its weight and
redirection records first, so `&dlfind` follows those records in the input.
Attempt 1 failed before gradients because of the original namelist ordering
and remains preserved. Attempt 2 reached L-BFGS energy/gradient evaluations.

ntwx=5 records DL-FIND evaluation callbacks, including trial evaluations;
these are not necessarily accepted optimizer iterations or physical times.
Completed windows also export the exact endpoint and a stripped DCD with
the fixed protein/substrate/three-water selection. Original full-system
trajectories remain unchanged. The construction scan is not itself a
validated reaction barrier. Keep restraint energies explicit in later plots.

Specification: configs/asp272_relay_pilot_20261025.json. Source, input hashes,
environments, logs and process/scientific records are retained per attempt.
Ruff passed and the first real native optimizer evaluations were checked;
test coverage was not measured for this one-off experimental controller.
