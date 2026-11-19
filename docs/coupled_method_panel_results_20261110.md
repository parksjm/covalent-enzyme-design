# Coupled enzyme Hamiltonian benchmark: results

All six states completed DFTB3/3OB, PM6-DH+, GFN2-xTB and
omegaB97X-D/def2-SVP evaluations on azalea. The wrapper exited zero;
analysis independently required Sander energy completion, normal Gaussian
termination and normal xTB termination. These are fixed-geometry comparisons,
not optimized method-specific pathways or stationary-point classifications.

Errors below are double differences relative to the original E--S geometry,
in kcal/mol. Positive values mean the approximate method overestimates the
relative energy versus DFT.

| State | DFTB3 | PM6-DH+ | GFN2-xTB |
|---|---:|---:|---:|
| E--S | 0.00 | 0.00 | 0.00 |
| Direct PT | -0.36 | -2.24 | +8.02 |
| Coupled PT (window 24) | -1.72 | -12.22 | +2.20 |
| Attack | -1.15 | -10.08 | -0.76 |
| Late attack | +2.17 | -13.21 | -2.43 |
| TI1-like | +7.42 | -14.17 | -1.57 |
| Six-state RMSE | 3.27 | 10.26 | 3.61 |

From coupled PT to TI1-like, DFTB3's error grows by 9.13 kcal/mol;
xTB's changes by -3.78 kcal/mol. Thus the earlier isolated-PT result did not
validate DFTB3 for attack or oxyanion formation. DFTB3 biases this branch against
TI1-like geometry relative to DFT. xTB is the better candidate for further
attack-region testing, but its +8.02 kcal/mol direct-PT error prevents blanket
promotion for the complete mechanism. PM6-DH+ is not promoted: its large
systematic stabilization of coupled states is inconsistent with the reference.
None passes every previously proposed energy gate.

Point-charge membership varies from 8,104 to 7,999 across these structures.
Raw finite-cutoff energy changes are therefore not accepted reaction energies.
Per-geometry method comparisons share the Amber embedding; double differences
reduce common electrostatic offsets but do not eliminate method-dependent
polarization or cutoff sensitivity. A fixed-membership embedding sensitivity
test and selected diffuse/larger-basis DFT checks are still needed. The reference
is an anionic QM region with link atoms and a modest basis, not an exact answer.

Next: evaluate reference/approximate gradients along attack and proton-transfer
directions, test whether the late-attack/TI1-like seed retains its chemical state
on unbiased DFT/MM relaxation, and validate xTB/MM integration before production
sampling. Hydrolysis does not require a stable tetrahedral minimum: a concerted
path remains possible. No barrier, stable TI1, or validated fast Hamiltonian is
claimed by this benchmark.

Artifacts: `qmmm/benchmarks/04_coupled_path_method_panel/analysis/` contains
`comparison.json`, `comparison.csv` and `method_errors.png`. The script is
`scripts/analyze_coupled_method_panel.py`. Ruff passed. The first execution under
uv failed because matplotlib is absent; the existing analysis pixi Python ran
successfully without installing dependencies. No new unit coverage was measured.

User preference: parallelize independent benchmark states whenever this reduces
time to solution. Azalea has 28 physical cores/56 threads; use physical-core
budgets (for example two 12-core Gaussian workers) and respect per-state xTB
embedding dependencies. Preserve completed work when changing scheduling.
