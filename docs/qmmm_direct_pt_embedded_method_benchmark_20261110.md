# Direct-proton-transfer embedded method benchmark

## Scientific question

The direct Thr233 `OG1-HG1` to N-terminal-N coordinate drive is uphill by
15.28 kcal/mol after subtracting its small driving bias, and its unrestrained
endpoint returns to E-S. The completed gas-phase small-model benchmark suggests
that DFTB3/3OB can overestimate the transferred-state energy by 8.6--17.7
kcal/mol relative to omegaB97X-D. This benchmark tests that hypothesis on exact
enzyme-path geometries rather than transferring a gas-phase correction.

All 15 requested evaluations completed normally on 10 November 2026. In the
actual enzyme QM region and electrostatic environment, DFTB3 follows the DFT
relative-energy changes more closely than GFN2-xTB does. The result does not
support replacing DFTB3 with GFN2-xTB for this path.

## Fixed geometries and chemical state

Five accepted coordinate-drive endpoints are compared: `pt_01`, `pt_05`,
`pt_10`, `pt_14`, and `pt_19`. They span proton-difference coordinates of
-1.003, -0.676, -0.184, +0.242, and +0.718 A, respectively. The copied restart
hashes are recorded in `attempt1/work/benchmark_manifest.json`.

The established QM mask contains 98 real atoms. Amber generates eight link
hydrogens, so the Gaussian and xTB calculations contain 106 atoms. This direct
count corrects the older window-24 DFT/MM manifest, which incorrectly stated
six links and 104 total QM-plus-link atoms. The QM charge is -2 and the spin
multiplicity is one.

## Matched electrostatic embedding

All three methods use the same Amber QM mask, link-atom construction,
`adjust_q=2` charge redistribution, and periodically imaged MM charges within
20 A of the QM region. QM Ewald and QM PME are disabled for this comparison
because Amber's external-Gaussian interface does not support them. DFTB3 is
therefore reevaluated with the same finite embedding rather than compared
directly with the original Ewald coordinate-drive energies.

Amber writes each omegaB97X-D Gaussian input. The retained input is the source
of truth for both the 106-atom QM/link geometry and the adjusted external-charge
list. The xTB converter changes Gaussian point-charge records from
`x y z q` in angstrom to xTB legacy records `q x y z 999.0` in bohr. The
hardness value 999.0 requests the point-charge limit instead of xTB's default
dummy-atom hardness. A definition-only smoke test parsed all 106 atoms with
xTB 6.7.1 on azalea.

The actual 20-A charge counts are 8,104, 8,105, 8,107, 8,095 and 8,088 for
`pt_01`, `pt_05`, `pt_10`, `pt_14` and `pt_19`, respectively. Therefore, the
raw finite-cutoff relative profiles are not physical reaction profiles: changing
the point-charge membership changes the electrostatic-potential zero seen by
the net -2 QM region. The scientifically interpretable comparison is the
per-window DFTB3-minus-DFT or xTB-minus-DFT double difference. This cancels the
common window-specific offset to first order.

## Methods and interpretation

Each fixed geometry receives one energy-and-gradient evaluation with:

- Amber DFTB3/3OB/MM;
- standalone GFN2-xTB 6.7.1 with the exact Amber-generated point charges; and
- Amber external omegaB97X-D/def2-SVP/MM through Gaussian 16.

The primary comparison is the method-relative embedded-QM profile. DFTB3 and
DFT also provide the complete Sander hybrid potential. Until Sander is rebuilt
with the xTB library interface, the analysis constructs an xTB hybrid profile
by adding the DFTB3-derived relative MM remainder to the embedded xTB energy.
That constructed profile is explicitly approximate; it is not labeled as a
native xTB/MM Sander energy.

This is a five-point fixed-geometry Hamiltonian diagnostic, not a reaction
barrier, stationary-point calculation, or free-energy profile. A 15/20/25 A
embedding-cutoff sensitivity check with a fixed atom/image set remains necessary
before treating absolute method differences as quantitatively converged.

## Results

Relative to omegaB97X-D/def2-SVP, the five-window embedded-QM DFTB3 profile has
an RMSE of 1.49 kcal/mol and an endpoint error of -0.57 kcal/mol. GFN2-xTB has
an RMSE of 3.91 kcal/mol and an endpoint error of +7.93 kcal/mol. The DFTB3
errors at `pt_01`, `pt_05`, `pt_10`, `pt_14` and `pt_19` are 0.00, +2.07,
+2.34, +1.01 and -0.57 kcal/mol. These double differences were subtracted from
the original periodic-Ewald DFTB3/MM profile to obtain a diagnostic DFT-corrected
estimate.

The original periodic non-restraint QM/MM energies at the five windows are
0.00, 2.49, 9.61, 13.58 and 15.28 kcal/mol. The corresponding DFT-corrected
estimates are 0.00, 0.42, 7.27, 12.57 and 15.85 kcal/mol. Thus, the sparse
higher-level correction remains monotonically uphill and does not reveal a
transferred-proton minimum. This finding argues against the hypothesis that a
large DFTB3 energy error alone hid such a minimum. It does not prove that no
minimum exists on the full DFT/MM surface because gradients, relaxation and
Hessian classification were not performed at DFT/MM.

Mechanistically, an isolated Thr-to-N proton-transfer coordinate is probably
incomplete for this enzyme conformation. The next path work should test coupled
proton transfer and OG1 attack, or environment relaxation around a chemically
controlled coupled state, rather than rebuilding Sander for xTB solely on the
basis of the earlier gas-phase model.

## Reproducibility and status

- Benchmark root:
  `qmmm/benchmarks/03_direct_pt_embedded_method_comparison/attempt1/`
- Calculation host/session:
  `azalea`; `nyl12_direct_pt_embedded_methods_attempt1_20261110` (complete)
- Gaussian concurrency: three windows, 16 threads per window.
- Scratch:
  `/data/scratch/gaussian_scratch/jpc/nyl12_direct_pt_embedded_methods_attempt1_20261110`
- Job provenance:
  `attempt1/job/{command,status}.json`, `run.sh`, and `job.log`.
- Hashed source snapshot:
  `attempt1/source_snapshot/`.
- Parser-only retry result: `analysis_retry1/result/`.
- Final periodic-path correction and plot: `final_analysis/`.

All 15 engine calculations terminated normally. The first automatic analysis
failed after the calculations because it requested an `EAMBER (non-restraint)`
field not printed by these `nmropt=0` one-step evaluations. The separately
hashed parser-only retry correctly used `EPtot` and completed without rerunning
scientific calculations. Process and scientific outcomes remain separate in
the append-only action journal.

Final focused verification: Ruff lint and format checks passed; six focused
pytest tests passed; combined coverage of the two benchmark modules is 81%;
Bash syntax passed; and the xTB point-charge definition smoke test terminated
normally.
