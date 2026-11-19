# Nyl12 QM/MM pathway registry

This is the single decision register for competing PA6-hydrolysis mechanisms.
New mechanism directories require a new registered hypothesis and a declared
promotion gate. Failed branches are retained for provenance but are not extended.

| ID | Hypothesis | Starting-state gate | Current status | Next promotion gate |
|---|---|---|---|---|
| M01 | Direct Thr233 HG1 transfer to the neutral N-terminal N, followed by OG1 attack | Accepted frame-1603 DFTB3/MM E-S endpoint; 12 Å movable shell | The 19-window proton-transfer drive completed, but unbiased endpoint optimization returned HG1 to OG1; a 20-window direct-PT umbrella pilot is running on violet | A converged PMF and overlap analysis must show whether a transferred-proton basin exists |
| M02 | Asp272 accepts HG1, then transfers it to the N-terminal N | Tier A/B E-S geometry plus a relay-ready Thr–Asp contact in the same structure | Preorganization screen complete on azalea; no coordinate drive launched | An unrestrained MM and 12 Å DFTB3/MM E-S optimization must retain Tier A/B geometry and OG1–Asp O = 2.5–3.2 Å with HG1–Asp O = 1.5–2.3 Å |
| M03 | Water-mediated proton relay | No accepted seed | Deferred | Literature support and a reproducible water placement that survives unbiased relaxation |
| M04 | Direct E-S to TI1 path without prescribing whether proton transfer or nucleophilic attack occurs first | Accepted frame-1603 DFTB3/MM E-S endpoint; 12 Å movable shell | The 24-window prefix reached proton-transferred geometry, but the continuation was invalid because H1 transferred from N back to OG1. The chemically contaminated US campaign is stopped. A corrected attack-only branch from valid window 24 is active on azalea with explicit HG1/H1/H2 identity gates. | Require the corrected branch to reach attack/TI1-like geometries without proton exchange, then pass the coupled DFTB3/PM6-DH+/GFN2-xTB versus omegaB97X-D benchmark before production 2D sampling |

## Shared controls

- Use the same PA6 RESP/GAFF2 topology, reviewed QM region, DFTB3 settings,
  12 Å DL-FIND movable shell, convergence thresholds, and endpoint diagnostics.
- Coordinate driving changes only the declared reaction coordinate. Do not add
  positional restraints or an OG1–carbonyl restraint during proton-transfer-only
  phases.
- After each driven phase, use restrained QM/MM MD only for environmental
  reorganization, then release every reaction-coordinate restraint for the
  stationary-point test.
- A biased endpoint or restrained-MD ensemble is not a PES minimum. Promotion
  requires an unbiased converged optimization that retains the intended state.
- Generate the same endpoint trajectory, reactive distances, and relative-energy
  plot for every promoted pathway.

## M04 asynchronous composite coordinate — 2026-11-09

The direct E-S-to-TI1 exploratory path uses one generalized Amber NMR
coordinate rather than independently scheduled proton-transfer and attack
restraints:

```text
xi = d(Thr233 OG1, HG1) - d(Thr233 N, HG1)
     - d(Thr233 OG1, substrate carbonyl C)
```

The mapped Amber atom indices are N = 13905, OG1 = 13916, HG1 = 13917 and
substrate carbonyl C = 19916. The corresponding generalized restraint uses
`rstwt = 1,-1,-1,0`. Because only the sum is biased, proton transfer and
nucleophilic attack may advance asynchronously; neither component has its own
window-by-window target. The initial value is approximately -4.317 Å and 24
windows advance it to -0.800 Å with increments no larger than 0.15 Å. The run
uses the accepted 12 Å DL-FIND movable shell and no positional restraints.

This coordinate intentionally does not determine the reaction mechanism. It
can permit stepwise or concerted motion, but it may also hide orthogonal
barriers or hysteresis. Component distances, energy and geometry are therefore
recorded for every window. The resulting path supplies seeds for umbrella
sampling and stationary-point searches; it is not itself a PMF, activation
barrier or transition-state validation.

An earlier trial used separately scheduled proton-transfer and attack
restraints. It was stopped after two accepted windows because it imposed an
unvalidated synchronous timing assumption. Its complete files and job
provenance are retained under the run's `clutter/superseded_synchronous_drive/`
directory and will not be extended.

The post-window-24 continuation exposed a degeneracy not controlled by this
scalar coordinate: H1, rather than the tracked HG1, transferred from the
terminal N back to OG1. Therefore windows after 24 and their umbrella samples
are not an E--S-to-TI1 path. The corrected branch starts from the exact window
24 restart, drives only OG1--carbonyl-C attack, and applies flat-bottom identity
gates to all three terminal N--H bonds. Its purpose is to create chemically
valid coupled benchmark geometries; it does not declare a TI1 stationary point.
The production free-energy calculation will use the explicit two-coordinate
space `xi_PT` and `d(OG1,C)` rather than the degenerate scalar.

The completed biased endpoint had OG1-HG1 = 2.275 Å, N-HG1 = 1.017 Å and
OG1-carbonyl C = 2.104 Å. Fully unbiased optimization converged but returned
to E-S: OG1-HG1 = 0.977 Å, N-HG1 = 2.073 Å and OG1-carbonyl C = 3.239 Å.
The coordinate-drive potential energy increased monotonically to 21.01
kcal/mol relative to its first window (20.81 kcal/mol after subtracting the
reported restraint term); no downhill approach to a TI1 basin was observed.

## M02 preorganization finding — 2026-11-09

The accepted frame-1603 DFTB3/MM endpoint is not relay-ready: its shortest
HG1–Asp272 oxygen distance is 4.718 Å. An azalea screen of all 2,000 WT
replica-3 MD frames found a relay-ready Asp968 OD1 contact in 1,927 frames by
the predeclared raw contact gate. However, replica 3 contained no frame passing
the existing native-Thr displacement gate. Representative relay-ready frames
1150 (Tier A) and 1476 (Tier B) have HG1–OD1 distances of 2.178 and 2.105 Å and
OG1–OD1 distances of 2.977 and 2.999 Å, respectively, but their Thr233 OG1
native displacements are 2.215 and 2.171 Å. Thus M02 currently represents a
distinct conformational basin and must not be forced from the M01 endpoint.
