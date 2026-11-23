# Native-Thr233 gate calibration — 2026-11-23

The unchanged gate is reproducible and compatible with every available
crystallographic catalytic site. It is **not validated as an activity classifier**.
No threshold, protonation state, atom mapping, MD trajectory or existing design
score was changed. All structural counts below are descriptive.

## Frozen scientific definition

After Kabsch fitting backbone N/CA/C/O of Phe113, Asn186, Asp272 and Asp274
(16 atoms), hold out Thr233. Require OG1 displacement ≤0.75 Å, Nα displacement
≤1.25 Å, |wrapped ψ difference| ≤35°, and |wrapped χ1 difference| ≤20°.
The reference is **auth chain A** of the frozen `source_inputs/n12c.cif` in
`qmmm/pa6_flipped_model30/10_wt_protonation_md/analysis_after_md/`.
Thr233 is Amber929; adjacent Thr234 is Amber930. ψ uses N233–CA233–C233–N234;
χ1 uses N233–CA233–CB233–OG1. This distinction excludes the neighboring Thr
from the held-out identity. We load the frozen
`catalytic_alignment_geometry.py` directly, and retain the source files
`nalpha_water_relay.py` and `catalytic_alignment.py` for the mapping/gate audit.

Hypothesis: the fixed geometric envelope should accept equivalent crystal
sites across independently refined room/cryo structures. This tests reference
compatibility and reproducibility, not whether an active substrate-bound
state must remain within that apo/reference envelope.

## Crystal calibration and copy identity

Original coordinates were downloaded from the primary RCSB files
[9Z0A](https://files.rcsb.org/download/9Z0A.cif) (room temperature, 1.90 Å) and
[9Z0B](https://files.rcsb.org/download/9Z0B.cif) (cryo, 1.75 Å). No coordinates
were minimized or regenerated. Author residue numbers define mature Thr233;
label numbers are retained explicitly, never substituted for author numbers.
All required residues are present in each counted author chain; none of these
protein chains was excluded. Residues231/232 are absent from both original
structures, while Thr233/234 are present. Each counted author chain contains
both α-domain Phe113/Asn186 and β-domain Thr233/234/Asp272/274. These are full
mature-site author chains, not eight separate α/β fragments.

| Coordinate source | Author chains; label chains | Frozen gate passes |
|---|---|---:|
| Supplied n12c derivative | A,F,H,V; A,B,C,D | 4/4 |
| Original cryo 9Z0B asymmetric unit | A–H; A–H | 8/8 |
| Original room 9Z0A asymmetric unit | A,B; A,B | 2/2 |

There are **two diffraction datasets, ten original asymmetric-unit sites**;
these are not ten independent biological experiments. Assembly-generated
symmetry images were not expanded or counted. The four supplied n12c sites
are derived, not a third experiment or four additional independent sites.
The supplied V site corresponds geometrically to original cryoG, with symmetry
placement absorbed by Kabsch alignment; it is not a novel independent copy.
Original assembly-generation operators remain in the retained RCSB files.

| Maximum over the ten original crystal sites vs n12cA | Observed | Fixed limit |
|---|---:|---:|
| Fitted backbone RMSD | 0.1254 Å | no gate |
| Held-out OG1 displacement | 0.1232 Å | 0.75 Å |
| Held-out Nα displacement | 0.0841 Å | 1.25 Å |
| Absolute ψ difference | 3.244° | 35° |
| Absolute χ1 difference | 6.834° | 20° |

The supplied reference is an edited derivative, not an original deposited
structure: n12cA/F/H/V map to cryoA/F/H/G with local backbone RMSD
0.00554/0.00548/0.00437/0.00515 Å. Sequence audit across shared observed
residues gives 313/313,309/309,313/314 and312/312 identities, respectively.
The sole shared-residue mismatch is **n12cH ALA228 versus original 9Z0BH PRO228**.
It lies outside the gate residues, but prevents describing the supplied
coordinates as an unmodified experimental reference. Full author/label maps,
observed sequences and selected alternate conformers are retained in JSON.
Alternate conformers use highest occupancy, then blank/A/name to break ties;
this explicit extraction policy does not alter the frozen gate.

Paired sensitivity uses the same 14 sites (including derived n12c) and the
same thresholds. All14 pass with either n12cA or room 9Z0AA reference under
(i) original 16-atom fit, (ii) 20-atom fit including Thr233 backbone, or
(iii) 8-atom α-domain-only Phe113/Asn186 fit. Largest OG1/Nα shifts are:

| Reference | Original fit | Including Thr backbone | α-domain-only fit |
|---|---:|---:|---:|
| n12cA | 0.123/0.084 Å | 0.114/0.074 Å | 0.264/0.202 Å |
| Room9Z0AA | 0.147/0.081 Å | 0.139/0.072 Å | 0.412/0.271 Å |

This bounded sensitivity supports crystal compatibility despite reference/fit
choice. Including held-out Thr in fitting is a diagnostic only: it would bias
its displacement downward and does **not** replace the original fit.
No corresponding sensitivity claim is made for substrate-bound MD frames.

## Known-active controls: existing static endpoints

The 55-model occupied-site CSV was used directly: four assigned sites per
model, five seeds per sequence, raw and relaxed stages. All 440 stored native
flags reproduce exactly; no ligand×site duplication from `all_sites.csv` was
introduced. Control activity is established by the project assays, not these
models. Seeds/sites are correlated structural outputs, not activity replicates.

| Control | Raw native | Relaxed native | OG1 fails | Nα fails | ψ fails | χ1 fails | Relaxed joint TierA/B+native |
|---|---:|---:|---:|---:|---:|---:|---:|
| WT | 20/20 | 5/20 | 15 | 0 | 1 | 0 | 1/20 |
| K59D | 20/20 | 6/20 | 14 | 0 | 0 | 0 | 0/20 |
| V111C | 20/20 | 16/20 | 4 | 0 | 0 | 0 | 7/20 |

WT's single ψ failure overlaps an OG1 failure. Relaxed OG1 medians are
0.794/0.842/0.699 Å (WT/K59D/V111C); ranges are 0.295–1.061,
0.184–1.261 and0.512–0.900 Å. Raw OG1 medians are 0.141/0.125/0.133 Å.
Thus much control rejection appears after relaxation and is specifically
held-out OG1 motion. K59D is a direct warning against treating a failed
native/joint filter as inactivity. Raw agreement also does not validate
activity: prediction/reference similarity can share structural priors.
Full continuous distributions and failure combinations for all 11 sequences
are in `results/control_summary.csv` and authoritative `calibration.json`.

## Existing matched WT MD: local scaffold breathing versus Thr displacement

Twelve previously produced occupancy CSVs were summarized; trajectories were
not read. All 120000 stored native flags reproduce. Each state has 40000
correlated sampled frames across four replicas. Counts are pooled descriptive
occupancy, not independent Bernoulli trials or confidence bounds.

| State | Native passes | FitRMSD 5/50/95% (Å) | OG1 fails | Nα fails | ψ fails | χ1 fails |
|---|---:|---|---:|---:|---:|---:|
| NH2/Asp272/Asp274 | 7 (0.0175%) | 0.242/0.383/0.555 | 39959 | 39976 | 21961 | 429 |
| NH2/Asp272/Ash274 | 98 (0.245%) | 0.295/0.461/0.948 | 39472 | 37088 | 17244 | 9947 |
| NH3/Asp272/Asp274 | 0 | 0.311/0.447/0.903 | 39822 | 39906 | 18072 | 361 |

NH2/Asp272/Asp274 OG1 and Nα median shifts are 2.450 and3.019 Å;
NH2/Asp272/Ash274 gives 2.722 and3.403 Å; NH3 gives 2.104 and2.738 Å.
The gate displacements are already measured **after a local rigid-body fit**,
so global translation/rotation does not explain them. Local scaffold breathing
exists: fitted-backbone RMSDs exceed crystal values. It is not equivalent to
Thr motion, however: median fitRMSD for native failures is 0.383/0.461/0.447 Å;
rare native passes give 0.459/0.460 Å for the two NH2 states. FitRMSD correlates
negatively with OG1 displacement (r≈−0.318/−0.348/−0.670 in the table's state
order), so the scalar fit residual does not explain rejection as a simple
larger-fit-error artifact. These correlations are descriptive and mix replicas;
they do not demonstrate a causal motion mechanism or validate the cutoff.

## Decision and reproducibility

Retain all four thresholds. Report this as a **crystal-reference-calibrated
geometry descriptor, activity-unvalidated**. Keep native and joint gates
separate, retain continuous metrics and individual failure reasons, and keep
known-active failures visible. Do not fit new cutoffs to preferred designs.
Room/cryo agreement does not establish that substrate-bound, protonated,
solution-state active conformations must resemble these reference sites.
Before using rejection for activity triage, require a prospectively specified,
matched known-active/inactive substrate-bound validation set and assay labels.
No such campaign is launched or authorized by this report.

Artifacts under `analysis/native_thr_gate_calibration_20261123/`:

- `inputs/`: original downloaded CIFs, exact frozen n12c and occupied-site CSV.
- `source_snapshot/`, `diagnostic_source_snapshot/`: immutable run-specific code.
- `process/`, `process_diagnostics/`: exact commands, hashes, environments,
  runnable scripts, logs and exit0 statuses from two bounded tmux analyses.
- `results/calibration.json`: authoritative crystal/reference/control analysis;
  CSV projections, `crystal_calibration.png`, separate scientific acceptance.
- `diagnostics/diagnostics.json`: authoritative existing-MD scalar summary and
  sequence audit; CSV projections and `known_active_raw_to_relaxed.png`.
- `environment.json`: explicit software versions, no new production packages.

Code is isolated in `native_thr_calibration.py` with two script entry points.
Ruff check/format passed; focused analysis-environment pytest **13 passed**, coverage **97% total**
(module 96%, calibration CLI 100%, diagnostics CLI 98%). Tests check original gate
boundaries, invalid metrics, auth/label mapping, rigid-motion invariance,
held-out perturbation, exact stored flags, duplicate detection and both CLI
artifact/acceptance outputs. Full-suite testing was not repeated for this
isolated analysis addition. Original artifacts and the first run remain intact.

Biopython is an explicit prerequisite for the CIF extraction entry point and
is already installed in the recorded analysis environment; it is not a new
project production dependency. The module imports it only when `read_sites`
is called and fails explicitly if unavailable. Default-environment tests skip
only the three Bio-dependent CIF cases; scalar gate/control/MD tests still
run (**10 passed, 3 skipped**). The full analysis-environment run passes all
**13 tests** with **97% total new-code coverage**. No replacement parser or silent scientific fallback is used. Immutable
completed-run source snapshots retain the parser originally executed.
