# Existing-MD leaving-group relay screen — 2026-11-18

Completed bounded read-only search of existing mature Nyl12–PA6 L4 E–S MD.
The primary state was ASP272/ASH274 with neutral catalytic Thr233 NH2;
ASH272/ASP274 and ASP272/ASP274 were controls with the same unchanged gates.
No new trajectories, protonation changes, invented water coordinates or
QM/MM calculations were used. No control microstate was promoted.

## Result

| Fixed microstate | Replica | Frames | Tier A/B, φ > 0 | All structural gates | One-water bridges |
| --- | ---: | ---: | ---: | ---: | ---: |
| ASP272/ASH274 | 1 | 10,000 | 5,865 | 8 | 0 |
| ASP272/ASH274 | 2 | 10,000 | 30 | 0 | 0 |
| ASP272/ASP274 | 1 | 10,000 | 819 | 0 | 0 |
| ASP272/ASP274 | 2 | 10,000 | 1,296 | 31 | 0 |
| ASH272/ASP274 | 1 | 10,000 | 7,253 | 12 | 0 |
| ASH272/ASP274 | 2 | 10,000 | 3,695 | 0 | 0 |
| Total | — | 60,000 | 18,958 | 51 | 0 |

There were 18,958 Tier A/B positive-φ frames and only 51 frames satisfying
all structural gates. None had the heavy-atom bridge, so leaving-leg and
strict two-leg hydrogen-bond counts were also zero. Primary `selected.json`
and `controls_attempt1/selected.json` report null primary/reserve selections.
No `selected.ncrst` or reserve restart exists because no candidate passed.
The primary selection wrapper was retained unchanged after controls.

Native displacement was the main bottleneck. Primary replica 1 failed OG1
and Nalpha displacement in 5,851 and 5,747 of its 5,865 Tier-eligible frames,
respectively (overlapping counts); medians were 2.033 and 2.348 Å. Its eight
fully eligible frames were 47, 58, 62, 75, 93, 95, 147 and 227.
Their Nalpha–Nleave distances were 4.061–4.563 Å. The best minimax observed
water was WAT22725 at frame 93: Nalpha–O 7.089 Å, O–Nleave 4.994 Å,
H–Nleave 5.137 Å, O–H···Nleave angle 76.10°. This is diagnostic only.

## Exact unchanged gates

The existing geometric thresholds were retained: Tier A required d1 ≤3.2 Å,
d2/d3 ≤3.8 Å, Bürgi–Dunitz angle 90–120°, and |φ| 70–110°; Tier B allowed
d1 ≤3.5 Å, d2/d3 ≤4.0 Å, angle 85–125°, and |φ| 60–125°. Only positive
signed φ frames were retained. Here d1 is Thr233 OG1–scissile C, d2 is
Phe113 N–carbonyl O, and d3 is Asn186 ND2–carbonyl O. Native Thr233 was
held out from a Kabsch fit to backbone N/CA/C/O of Phe113, Asn186, Asp272
and Asp274 against `cif/n12c.cif` chain A. The unchanged native limits were
OG1 displacement ≤0.75 Å, Nalpha displacement ≤1.25 Å, |Δψ| ≤35°, and
|Δχ1| ≤20°. The intact-PA6 check required every topology bond within the
four PA6 residues (Amber 1313–1316) to exceed 0.65 Å and remain ≤1.30 Å
when either atom was H, or ≤1.85 Å otherwise. All Tier A/B positive-φ
frames passed this bond-retention check.

An observed one-water bridge required the same topology-owned water O
within 3.5 Å of both Thr233 Nalpha and the scissile leaving amide N.
Productive leaving-leg geometry additionally required an O–H···Nleave
contact with H···N ≤2.5 Å and O–H···N angle ≥120°. An observed Nalpha–H···O
contact used the same donor–acceptor/H–acceptor/angle criteria and defined
a separate strict two-leg class when both legs passed. Absence of that
incoming contact in neutral NH2 E–S did not alone reject a leaving-leg
candidate: the future proton transferred from OG1 is not present on Nalpha
in the screened state. Distances used exact triclinic minimum images,
including the truncated-octahedral cells; fractional rounding alone was not
used as the final distance criterion. Atom names, covalent hydrogen owners,
neutral NH2 inventory and declared Asp state were read from each actual MD
topology. No water coordinates, hydrogens or protonation states were edited.

## Authoritative artifacts

- `bounded_screen_acceptance.json`: combined scientific/process outcome.
- `bounded_screen_counts.csv`: six-case count table; JSON remains authoritative.
- `primary_config.json`, `primary_attempt1/`: original primary protocol, source
  snapshot, full source/input/trajectory hashes, per-frame gates and zero candidates.
- `controls_config.json`, `controls_attempt1/`: separate four-replica controls
  with explicit per-case residue-state validation and unchanged thresholds.
- `diagnostics_attempt1/diagnostics.json`, `primary_nearpasses.csv`: all eight
  primary nearpasses and overlapping native-failure counts for all six replicas.
- `primary_job/`, `controls_job/`, `diagnostics_job/`: tmux commands,
  nonsecret environment overrides, source/input hashes, logs and exit-zero status.
- `focused_coverage.json`: nine focused tests; screen 94%, geometry 99%,
  diagnostics 96%, combined 96% statement coverage. Ruff passed.

The original primary attempt must be reproduced with its retained source
snapshot and original config. The subsequent control screen requires the
explicit per-case `residue_names` declarations in `controls_config.json`.
The analysis interpreter was
`/data/jpc/code/enzyme-substrate-modeling/.pixi/envs/analysis/bin/python`;
SciPy Amber NetCDF I/O and ParmEd were already installed. No dependencies
were added. Exact topology references and hashes are in each attempt.

## Interpretation limits

This is a negative preorganization result in the sampled fixed-charge E–S
ensembles. It does not exclude transient water reorganization, multiple-water
networks or other conformations. Frames are correlated, not independent
chemical trials. No ensemble-averaged acidity, protonation population,
proton-transfer rate or catalytic free-energy barrier is inferred. The
screen concerns substrate hydrolysis by mature Nyl12, separately from
maturation/autocleavage. Successful process completion does not turn this
negative screen into an accepted starting structure.
