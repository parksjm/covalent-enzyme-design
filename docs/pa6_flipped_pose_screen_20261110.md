# PA6 L4 flipped-pose screen — 2026-11-10

## Mechanism-based orientation requirement (2026-11-12)

Retain signed Herschlag phi as a primary mechanism-based filtering metric.
Under our fixed scissile-amide atom mapping, positive phi identifies the flipped
pose and negative phi the original pose. Raw distributed model30 has +72.24
degrees and its MM-minimized endpoint +97.948 degrees. Use positive phi for the
flipped-pose direct acylation hypothesis together with existing Tier A/B and
native-Thr gates; absolute phi alone loses this orientation information. This
is a mechanism-specific hypothesis, not proof that negative-phi poses cannot
react through other routes. Record the atom mapping and monitor sign throughout
relaxation and MD. Near-zero values require inspection.

The approved comparison is Asp272-COOH/Asp274-COO- with neutral Thr233
N-termini in all four sites. Measure all four OD1/OD2 cross-residue distances
and their minimum to avoid oxygen-label artifacts. Short O-O distances support
but do not establish a protonation assignment. Calculations remain pending.

Current frame1603 E–S projected Herschlag phi is -114.341 degrees. With the
same atom mapping/convention, positive catalytic phi identifies the flipped
local amide binding pose. Earlier statements suggesting phi could not establish
that pose were overbroad. Distances assess catalytic usefulness, not orientation
identity. Existing tier filters already use abs(phi); neither sign was excluded.

Screen executed in tmux `nyl12_orientation_screen_20261110` and completed with
exit0. Script: `scripts/screen_pa6_binding_orientation.py`. Inputs, hashes and
all measurements: `analysis/pa6_binding_orientation_20261110/`.

Scope: 82 coordinate artifacts (51 historical OpenMM-relaxed WT models0–50,
5 independent CCD OpenMM models,1 CCD ds1 model,10 locally available distributed
raw models,12 template-conditioned raw models,3 AF3-contact raw models).
Raw/relaxed versions are separate artifacts, not independent statistical
replicates. Distributed raw availability is partial. All four protein sites and
all mapped central-amide ligand copies were screened; other oligomer lengths,
mutants, capped models and incomplete heterodimer contexts were excluded.

TierA: d1<=3.2, d2/d3<=3.8 A, BD90–120 degrees, abs(phi)70–110 degrees.
TierB: d1<=3.5, d2/d3<=4.0 A, BD85–125 degrees, abs(phi)60–125 degrees.
Native Thr gate after 16-atom scaffold fit: OG1 shift<=0.75 A,N shift<=1.25 A,
abs(delta psi)<=35 degrees,abs(delta chi1)<=20 degrees. No new acceptance
threshold for terminal-N/leaving-N distance was introduced.

15 site/ligand pairs pass TierA/B;8 also pass native-Thr gates. Three joint
passes have positive phi (flipped relative to current E–S):

| Ensemble/model | Tier | phi (degrees) | Nterm–leavingN (A) | OG1–leavingN (A) | N shift (A) | OG1 shift (A) |
|---|---|---:|---:|---:|---:|---:|
| historical_openmm/pa6_l4z_l2_nyl12_wt_model_23 | A | 99.01 | 4.14 | 2.98 | 0.70 | 0.19 |
| distributed_raw_partial/pa6_l4z_l2_nyl12_wt_model_30 | B | 72.24 | 3.77 | 3.31 | 0.39 | 0.40 |
| distributed_raw_partial/pa6_l4z_l2_nyl12_wt_model_25 | B | 96.64 | 4.39 | 2.77 | 0.41 | 0.35 |

Historical model40 has positive phi and TierA but fails native-Thr displacement.
Model23 is NOT new: its previous QM/MM work is in `qmmm/pa6_model23/` and did
not establish a competent converged E–S. A flipped sign alone therefore is not
a cure. Distributed raw model30 is the best NEW geometry candidate for testing
terminal-N approach (3.77 A vs4.44 A in current E–S); model25 is a backup.
Neither is MM-minimized or MD-equilibrated by this screen, and 3.77 A alone
is not a transfer-ready hydrogen bond. Next stage, if pursued: build with the
validated RESP/GAFF2 substrate and ff19SB; unrestrained MM minimize, remeasure,
then MD only if useful geometry survives. No new MD or QM/MM jobs launched here.

NH3+ microstate hypothesis retains OG1–HG1 as proposed leaving-group donor.
See `docs/nyl12_alternative_proton_routes_20261110.md` for Asp/water comparisons.
Original NEB left running unchanged; cycle329 was observed active/unconverged.

## Authorized next stage: model30 (2026-11-10)

User requested starting one flipped model, leaving the current NEB running,
and stopping direct proton-transfer US if further sampling is unproductive.
Selected distributed raw model30, not historical relaxed model30: these are
distinct ensembles. The source path and hash are pinned in the screening
manifest. This is an opposite-phi local amide pose with native-like Thr233.

Prepare with released modular PA6 L4 multiconformer RESP/GAFF2 parameters,
ff19SB protein, TIP3P solvent and the established 150mM salt/12A solvent buffer.
The initial comparison retains the established neutral beta-subunit Thr
N-termini; this isolates orientation rather than changing orientation and
protonation simultaneously. The NH3+ microstate, with OG1-HG1 as proposed
leaving-N donor, remains a separate mechanism test, not silently assumed here.

Raw heavy-atom coordinates require explicit name mapping to modular residues
and hydrogen addition; no legacy AM1-BCC charges should enter the new topology.
Check residue links, released atom charges/types, protein chain breaks and
preservation of supplied heavy atoms. Run fully unrestrained periodic MM
minimization (ntr=0,nmropt=0,ntc=1,ntf=1), maximum10000cycles, target RMS
force0.001 kcal/mol/A. A cycle limit is not convergence. Re-measure phi,
TierA/B, native Thr displacement/torsions and proton-transfer contacts before
promoting the model to MD or QM/MM. No automatic long MD campaign is justified
solely by the raw static screen.

The direct OG1-to-terminal-N US campaign tests a different, isolated proton
transfer. Its existing analyzed PMF is uphill without a transferred-proton
minimum; extending it would refine that coordinate but would not establish the
coupled attack/leaving-group route. The useful data remain mechanistic evidence,
not a converged activation free energy. Stopping it is a resource-prioritization
decision, not proof that transient OG1-to-terminal-N transfer cannot occur during
another coupled route. Final stop inventory is recorded separately.

Replica-plan discussion: "distributed" denotes the Boltz batch scheduled
across workstations, not a different modeling method. This model30 comes from
the aster portion of that batch. Current system preparation keeps all four beta
N-termini neutral. Recommended next pilot is two independent MD replicas of
that neutral system and two of an all-four-NH3+ system, after each state passes
its own unrestrained MM screen. This is a recommendation, not a launched MD
campaign. Keep protein/substrate parameters and protocol matched, adjust ions
for the +4 protein-charge difference, use independent velocities, and assess
replica agreement before extending or adding a third. Two replicas do not prove
convergence. Substrate remains PA6 L4 zwitterion in both cases.

Actual launch: violet tmux `nyl12_model30_mm_20261110`; root
`qmmm/pa6_flipped_model30/`. Mapping smoke accepted33heavy atoms and32bonds.
Immutable source snapshot and job input hashes retained. Build/minimization
completion must be checked in job/status.json and per-stage acceptance files.
DirectPT US stopped with22ps/window completed rawdata retained; latest analyzed
prefix12ps. See docs/direct_pt_umbrella_stop_20261110.md. NEB unchanged.

## First model30 MM endpoint

Completed10000 fully unrestrained pmemd.cuda minimization cycles; exit0,
RMSforce0.042036 kcal/mol/A exceeds requested0.001. Not force-converged.
Topology checks passed104922 atoms; all4 beta NH2 retained. Parameter release
and source hashes retained. Scripted endpoint audit passes TierA and nativeThr:
phi+97.948degrees,d1/d2/d3=2.882/2.978/2.871A,BD93.082degrees,
Nterm–leavingN3.940A,OG1–leavingN3.040A. Native OG1/N shifts0.263/0.698A,
psi/chi1 differences-21.662/+5.626degrees. Flipped pose survives this relaxation.
Artifacts: `qmmm/pa6_flipped_model30/01_mm_minimization/endpoint_complex.pdb`,
`endpoint_geometry.json`, `measure_endpoint.py`, and scientific_acceptance.json.
Actual-input heavy mapping smoke and topology audits passed; Ruff passed.
New converter unit-test coverage was not measured. No MD/QMMM launched.
