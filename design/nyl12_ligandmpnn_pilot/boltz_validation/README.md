# Boltz-2 validation of the eight LigandMPNN variants

PI-approved scope, 2026-11-22: WT plus the eight deduplicated Model23 variants,
five matched prediction seeds per sequence, **four zwitterionic PA6 L4 molecules
per tetramer**. The 45 models provide 180 occupied-site observations per stage;
sites in one model are nested, not independent replicates. No activity ranking
is implied by geometry or generator frequency.

Prediction uses the existing Boltz-2 runtime and frozen source, native n12c.cif
protein template, three recycles, 200 diffusion steps and one sample per process.
WT homolog MSAs are reused identically, changing only the query row for each
variant. Two inherited generation contacts per site establish ligand placement.
Those contacts are generation aids and are absent from Amber relaxation. No
multiligand affinity calculation is requested.

## Verification snapshot — 2026-11-22

The campaign status receipt records **21/45 normalized** and **15/45 relaxed**
models; status remains `running`. This is a campaign-specific status snapshot,
not a fresh host-wide liveness audit. The exact counts and source hash are in
[`verification/checkpoint_receipt.json`](verification/checkpoint_receipt.json).

WT seed 01 completed minimization with process exit 0 and force convergence:
final RMS gradient 0.048767 against a 0.05 kcal/mol/Å target. Its separate
four-site geometry preflight reports raw 0/4 Tier A/B, 4/4 native-pass and 0/4
joint-pass; after relaxation it reports 1 B, 2/4 native-pass and 0/4 joint-pass.
This is one model seed with four nested sites, and carries no
activity or variant-ranking interpretation. Tables, report and plots are under
`checks/WT_seed_01_analysis/`.

The existing focused test run reported 55 passes plus one laptop-pull test and
88.90% total coverage (89% rounded). All ten Python modules display at least
82% when rounded to whole percentages; the raw minimum is 81.61% in
`design_component_normalization.py` (82% rounded). Coverage JSON was exported
from the saved `/tmp/nyl12_design_validation_coverage` data without rerunning
tests. See `verification/coverage.json` and the checkpoint receipt for exact
values and source hashes.

Violet and tulip run disjoint prediction lists on four GPUs. Geranium consumes
hashed, graph-checked models using its two existing Amber GPUs. PI explicitly
approved both shared hosts after workload inspection; other users' login
sessions and all earlier QM/MM/US jobs remain untouched. The preceding WT MD
campaign completed naturally before these GPU jobs started. Azalea's resident
model servers are unchanged. No environment/dependency installation was needed.

The first live WT component exposed different Boltz atom labels from the
parameterization input. Every model is therefore matched to the pinned PA6
chemical graph, including bond orders, formal charges and the central amide.
The actual WT graph has one exact 33-heavy-atom isomorphism per substrate;
C54/O49 retain their intended roles, native N57 maps to canonical N59, and
native C61 maps to canonical C62. Original predictions and component pickles
are retained. Original coordinates and protein lines are unchanged by naming
normalization. Geometry and FF transfer use these normalized names.

Relaxation uses ff19SB plus the MD-ready release
`pa6_l4_modular_resp_gaff2_v1_20261031`, including its improper corrections,
TIP3P with a 12-Angstrom solvent buffer and approximately 150 mM NaCl. All
four catalytic Thr N-termini are neutral NH2; Asp272 and Asp274 are anionic.
Other termini and standard residue states follow the matched WT baseline.
The old whole-ligand MOL2 supplies graph/naming correspondence only; its
Gasteiger values are never loaded as force-field charges. Each PA6 copy's
Amber atom types/RESP charges, three peptide links and zero formal charge are
verified, as are all neutral Thr termini and system/counterion charges.

All solute and solvent atoms move freely during minimization: ntr=0,
nmropt=0, ntc=1, ntf=1; no contact, positional or reaction-coordinate restraints.
The limit is 10,000 cycles with RMS force target 0.05 kcal/mol/Angstrom.
Process completion and force convergence have separate receipts; an endpoint
that reaches its cycle limit is retained as diagnostic, not called converged.
These are structure-relaxation tests, without new dynamics or QM/MM refinement.

Useful locations, identical on workstations and in the curated laptop pull:

- `status.json`: planned, normalized and minimized counts from the coordinator.
- `prepared/candidate_manifest.json`: sequence identities, mutations, lineage,
  all five seeds, and the 45 unique model identities.
- `predictions/canonical/VARIANT/VARIANT_model_N.pdb`: original Boltz models,
  accompanied by native confidence JSON; indices N=0..4.
- `models/VARIANT/seed_NN/raw/normalized.pdb`: input with verified canonical
  ligand atom names; NN=01..05. `mapping.json` records the original labels.
- `models/VARIANT/seed_NN/relaxation/relaxed.pdb`: canonical heavy-atom endpoint;
  `acceptance.json` records force convergence, protonation and topology checks.
  Full Amber inputs/restarts and a solute optimization trajectory stay in the
  same model's relaxation directory for reproducibility.
- `analysis/occupied_sites.csv`, `all_sites.csv`, `model_summary.csv`: generated
  after all 45 relaxations; include d1/d2/d3, Burgi-Dunitz angle, signed phi,
  native-Thr displacement/torsions, Tier A/B, and joint native/Tier criteria.
- `analysis/relaxation_acceptance.csv`: force convergence beside every model.
- `analysis/plots/geometry_before_after.png` and
  `tier_native_comparison.png`: paired seed/model comparisons across all sites.

Identical ligand labels can exchange sites. Primary occupied-site assignment
minimizes the sum of four OG1-C54 distances over all one-to-one permutations;
all 16 ligand-site pairs and the original intended assignments are preserved.
Do not equate four passes in one model with four independent successful designs.
Compare candidates to WT under this exact protocol; retain every seed in the
denominator and inspect force/geometry failures before choosing follow-up MD.

Persistent sessions: violet `nyl12_design_boltz_violet_gpu0_20261122`,
`nyl12_design_boltz_violet_gpu1_20261122`,
`nyl12_design_validation_checked_20261122`; tulip
`nyl12_design_boltz_tulip_20261122`; geranium
`nyl12_design_relax_checked_gpu0_20261122` and
`nyl12_design_relax_checked_gpu1_20261122`.

Initial setup failures remain recorded: a frozen-script import path and a
system-wide charge check tighter than the approximately 0.000017 e accumulated
rounding of four PA6 copies. The checked runtime fixes only import discovery
and that numerical gate (1e-4 e globally); local chemical checks and all charges
are unchanged. The original source snapshot, failed receipts and Geranium WT
preparation history are retained. Checked runtime imports and the real WT dry
topology passed before restart. No failed stage is silently called successful.


If WT fails to retain the geometry under this matched protonation/relaxation
protocol, do not treat the same failure in a variant as evidence of inactivity.
Review the control, protonation and preparation before setting a rejection rule.
The same frozen Boltz source is used on both prediction hosts; their existing
PyTorch/CUDA builds differ and are recorded as a potential runtime/batch factor.
No quantitative activity conclusion is made from this exploratory ensemble.


`dispatch/runtime_versions.json` records the two existing prediction runtimes;
`prediction_host_assignments.csv` maps every seed/model to its host. The Boltz-2
checkpoint SHA256 is identical on both hosts. Cross-version RDKit pickle warnings
are retained in the coordinator log; only graph/name metadata are read, and
atom count, full bond orders, charges and central-amide roles are checked after
decoding. All structural coordinates come from the untouched PDBs. For a later
quantitative paired-seed comparison, check shortlisted variants and WT together
on one common runtime; no such repeat is automatically queued here.

PA6 L4 is a tractable local substrate model. L4-to-two-L2 turnover remains an
unproven modeling hypothesis; it is not the experimentally established dominant
product route or a general model of a long polymer's terminal chemistry.
