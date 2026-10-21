# QM method and substrate checkpoint — 2026-10-20

## Agreed calculation protocol

Use Gaussian 16 `wB97XD/def2TZVP` for the first neutral free-NH2
L-threonine methylamide model optimizations and harmonic frequencies.
No additional empirical dispersion correction is applied. Use UltraFine
integration, tight optimization/SCF, explicit charge 0 and multiplicity 1.
The model is H2N–CH[CH(OH)CH3]–C(=O)–NHCH3, C5H12N2O2,
with CA S / CB R stereochemistry. Two deterministic ETKDG starts have initial
psi values of -40 and +140 degrees; the Gaussian optimizations are unconstrained.
These two starts are an initial sample, not an exhaustive conformational ensemble.

The separate Link1 calculation uses HF/6-31G(d) Merz–Kollman ESP with
`Pop=MK IOp(6/33=2,6/42=6)`. Geometry optimization is not performed at HF.
The runner checks optimization completion, all 57 vibrational modes positive,
normal termination, atom inventory/order, stereochemistry, conservative bond
and short-contact screens, electronic state, and ESP extraction. Only accepted
QM artifacts enter the content-addressed cache; these are not fitted parameters.

The literature precedent is Robertson, Tirado-Rives, and Jorgensen,
[JCTC 2015, DOI 10.1021/acs.jctc.5b00356](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00356).
Their OPLS-AA/M development used relaxed peptide scans at
omegaB97X-D/6-311++G(d,p), with B2PLYP-D3BJ/aug-cc-pVTZ single-point energies
for the principal fitting reference. Thus the functional has direct peptide
parameterization precedent, but our def2-TZVP choice and neutral terminal model
are not a reproduction of that protocol. The
[QUBE protein-force-field study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6740169/)
also used those relaxed peptide scans. This supports the geometry choice;
it does not establish ff19SB compatibility by itself.

RESP fitting, cap constraints/charge transfer, bonded terms and torsional
validation remain to be implemented and reviewed. Retain AMBER conventions
for protein atom types, Lennard-Jones terms and 1–4 scaling; do not import OPLS
parameters merely because its QM geometry protocol is useful. Validate modified
terminal and covalent-state parameters against QM conformations/energies and
the intended ff19SB environment before calling them production-ready.
No standard protein residues are being reparameterized by these two jobs.

DFTB3/3ob remains the proposed exploratory method for local TI1 optimization
and reverse reaction paths. It must be checked against higher-level reference
calculations for the oxyanion and proton-transfer chemistry before interpreting
barriers quantitatively. The new Gaussian jobs optimize small parameter models,
not enzyme coordinates or TI1.

## Actual jobs and retained evidence

Output root:
`/data/jpc/covalent-enzyme-design/output_nyl12_resume_20261020/neutral_thr_nme_wb97xd/`.
Input preparation finished successfully in tmux. Both Gaussian jobs started on
violet and produced `E(RwB97XD)` SCF output. Sessions:

- `nyl12_thr_wb97xd_01_20261020`
- `nyl12_thr_wb97xd_02_20261020`

Each has four shared cores and 4 GB memory. Exact commands and terminal status
are under `jobs/conformer_01/` and `jobs/conformer_02/`; live Gaussian output is
under `calculations/conformer_01/gaussian.log` and its conformer_02 counterpart.
Sources and lockfile are snapshotted under `source/`; preparation inputs are
retained under `prepared/`. Each calculation has `calculation.json` recording
input/tool hashes, versions and its separate scratch directory. Scratch is under
the authorized `/data/scratch/gaussian_scratch/jpc/`. Cache lookup occurs before
execution. A running wrapper's null exit status is not convergence or acceptance.

## Substrate identity and confidence

The polymer family label is poly(X,AMA). The modeled eight-monad sequence is
**MAXAMAXA**, with scission at **MAXA|MAXA** (A4–M5). Existing target identifier
`poly_xama_maxa` is retained for provenance. The graph yields two identical MAXA
hydrolysis fragments, C26H42N4O5, exact mass 490.315520444 Da. No substrate atoms,
products, or coordinates were changed by the audit.

The successful tmux report is
`output_nyl12_resume_20261020/confidence_audit_v2/report.json`.
It retains atom-level pLDDT, monad assignments, directional interchain PAE,
original ranking values, and source/input hashes. The initial audit failed
because the summary ranking score is rounded to two decimals while the CSV
retains full precision; that attempt remains under `confidence_audit/`.
The corrected implementation checks explicitly for two-decimal agreement.

| Model | Ligand mean pLDDT | Atoms below 50 | Reactive C / O / N pLDDT |
|---|---:|---:|---|
| MAXA\|MAXA | 29.6412 | 59/69 | 83.84 / 54.43 / 60.16 |
| PA6 L2 | 60.8382 | 11/33 | 92.53 / 85.86 / 87.64 |

Terminal monads 1 and 8 average 13.3582; internal monads 2–7 average 34.9644.
Uncertainty extends beyond the ends. The high whole-complex ranking score must
not be used to certify the ligand pose. Retain full MAXA|MAXA as a reference
consistent with the intended products, explore alternative tail conformations,
and use PA6 as the stronger initial local-geometry case. Low confidence alone
does not demonstrate excessive substrate length. Any truncation must be a
separately identified sensitivity model with explicit caps/charge and product
mapping; no truncation is currently implemented or silently substituted.

## TI1 to ES and restraint effects

The intended route is to repair and validate TI1 first, then trace backwards
towards ES. Reverse both the Thr OG1–substrate C bond formation and the proton
transfer: the intended ES has neutral terminal NH2 and Thr OG1H, with the
substrate carbonyl restored to trigonal-planar geometry. Removing a bond alone
does not produce the correct ES state. Preserve atom mapping and record d1,
d2, d3, a_BD and the existing projected/signed Herschlag phi convention along
the path. TI1 and ES have different geometric targets; do not constrain ES to
the TI1 bond length or tetrahedral carbonyl geometry.

The previous backbone positional restraint (5.0 kcal/mol/angstrom^2 in the
Amber restraint term) can retain strain from the starting pose and obstruct
local rearrangement. After repairing initial clashes, compare a documented
restraint taper and a less-restrained endpoint, with the same input chemistry,
monitoring active-site geometry, contacts and proton/connectivity inventories.
These are planned sensitivity checks, not evidence that a relaxed state exists.

Both canonical starting TI1 models still fail the water-contact gate described
in the standalone handoff. No new Sander calculation or accepted TI1 minimum
has been produced. Local conformer/geometry repair, calibrated reaction paths,
separate ES topology construction and validated force-field fitting remain.

## Design and multistate screening clarification

The user's follow-up asks whether TS1 guides RFD3/LigandMPNN, whether redesigned
sequences receive covalent AF3 models, and whether independent ES predictions
and later catalytic states should be included. Recommended staged workflow:

1. Establish a reference TI1 and reverse path to ES. A scan maximum is only a
   TS candidate: refine the saddle, verify the intended unstable mode and
   endpoint connectivity (IRC or equivalent downhill path). Do not assume TI1
   exists as a metastable minimum or that attack and proton transfer comprise
   a single saddle. Map the reference enzyme's remaining acylation and
   deacylation paths coarsely before committing a campaign to TS1 alone.
2. Extract a small TS1 atomic motif/ensemble for geometry-conditioned design,
   retaining catalytic Thr, oxyanion donors and proton-relay requirements.
   For Nyl12, start with the existing symmetry-tied focused LigandMPNN mask;
   use RFD3 local backbone/loop diversification if sequence redesign cannot
   supply suitable geometry. Preserve maturation residues, processed chain
   architecture and oligomer interfaces. RFD3 can condition on catalytic atom
   coordinates; LigandMPNN proposes sequences in their atomic context. Neither
   computes transition-state stabilization or a reaction barrier. A TS is
   not encoded as a stable chemical adduct simply by inventing a full bond.
3. Predict redesigned sequences independently in both covalent TI1 and
   noncovalent ES forms using multiple seeds and matched WT controls. ES has
   its own planar-carbonyl graph and proton inventory, no enzyme–substrate
   bond, and no externally imposed active-site contact restraints. Avoid
   supplying the designed pocket pose as a template in the independent arm.
   The reverse-path-derived ES answers connectivity to TI1; fresh ES
   predictions test pose recovery independently. Both are needed. Rebuild
   hydrogens/protonation and refine covalent AF3 outputs physically; AF3 does
   not validate a saddle or its chemistry.
4. For candidates passing fold/assembly checks, screen ES, TI1, acyl enzyme
   with a hydrolytic water, and a deacylation TI2-compatible local geometry.
   Check the scissile register, state-specific geometry, clashes, oxyanion
   interactions, departing-group protonation route, water access, and room
   for products to exit. Apo checks detect pocket/fold collapse. Do not demand
   identical tail conformations or identical d1/angles across states. A single
   failed AF3 seed is not a reliable rejection, especially for MAXA|MAXA.
5. On a small diverse shortlist, recalculate the relevant acylation and
   deacylation barriers for each sequence, referenced to its corresponding
   reactant basin with balanced atoms/protons. Include TI1 collapse and
   acyl-enzyme hydrolysis, not only formation of TI1. Compare changes in
   barriers, not absolute energies across different sequences or unbalanced
   states. Screening energy paths are not free-energy/rate predictions.
   Deeper binding/intermediate wells can raise outgoing barriers if their
   transition states do not improve correspondingly. Product-release checks
   begin geometrically; investigate energetics if a product trap is indicated.

This adds multistate validation around a simple initial design objective;
it does not require multistate generative-model training or full free-energy
profiles for every proposed sequence. Unrestrained pose recovery and subsequent
restraint-free stability are supporting structural evidence, not proof of
binding affinity, productive association kinetics, or turnover. Calibrate the
filters against WT and available experimentally characterized variants.

Relevant primary/official sources:

- [RFD3 enzyme-design tutorial](https://rosettacommons.github.io/foundry/models/rfd3/tutorials/enzyme_design_tutorial.html)
  for catalytic motif conditioning.
- [LigandMPNN, Nature Methods 2025](https://www.nature.com/articles/s41592-025-02626-1)
  for atomic-context sequence design and its geometry-based features.
- [AF3 input specification](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)
  for custom chemistry and explicit bonded atom pairs.
- [Proteasome QM/MM reaction-path study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3851579/)
  as evidence that an N-terminal Thr hydrolase can have a later acylation
  bottleneck; its exact mechanism/rate-limiting step is not assigned to Nyl12.

Implementation status: this section records the recommended design/screening
plan, not implemented TS searches, RFD3 runs, candidate AF3 predictions or
multistate energy results.
