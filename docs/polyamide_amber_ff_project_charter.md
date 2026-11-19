# Polyamide Amber force-field project charter

Status: proposed separate project, 2026-10-28. This document does not certify
parameters or authorize a new repository by itself.

## Purpose

Develop a reproducible, enzyme-ready Amber parameter library for PA6, PA66 and
selected covalent states in polyamide hydrolysis. The Nyl12 redesign project
will consume versioned releases from this library rather than own the general
parameterization machinery.

The working repository name is `polyamide-amber-ff`. A separate repository is
recommended because the reusable chemical library, QM training data, mapping
software, validation systems and publication have a broader scope and release
cycle than Nyl12 redesign.

## Novelty requirement

The publication cannot be based only on providing another nylon force field.
Prior work includes the ab initio-derived MSXX nylon force field, PA6 models
based on AMBER99SB charge variants, reactive INTERFACE-force-field PA6 models,
and GAFF-based polyamide simulations. A 2025 workflow review also emphasizes
the lack of standardized and reproducible polymer parameter benchmarks.

The differentiating contribution should be the combination of:

- modular, length-independent Amber residue templates for PA6 and PA66;
- charged, neutral and capped terminal states;
- enzyme-reaction-state extensions with explicit chemical definitions;
- graph-based conversion of Boltz/AlphaFold 3 ligand coordinates to Amber;
- reproducible QM, electrostatic, torsional, condensed-phase and enzyme-bound
  validation;
- versioned, tested templates and complete provenance.

Relevant starting references include:

1. Dasgupta et al., JACS 1996, DOI `10.1021/ja944125d` (MSXX nylons).
2. Lukasheva et al., Soft Matter 2017, DOI `10.1039/C6SM02169G` (PA6 with
   AMBER99SB-derived charge models).
3. Multiscale PA6 modeling with reactive INTERFACE FF, Ind. Eng. Chem. Res.
   2021, DOI `10.1021/acs.iecr.1c02440`.
4. Turney and Matta, Macromolecules 2025, DOI
   `10.1021/acs.macromol.5c01166` (polymer-building/parameterization review).

## Layered deliverables

### 1. General polyamide library

- PA6 internal repeat, charged/neutral termini and neutral caps.
- PA66 diamine- and diacid-derived repeat representation, terminal variants
  and neutral caps.
- Amber OFF/LIB, FRCMOD and canonical structure files.
- Builders for declared oligomer length, sequence and terminal chemistry.
- Atom/residue mapping from SMILES, SDF, PDB and mmCIF.

The detailed PA6 residue boundary and charge strategy are in
`docs/modular_polyamide_residue_parameterization_plan.md`.

### 2. Enzyme-state extensions

Each chemical minimum receives a separate topology and parameter set. The
library does not claim that one classical topology models chemical reactions.

| State | Required connectivity | Current readiness |
|---|---|---|
| E--S | noncovalent enzyme plus ordinary modular substrate | can start after core residue charges pass |
| TI1 | Thr OG1--scissile carbon bond; tetrahedral carbon and oxyanion | blocked on validated protonation/charge minimum |
| acyl enzyme | Thr O--acyl carbon ester; scissile C--N bond broken | blocked on QM/MM path and leaving-group protonation |
| TI2 | water-addition tetrahedral intermediate during deacylation | optional initially, required for a complete catalytic-state library |
| products | two oligomer fragments with newly formed termini | requires declared pH/protonation state |

The reaction-state layer should be mechanism-family-specific, initially an
N-terminal-threonine hydrolase extension. It must not be presented as a
universal polyamide force field.

## Chemical gates before TI1/acyl fitting

For each state, freeze the following in a machine-readable specification before
QM charge or bonded-term work:

1. scissile repeat indices and retained product fragments;
2. complete bond-order/connectivity graph;
3. proton owner for catalytic Thr N/OG1, oxyanion, leaving-group N, catalytic
   Asp residues and water;
4. total formal charge and multiplicity of the capped QM reference;
5. enzyme/substrate residue boundaries and cross-residue bonds;
6. source stationary-point coordinates and QM/MM method;
7. integer-charge constraints for reusable residues.

For the proposed PA6 L4-to-L2 reaction, the working scissile bond is between
repeat 2 carbonyl C and repeat 3 N. In the acyl enzyme, the carbonyl-containing
L2 fragment is expected to remain attached to Thr, while the N-containing L2
fragment is the leaving product. This remains a hypothesis until the QM/MM
path and proton transfers validate it.

## Parameterization and validation

- Geometry references: validated DFT or DFT/MM stationary points.
- Charges: multiconformer RESP with gas-phase HF/6-31G(d) MK ESPs unless a
  predeclared benchmark establishes a better compatible method.
- Bonded terms: GAFF2 transfer where it passes; targeted refitting against
  bidirectional DFT scans where it fails.
- Protein: ff19SB outside explicitly declared state-specific catalytic
  residues.
- Solvent: a declared Amber-compatible water/ion model.
- Validation: QM conformer/torsion energies, ESP residuals, oligomer-length
  transfer, hydration/condensed-phase properties where data exist, explicit
  water MD, enzyme-bound stability and mapping round trips.

The current full-L4 RESP/GAFF2 model and active omegaB97X-D scans are reference
data for the separate project; they are not silently relabeled as modular
production parameters.

## Proposed repository layout

```text
polyamide-amber-ff/
  chemistry/          # declarative PA6/PA66 and state graphs
  qm/                 # input generators and retained reference manifests
  parameters/         # staged candidates; only releases are certified
  amber/              # OFF/LIB/FRCMOD and assembly examples
  mapping/            # Boltz/AF3 graph mapping and canonical naming
  validation/         # QM, MD, condensed-phase and enzyme-bound gates
  examples/           # L2/L4/L6, PA66 and catalytic-state systems
  docs/               # methods, provenance and publication material
  tests/
```

Large QM/MD outputs should remain external artifacts with hashes and manifests;
the repository should retain inputs, code, summaries and selected reference
structures.

## Existing assets and exclusions

- Historical PA6 L4 zwitterion-to-L2 parameters are under
  `/data/jpc/enzyme-substrate-modeling/output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/1_relaxation/amber_params_pa6_l4z_l2/`.
  They are a 79-atom AM1-BCC/GAFF baseline, not a modular production library.
- The current 79-atom multiconformer RESP/GAFF2 PA6 L4 candidate is under
  `/data/jpc/covalent-enzyme-design/parameters/pa6_l4_zwitterion_resp_gaff2/03_multiconformer_resp/results_attempt3/`.
  Its charge sum is -0.000004 e. It is the preferred whole-L4 comparator but
  remains uncertified pending the active bonded-term benchmark.
- The corrected all-carbon PA66 L3 zwitterion-to-L1 parameters are under
  `/data/jpc/enzyme-substrate-modeling/output_nyl12_pa66_md_ligand_topology_repair_all_carbon_20260909/amber_params_pa66_l3z_l1_all_carbon_zwitterion/`.
  They are a 117-atom AM1-BCC/GAFF baseline used in the corrected PA66 MD
  campaign. They are not modular RESP/GAFF2 parameters.
- Older PA66 artifacts with OpenBabel-inferred alkyne-like bonds are invalid
  and must not be migrated or used as training references.
- The TI1 archive
  `/data/jpc/enzyme-substrate-modeling/output_nyl12_ti1_parameters_20261022/pa6_local_bonded/`
  contains capped QM and curvature references but no accepted native TI1 RESP
  fit, MOL2, FRCMOD, LIB or validated topology. Current TI1 bonded terms remain
  placeholders.
- No acyl-enzyme parameter artifact was found. This state must be built only
  after its stationary-point chemistry is defined.

## Compute allocation

- **violet:** finish the active PA6 torsion/improper scans and their analysis.
- **azalea:** currently reachable with 56 CPU cores, 503 GiB RAM, two idle
  98-GB RTX PRO 6000 GPUs, Gaussian 16 group access and AmberTools 26 Sander.
  It is suitable for parallel future Gaussian references and GPU modeling/MD
  after the exact GPU Amber/OpenMM environment is verified.
- **tulip:** reported free, but SSH from violet currently fails public-key
  authentication. Do not schedule work there until access is repaired.

Do not migrate live Gaussian scans between hosts. Use azalea for independent
new references after chemistry specifications and input hashes are frozen.

## Near-term work while current scans run

1. Inventory existing PA6/PA66 topologies, charges and MD data.
2. Freeze the PA6 residue/cap graph schema and canonical atom names.
3. Define PA66 monomer boundaries and termini explicitly.
4. Build the E--S modular state first.
5. Use the reaction-path project to supply validated TI1 and acyl-enzyme
   stationary points and protonation states.
6. Benchmark modular assemblies against whole-oligomer QM before production MD.

Repository creation and any independent publication branch should begin only
after the user confirms the scope, working name and private GitHub destination.
