# Optional design and interaction resources — 2026-11-20

These are resources for expert recommendations, **not a required software
stack or a new execution queue**. The PI clarified that we need not use them all.
The current mechanism/US/pKa work continues unchanged. Installation, model
download and production use are separate decisions.

## Resource inventory

| Resource | Intended optional use | Verified local state and readiness |
| --- | --- | --- |
| [RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina), [project](https://rfd4-proteina.org/) | Preferred candidate over RFD3 when backbone or joint sequence/structure redesign is justified; motif conditioning, ligand context and partial diffusion | Not found in bounded violet/azalea checkout/environment checks; no Nyl12 qualification yet |
| [Foundry / RFdiffusion3](https://github.com/RosettaCommons/foundry) | Retain the existing implementation and checkpoints for reproducibility and optional comparison | Violet checkout `/data/jpc/code/foundry`; rc-foundry 0.1.11.dev34+gcee116dc5. RFD3/RF3/MPNN checkpoint files exist in the default cache on both hosts; inference readiness was not retested |
| [tmol](https://github.com/RosettaCommons/tmol) | Standalone PyTorch scoring, packing, minimization/FastRelax and possible future design interfaces | Not found in bounded checks; adapter and PA6/covalent-state compatibility unvalidated |
| [ProLIF](https://github.com/chemosim-lab/ProLIF) | Preferred candidate for trajectory interaction fingerprints and residue-contact occupancies | Not found in bounded checks; native Amber/PA6 conversion and agent adapter still to qualify |
| [PLIP](https://github.com/pharmai/plip) | Complementary selected-structure interaction reports and an independent contact-definition comparison | Not found in bounded checks; local preparation and adapter still to qualify |
| [polyamide-amber-ff](https://github.com/parksjm/polyamide-amber-ff) | Reuse versioned Amber-compatible PA6 RESP/GAFF2 libraries in structural preparation, MD and comparisons | Existing `/data/jpc/polyamide-amber-ff`; scoped zwitterionic L4 release described below |
| [AF3 contact restraints](https://github.com/parksjm/alphafold3-contact-restraints) | Existing optional contact-guided E-S/adduct model generation, alongside Boltz and unrestrained prediction controls | Checkout `/data/jpc/code/alphafold3-contact-restraints` on both hosts; identified azalea environment has alphafold3 3.0.1. Local checkout and current GitHub source are separate revisions |
| [Biohub ESM family](https://github.com/Biohub/esm) | Optional ESMC sequence priors/embeddings, ESM3 generation and ESMFold2 ligand-aware structure comparison | Source-reviewed options; no new installation or inference. Pin exact released weights and artifact-specific license before use |
| [ESM-IF1](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding) | Protein-backbone inverse-folding comparison | No default explicit ligand atom conditioning; retain ligand-aware LigandMPNN as baseline |

“Not found” refers to the targeted directories and identified environments,
not an exhaustive filesystem search. Package/checkpoint presence does not prove
inference readiness. An installed research tool is not automatically an enabled
mechanism-agent tool. Both hosts report x86_64; native libraries, GPU architecture,
CUDA and driver compatibility still require host-specific qualification.

## RFD4 and tmol disposition

The PI's PLM questions and the known-active calibration plan are documented in
[`plm_and_filter_calibration_20261121.md`](plm_and_filter_calibration_20261121.md).
ESMFold2 is a current released, ligand-capable all-atom predictor; it should not
be confused with historical ESMFold. These options supplement the current
workflow only if they add discrimination on WT/known-mutant controls.

RFD4's official README describes joint sequence/structure generation, catalytic
motif conditioning, small-molecule pockets and partial diffusion. These are
relevant capabilities, not evidence that the combined Nyl12 configuration works.
Evaluate RFD4 first if we need backbone remodeling. Keep the inexpensive
fixed-length, symmetry-tied LigandMPNN baseline; RFD4 can co-design sequences,
so LigandMPNN after RFD4 is an optional comparison rather than a mandatory step.

Before adopting it, demonstrate preservation of the mature alpha/beta assembly,
four equivalent active sites, catalytic anchors, chain/residue correspondence,
PA6 atom mapping, state-specific protonation/connectivity, and accessible loops.
Native symmetry enforcement for this processed assembly has not been verified.
Use exact WT lengths initially; indels require the explicit variant mapping
contract in `redesign_backbone_symmetry_and_access_20261022.md`. A nonstationary
TS geometry must not silently become an ordinary stable ligand bond graph.

RFD4 code and the official model card identify BSD-3-Clause licensing; optional
third-party assets retain separate terms. Pin the paired flow/autoencoder
revision and runtime before a pilot. The checked RFD4 pyproject requires
Python >=3.12 and torch >=2.11, whereas the current Foundry environment has
torch 2.7.1: do not upgrade that working environment in place.
RFD4's checked CUDA extras constrain tmol to >=0.2,<0.3, while the standalone
tmol README describes 1.0. Keep their environments/version contracts explicit.
tmol is a possible standalone design component, not an automatic substitute for
Amber relaxation or QM/MM. Amber OFF/frcmod files are not thereby tmol parameters;
nonstandard PA6 and reactive-state scoring need separate qualification. Neither
an interaction score nor a tmol score is a chemical activation free energy.

## Interaction-analysis policy

ProLIF explicitly supports interaction fingerprints from MD, docking and
experimental complexes. PLIP profiles noncovalent contacts in selected structures.
The proposed division is ProLIF for ensemble occupancies and PLIP for selected
static cross-checks; one tool may suffice. Different contact definitions need not
give identical counts. Preserve explicit hydrogens/protonation, formal charges,
bond orders, waters, atom mapping and all four site identities. PLIP documents
`--nohydro` for already hydrogenated inputs; do not let a convenience preparation
step silently change the microstate. Qualify conversion of modular PA6 residues
and covalent adducts before interpreting fingerprints.

Track hydrogen bonds, charged contacts, substrate retention, oxyanion-hole
contacts and candidate water networks where supported by the selected tool.
Keep the validated d1/d2/d3, signed Herschlag phi, Burgi–Dunitz and native-Thr
metrics as separate catalytic-geometry evidence. Fixed-topology contacts alone
do not establish proton-transfer mechanisms, affinities or catalytic activity.

## Existing PA6 and AF3 assets

Use the authoritative manifest for
`pa6_l4_modular_resp_gaff2_v1_20261031`:
`/data/jpc/polyamide-amber-ff/parameters/releases/pa6_l4_modular_resp_gaff2_v1_20261031/manifest.json`.
It covers the zwitterionic, net-neutral `P6N-P6I-P6I-P6C` L4 model with RESP,
GAFF2 and three accepted improper overrides. It is not blanket validation of
L2, arbitrary polymer lengths, bulk nylon, TI1 or acylenzyme parameters.
The force-field project records remaining compact-conformer transferability
checks. Compare modular RESP, whole-L4 RESP and historical AM1-BCC on identical
independent grids before deciding whether a targeted refit is warranted. L2/L4 accuracy acceptance
and final-results manuscript integration remain incomplete. Reuse this work;
do not start an independent parameterization here.

The public AF3 fork documents noncovalent diffusion-time contact guidance and
separate validation. The audited local checkout has an older, upstream-style
README and a workstation remote; do not assume it matches current GitHub main.
Preserve both identities and reconcile them before another run. Do not equate
guided contacts with independent binding evidence. Continue to evaluate
unrestrained E-S recovery, relaxation and state-specific geometry. AF3 source,
weights and output terms remain separate prerequisites.

## Evidence and expert use

The source URLs, exact Git blob identities/content hashes, local checkout
revisions, package versions, audit scope and local status-file hashes are in
`../provenance/redesign_tool_inventory_20261120.json`.
The public PA6 GitHub README was not accessible through the connected GitHub
tool during this audit; release scope comes from the current local manifest
and November 20 checkpoint, not an assumed remote state.

The expert catalog is maintained at
`/data/jpc/mechanism-agent/docs/TOOL_CATALOG.md`. Supply the relevant subset
with source IDs in future review packets; these options were not retroactively
added to completed reviews. Experts should explain the scientific question a
tool would resolve and prefer validated existing resources. A catalog entry
does not authorize installation or model-issued execution.
