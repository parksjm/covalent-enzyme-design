# Curated data and environment inventory

`data_manifest.json` is the authoritative file allowlist. Every selected source has a size and SHA256. `core` contains all three existing WT systems and compact analysis inputs/results. `inspection` adds self-contained QM/VMD bundles to core: choose one archive, rather than downloading both. `raw_reference_only` records twelve excluded full production trajectories (about 75.55 GB); their hashes were not invented or measured by reading the bulk files.

## Existing WT systems

| Dataset | Four-site Asp-pair state | Starting coordinates |
|---|---|---|
| `wt_asp272_ash274` | ASP272 / ASH274 | Existing unrestrained minimized restart |
| `wt_asp272_asp274` | ASP272 / ASP274 | Existing unrestrained minimized restart |
| `wt_ash272_asp274` | ASH272 / ASP274, optional scientific control | Existing unrestrained minimized restart |

Each dataset supplies a matched topology, restart, explicit protonation/build inputs, two original MD protocols, equilibrated and 100 ns endpoint restarts, protected residue identities, and pinned source acceptance evidence. No new systems or simulations were created. The normalized `md_start_receipt.json` verifies existing process completion, readable coordinates, topology and coordinate hashes. It does **not** assert force convergence or a stationary point. The doubly anionic baseline has recorded RMS gradient 0.042036 kcal/mol/Å and a passing MD-start contact gate, while its original strict gradient gate is false.

`datasets/manifest.json` uses paths relative to its own directory. The runner independently checks all 24 protected residues, all four neutral beta Thr termini, and the declared four-site ASP/ASH inventory before rendering MD inputs. A changed protonation state requires rebuilding the topology and counterions; none is substituted silently.

## Numbering and original inputs

The original tetramer has alpha chains A–D, each 232 residues, and mature beta chains E–H, each 96 residues. Native precursor beta residues 233/272/274 correspond to beta PDB residues **1/40/42**. Amber catalytic beta starts are 929/1025/1121/1217; Asp272 indices are 968/1064/1160/1256 and Asp274 indices 970/1066/1162/1258. Alpha protected residues 113/156/186 are PHE/LYS/ASN. The manifest retains both precursor and actual structure numbering.

Original model30 PDB, modular conversion map, whole-L4 source/atom map, template CIFs and original PA6 L4 WT/eight-mutant Boltz YAMLs are included. Historical mutants are references, not newly prepared or validated ASH274 mutant systems. Current build modules that were removed from the development tree are represented by their actual retained executed-source snapshots under `preparation_source`; historical absolute source paths in receipts are provenance, not portable executable defaults.

## Force field and build recipe

The authoritative release is `pa6_l4_modular_resp_gaff2_v1_20261031` from `parksjm/polyamide-amber-ff`. All release files are included, with every declared release artifact hash checked. Its exact scope is the zwitterionic L4 sequence P6N–P6I–P6I–P6C, RESP charges plus GAFF2 and the accepted improper terms, with mbondi3 radii. It is not a classical TI1/acyl-enzyme parameter set.

The actual LEaP inputs source ff19SB and TIP3P. No separate custom neutral-Thr library is missing: all four beta N termini explicitly remove H3 and use N −0.4844, H1/H2 +0.1224. These values, expected HD2 inventory, each protected Asp/ASH residue charge and all 24 identities were measured directly from each source topology. The inherited small total charge residual is recorded, never renormalized. Standard Amber libraries and executables must come from the student's Amber installation; software binaries/installers are excluded. The current PA6 repository/release root has no LICENSE file; no license text was inferred. The separately retained legacy repository LICENSE identifies its own original provenance only.

## Environments and existing results

The actual six-replica analysis command used `/data/jpc/enzyme-substrate-modeling/.pixi/envs/analysis/bin/python3.12`. This is a **different** environment from `/data/jpc/code/enzyme-substrate-modeling/.pixi/envs/analysis`, which lacks MDAnalysis. Exact observed packages and the original pixi configuration/lock are supplied. The historical environment currently contains MDAnalysis 2.10.0, NumPy 1.26.4, ParmEd 4.3.0, SciPy 1.16.3, pandas 2.3.3, matplotlib 3.10.8 and Biopython 1.86. No dependencies were installed during this inventory.

Actual source MD logs identify Amber 26 PMEMD, and the inspection export log identifies CPPTRAJ V7.6.2. Original replica protocol/extension receipts retain the executable hashes, stage seeds and commands. Existing production is two independently seeded 100 ns replicas per WT microstate, with 2 fs steps and 10 ps production frame spacing. The six full per-frame CSVs and summary JSON are available without the excluded raw trajectories. Existing matched catalytic geometry comparison is for the AE pocket; it does not claim all four pockets were independently analyzed or independently sampled.

`analysis/raw_comparison_manifest.json` rewrites the original raw-analysis paths relative to the extracted bundle root. Its companion command receipt specifies prerequisites explicitly: the twelve optional raw trajectories must first be acquired separately at those paths. Core CSV inspection needs no trajectory transfer. The exact ASH272 SI generator, its two matching 10,000-frame CSVs and PNG/PDF panel are included as **ASH272**, not relabeled ASH274. The older WT SI example is marked a historical style reference. Producing a matched new ASH274/all-four-site panel remains student work, not a completed result.

The inspection profile pins concrete immutable `update_*` generations rather than a live `latest` pointer. It includes topology, DCD, VMD script and associated tables for the reverse-coordinate path and the original NEB band. These visualizations are inspection artifacts: a provisional path/maximum is not a validated transition state, barrier or PMF.
