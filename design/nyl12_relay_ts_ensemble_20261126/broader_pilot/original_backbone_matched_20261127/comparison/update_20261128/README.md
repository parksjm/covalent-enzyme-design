# Matched WT / relay 03 / relay 06 comparison

Sites ordered AE / BF / CG / DH in the canonical derived view. Failed-assembly before values are raw diagnostics; after values are unavailable. Accepted-case before values are the serialized Amber inputs. Signed positive-phi A/B assignment controls acceptance; nearest-d1 diagnostics can fail it.

| Case | Tiers before → after (all four) | +phi A/B before → after | RMS gradient | <1 Å interchain pairs | Native protein ipTM | Group protein-aligned ipSAE2 I/J/K/L |
|---|---|---|---|---|---|---|
| WT_seed_20261126 | A/A/B/B → A/A/A/A | 4 → 4 | 0.046729 | 0 | 0.4815 | 0.026397/0.016310/0.014456/0.018853 |
| WT_seed_20261127 | reject/reject/reject/reject → NA/NA/NA/NA | 0 → not run: invalid assembly | None | 591 | 0.3479 | 0.019571/0.014993/0.010802/0.016093 |
| relay_prot101_03_seed_20261126 | B/A/A/A → C/B/B/B | 4 → 3 | 0.049565 | 0 | 0.4634 | 0.016422/0.012646/0.009494/0.008297 |
| relay_prot101_03_seed_20261127 | reject/reject/reject/reject → reject/reject/reject/reject | 0 → 0 | 0.045007 | 0 | 0.4849 | 0.014738/0.011820/0.011013/0.015285 |
| relay_prot101_06_seed_20261126 | reject/reject/reject/reject → NA/NA/NA/NA | 0 → not run: invalid assembly | None | 535 | 0.3422 | 0.017145/0.011036/0.009317/0.013247 |
| relay_prot101_06_seed_20261127 | reject/reject/reject/reject → NA/NA/NA/NA | 0 → not run: invalid assembly | None | 641 | 0.3531 | 0.013981/0.011739/0.009676/0.011800 |

All three launched Amber calculations meet the saved force criterion. WT seed 20261126 retains 4/4 positive-phi A/B sites. Relay 03 seed 20261126 retains 3/4 (AE becomes C); seed 20261127 has 0/4 before and after. WT seed 20261127 and both relay 06 seeds fail assembly precheck and have no MM endpoint.

Group scores use the exact A+B+C+D+E+F+G+H rows and explicit protein-aligned ipsae12 direction. Native ligand identifiers persist through canonicalization, which relabels beta chains only. WT26 site ligands are I/J/K/L; 03-26 I/J/L/K; 03-27 I/J/K/L. Failed assemblies have repeated nearest ligands, so no four-site confidence join is asserted. Per-native-ligand scores and both interface directions are retained separately.

Detailed d1/d2/d3 (Å), alpha a_bd and Herschlag phi (degrees), selected carbonyl identity and diagnostic/accepted selection are in site_geometry.csv. Native raw 48-way measurements per case are in raw_all_amide_geometry.csv. Full input hashes and mapping evidence are in results.json. No predictions, minimizations or scientific jobs were rerun.
