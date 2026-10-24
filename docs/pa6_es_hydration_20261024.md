# Hydration of the occupied WT PA6 ES pocket — 2026-10-24

All 6,000 saved frames from the three existing 100 ns WT PA6 L4 trajectories
were analyzed without new MD. No water oxygen was detected within 4 Å of
Thr929 N or OG1 in any saved frame. No water simultaneously contacted both
atoms within 3.5 Å. The nine selected near-reactive seeds agree with this result.

| Replica | Saved frames | Any water O within 4 Å | Same O within 3.5 Å of both |
| --- | --- | --- | --- |
| 01 | 2000 | 0 | 0 |
| 02 | 2000 | 0 | 0 |
| 03 | 2000 | 0 | 0 |

CPPTRAJ watershell used explicit `:WAT@O` solvent selection with periodic
imaging. Separate 3.5 Å shells around N, OG1 and their union give the shared
water count by inclusion–exclusion. A separate 4 Å union shell measures local
availability. All four actions returned 2,000 records per replica; counts were
internally consistent. Source SHA256 values were inherited from the completed
reuse audit, with unchanged file sizes/timestamps verified for this analysis.
The source trajectories retain their 27,690 explicit, hydrogenated waters.

This is evidence about the existing bound-start, fixed-protonation **ES model**,
not proof that hydrolysis lacks water or that a relay is impossible. The output
cadence is 50 ps; shorter visits can be missed. Buried hydration can be slow,
and reaction/protonation changes may alter access. Frames are correlated and
are not 6,000 independent experiments. A classical trajectory cannot transfer
protons or sample chemical bond changes.

For the next pilot, prioritize the observed Thr/amine/carboxylate environment
and oxyanion donors; do not force a bridging water into the ES pocket. Treat a
water-assisted alternative as a separately prepared, physically justified
hypothesis. Acyl-enzyme hydrolysis still requires attacking water in its own
state; these ES data do not determine that water's placement or entry mechanism.

Evidence: `output_nyl12_pa6_hydration_20261024/hydration_summary.json`,
`scientific_acceptance.json`, retained sources and per-replica CPPTRAJ logs.
The nine portable inspection pockets are in
`output_nyl12_pa6_es_audit_20261024/seed_inspection_attempt4/`.
Command reference: https://amberhub.chpc.utah.edu/watershell/ .
