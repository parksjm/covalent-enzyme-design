# Matched Asp-pair MD persistence and water screen

Status: exploratory analysis, 2026-11-16. All six matched 100 ns replicas
were classified using the same Tier A geometry definitions and 10 ps frame
cadence. The source trajectories, topology hashes, 10 ns block counts and
all-frame geometry tables are recorded in
`qmmm/pa6_flipped_model30/08_asp_pair_md_comparison/analysis_attempt2/`.
The water scan used one Tier A-producing replica per microstate and is under
`qmmm/pa6_flipped_model30/09_propka_static_tier_a/water_screen_attempt1/`;
the tmux command, script hash and process log are in `water_job/`.

| Microstate | Replica | Tier A / 10,000 | Tier A in final 50 ns | Median site-E Asp O–O (Å) |
| --- | --- | ---: | ---: | ---: |
| ASP272/ASP274 | 1 | 37 | 0.0% | 4.89 |
| ASP272/ASP274 | 2 | 238 | 0.1% | 4.90 |
| ASH272/ASP274 | 1 | 3,886 | 37.1% | 2.56 |
| ASH272/ASP274 | 2 | 1,724 | 0.6% | 2.56 |
| ASP272/ASH274 | 1 | 1,677 | 12.9% | 2.73 |
| ASP272/ASH274 | 2 | 0 | 0.0% | 2.65 |

The protonated-Asp systems tend to pull the two Asp oxygens together and
have more Tier A frames in some replicas, but late Tier A persistence is
not reproducible across either protonated pair. The frames are correlated;
the six replicas, rather than 60,000 independent frames, set the strength
of this comparison. Fixed-charge MD cannot transfer the Asp proton and
does not determine its preferred location.

For the water screen, a **candidate** requires one water oxygen within
3.5 Å of both Thr233 OG1 and the named nitrogen. This is a proximity
criterion only, without hydrogen-bond angles or proton orientation.

| Replica | OG1/Thr-N candidate frames | Among Tier A frames | OG1/leaving-N candidate frames | Among Tier A frames | Longest same-water streak, OG1/Thr-N |
| --- | ---: | ---: | ---: | ---: | ---: |
| ASP272/ASP274 replica 2 | 7 / 10,000 | 4 / 238 | 996 / 10,000 | 3 / 238 | 1 frame |
| ASH272/ASP274 replica 1 | 1,804 / 10,000 | 500 / 3,886 | 262 / 10,000 | 101 / 3,886 | 19 frames |
| ASP272/ASH274 replica 1 | 5,500 / 10,000 | 754 / 1,677 | 96 / 10,000 | 35 / 1,677 | 32 frames |

At 10 ps spacing, the longest observed same-identity OG1/Thr-N candidate
streak spans at most roughly 0.3 ns. Many oxygen-proximity events may be
water exchange rather than one ordered relay. At the selected
ASP272/ASH274 Tier A frame 3492, water residue 1576 lies 2.91 Å from
OG1 and 2.97 Å from Thr N. OG1-HG1···O(water) has H···O 1.99 Å and
angle 160.2°, but neither water H is aligned to donate to Thr N
(best H···N 2.86 Å, O-H···N 87.5°). That snapshot does **not** show a
complete two-leg proton relay despite oxygen proximity.

The explicit donor–hydrogen–acceptor check below was therefore required.
Do not infer a pseudo-ring mechanism from an oxygen distance alone. The
separate protein-only pKa pilot is in
`docs/propka_asp_pair_pilot_20261116.md`.

## Two-leg hydrogen-bond check

The follow-up scripted scan is in
`qmmm/pa6_flipped_model30/09_propka_static_tier_a/hbond_screen_attempt1/`
with its tmux job in `hbond_job/`. It required **both**
Thr233 OG1–HG1···O(water) and O(water)–H···N to pass donor–acceptor
≤3.5 Å, hydrogen–acceptor ≤2.5 Å, and D–H···A ≥120°. The N acceptor
was either neutral Thr233 N or the PA6 scissile amide N. The former
assumes the prepared neutral N-terminus; the latter is a geometry screen,
not proof that an unbroken amide N accepts a proton in E–S.

| Replica | Complete OG1→water→Thr-N frames | Of those, Tier A | Complete OG1→water→leaving-N frames | Of those, Tier A | Longest same-water complete Thr-N streak |
| --- | ---: | ---: | ---: | ---: | ---: |
| ASP272/ASP274 replica 2 | 1 | 1 | 1 | 0 | 1 frame |
| ASH272/ASP274 replica 1 | 25 | 6 | 5 | 0 | 1 frame |
| ASP272/ASH274 replica 1 | 367 | 32 | 0 | 0 | 3 frames |

The 32 Tier A frames with a complete Thr-N bridge in the ASH274 replica
are 1.9% of its Tier A frames, and no same water retains both hydrogen
bonds for more than three consecutive 10 ps samples. No Tier A frame
has a complete two-leg bridge to the leaving-group N in these three
replicas. This is weak support for an *ordered, persistent* water relay
under the current neutral-Thr setup, though short-lived reactive
arrangements cannot be ruled out. The priority is to verify substrate
register and local protonation for one representative per microstate;
do not select a pKa or QM/MM starting state by water-oxygen proximity
alone.
