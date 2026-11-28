# Early projected-PMF diagnosis: exact displayed 70 ps + 43-history pool

This diagnostic holds all **82,000 samples and existing MBAR weights fixed**.
It does not refit, smooth, alter bins, change sampling or replace a publication.
Inputs are the immutable `review_20261128_01/checkpoint_001_main70ps_supplement43`
generation resolved through its `current` link. Input hashes, command,
one-thread environment, timestamps and exit status are retained in
`input_sha256.json` and `job/`. The job ran at nice 19 and exited 0.

## Observations

| Hard-coordinate bin center | PMF kcal/mol | Raw frames | Kish weight ESS | Largest-frame fraction | Exact-vertex mass fraction |
|---:|---:|---:|---:|---:|---:|
| 0.1125 | 4.819 | 399 | 110.08 | 2.60% | 0% |
| 0.1625 | 1.294 | 4,658 | 508.27 | 0.74% | 25.89% |
| 0.2125 | 1.349 | 1,786 | 47.83 | 9.30% | 0% |
| 0.2375 | 1.873 | 1,109 | 8.78 | 27.58% | 0% |
| 0.2625 | 2.900 | 864 | 1.066 | 96.84% | 0% |
| 0.2875 | 6.729 | 8,808 | 5.88 | 27.75% | 9.71% |

The 0.2625 bin is effectively a one-frame reweighting estimate. That frame is
`activation_grid33/activation_support_3_03`, retained index 78,151 in the pooled
array, at 7.60 ps of its own history, with
`(q_acyl,q_act,q_del)=(-2.00176264,-0.81377617,-2.38273186)` Å. It projects
to hard s=0.25498599 on segment 7 and is **not at a vertex**. Its lineage
supplies 97.09% of bin mass. Thus thousands of frames elsewhere do not provide
many effective unbiased observations for this bin.

For bins 0.2125, 0.2375, 0.2625 and 0.2875, respectively,
87.51%, 99.86%, 99.9987% and 99.95% of weighted mass comes from samples at
history age ≤10 ps. Main57 histories supply only 12.28%, 0.00261%, 0.00102%
and 0.0494%. Extending only the older main histories is therefore not equivalent
to substantially extending the histories that currently dominate this portion
of the displayed pooled estimate. These ages use recorded selected indices,
start times and frame intervals; they are not ages since a fully independent
physical starting ensemble. Small differences at exact time boundaries are
floating-point representations of the recorded ages.

## Hard-polyline geometry: a measured contribution, not the whole explanation

The hard projection clips each internal segment parameter to [0,1], then picks
the closest segment. A finite-volume region of continuous three-coordinate
space can therefore project to an **exact vertex**, producing an atom in the
one-dimensional probability measure. This is not a single duplicated frame.

The early corner s=0.1565755076 contains 1,391 frames and 25.8875% of its
0.15–0.175 bin's probability. Removing only that corner atom while holding
all other weights fixed would raise that bin by
`−kT log(1−0.258875)=0.179` kcal/mol (ignoring any reference renormalization,
which cancels for pairwise PMF differences). This is much smaller than the
observed **3.524 kcal/mol** descent from 0.1125 to 0.1625. The corner at
s=0.18787518 contributes only 0.1943% of its bin mass. At s=0.29329410,
7,048 frames project exactly to the corner, but contribute only 9.7096%
of their bin's weighted mass. Raw pile-up and weighted dominance are distinct.

## Coordinate comparison with unchanged statistical weights

The frozen BGP coordinate uses 153 equally spaced path nodes and
lambda=56.36835607868901 Å⁻². Its profile still has an early hump:
3.429 kcal/mol at named BGP s=0.1125 and 2.261 at 0.1875; the feature is
shallower than in the hard-coordinate density. BGP also retains the extremely
poor weighted support around 0.2625 (Kish ESS=1.0663). Changing projection
therefore does not fix concentrated reweighting.

The two named coordinates induce different density measures and have different
reference cells. Their equal numeric s values do not label identical surfaces.
The BGP gradient norm is strongly nonuniform: weighted within-bin means are
3.026 at BGP 0.1125 and 0.408 at 0.1625 (a 7.42-fold ratio), and 0.0448 at its first bin. This
supports substantial coordinate-density distortion. These gradient summaries
are diagnostics, **not a Jacobian-corrected free energy**. In several dimensions
the coarea surface integral matters; dividing the histogram by a single mean
gradient is not justified. A smoother BGP plot does not establish or disprove
a physical metastable basin.

The separately declared initial branch, −2≤q_acyl<−0.8 Å, gives a physical
q_act conditional profile that rises from approximately zero at −1.05 Å to
4.802 at −0.55, 9.329 at −0.25 and 13.160 kcal/mol at +0.25 Å. It does not
show the same early path-coordinate dip. It is a different marginal under an
explicit conditioning restriction; the dominating hard-bin frame at
q_acyl=−2.00176 Å is just outside that restriction. This comparison demonstrates
coordinate/conditioning sensitivity rather than proving the path dip spurious.
The top frames of the hard 0.1125 and 0.1625 bins still have q_act near
−1.02 and −0.91 Å, respectively: advancing hard s does not uniquely imply
the corresponding advance along the reference path's q_act values.

## Interpretation and limits

The evidence supports **three concurrent causes**: nonuniform multidimensional
projection, measured corner atoms, and severe importance-weight concentration
in short corrective histories. Corner atoms alone cannot explain the early
peak-to-dip amplitude. The initial 0.1125 peak itself is not a one-frame event
(ESS110, largest weight2.6%); repeated persistence is insufficient to promote
it to a physical barrier or intermediate while the coordinate measure and
conditional ensemble change its appearance substantially.

Kish ESS here diagnoses unequal weights only; temporal correlation can make
the independent information smaller. This analysis does not establish
stationarity, transition rates, a converged PMF or any accepted activation
barrier. It deliberately makes no new bootstrap confidence claim.

## Artifacts and verification

- `early_diagnosis.png` / `.pdf`: coordinate, atom, weighted-support and physical-coordinate panels.
- `early_bins.json`: complete per-bin counts, weights, age and group summaries.
- `early_vertices.json`: exact vertex locations, counts and weighted mass.
- `dominant_frames.json`: top five frames in each early bin with q and history age.
- `lineage_contributions.json`: every contributing state and mass fraction.
- `bgp_gradient.json`: gradient-range and weighted-mean diagnostics.
- `initial_branch_qact.json`: exact existing physical conditional profile.
- Matching CSVs are table conveniences; JSON is authoritative.

Verification: Ruff format/check passed; exact hard histogram counts and BGP
PMF reproduction assertions passed against the retained baseline. Figure
inspected. No new production dependencies or production changes. No new test
suite added for this isolated data-analysis artifact; code coverage unmeasured.
