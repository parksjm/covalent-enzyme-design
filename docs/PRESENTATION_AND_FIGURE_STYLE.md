# Scientific presentations, plots and molecular images

PI preferences, adopted 2026-11-25. These conventions apply to project PPTX
files, slide figures and publication images. Journal-specific requirements
govern final manuscript exports. For interactive trajectories, also follow
`VMD_MOVIE_CONVENTIONS.md`.

## Narrative and slide selection

- Begin with a **proper title slide**: scientific project title, short purpose
  and the current build date. Include authors/affiliations only when supplied.
  The presentation filename, title-slide date and build metadata must use the
  actual build date in **America/New_York**, with `YYYYMMDD` in the filename.
  Generate the date in the builder; never leave a hardcoded prior date in the
  output filename or footer. Evidence timestamps remain their actual dates:
  rebuilding a deck does not update the observations it describes.
- State one main point per slide. Use a **declarative sentence as the title**
  whenever it sounds natural in a scientific presentation. A proper cover and
  concise technical topic titles are appropriate; do not force awkward claim
  sentences merely to satisfy this preference.
- Establish the experimental motivation, enzyme/substrate model, reaction
  scope and ultimate goal of improved catalytic efficiency before detailed
  results. Present the methodological progression: force-field parameterization,
  WT MD/protonation studies, MM/QM-MM preparation, QM method/region checks,
  reaction paths, free-energy sampling/analysis and TS-based sequence redesign.
- Use computational-chemistry terminology. Avoid internal execution phrases
  such as “retained endpoints”, “relay seed connects” or “orientation test,
  not relay proof” as slide titles. Explain E–S, acylenzyme, PES, PMF, TS,
  geometric metrics and any screening categories before relying on them.
  Show both substrate directions and their actual proton-transfer geometry;
  “flipped”, “Tier A/B” and “N-terminal delivery” are not self-explanatory.
- Name the proposed channels consistently: **Mechanism A, direct proton
  transfer**, and **Mechanism B, N-terminal proton relay**. Show each proposal
  before its data, with clear placeholders for unfinished calculations. Do
  not use “older direct route” as its scientific identity. Distinguish direct-
  seeded mixed-channel sampling from a qualified mechanism-specific PMF.
- Explain why alternate protonation states were considered. Cite Negoro et al.
  (2023), DOI `10.1111/febs.16755`, for the proposed NylC Asp protonation/triad;
  distinguish mechanistic inference from direct protonation evidence and
  maturation-confounded mutagenesis from mature substrate hydrolysis.
- Use figures or tables to make the argument. Keep additional words to the
  minimum needed to interpret the display. Move explanation to speaker notes.
- Present the current scientific story, not a chronological execution log.
  Remove obsolete failed scans and irrelevant debugging plots from the main
  presentation, including the approximately 95 kcal/mol driving example.
  Preserve negative results and their provenance on disk. Retain a negative
  result in the deck only when it explains a material scientific decision;
  do not hide contradictory evidence or imply an unsuccessful route succeeded.
- Separate observations, hypotheses and proposals. The 16.3 kcal/mol relay PES
  feature is not a PMF or a certified TS barrier. A stable TI1 is optional.
  L4→L2 is a modeling hypothesis, not an established substrate/product mapping.
  Lysate/MS measurements and discordant K59E replicates require qualification.

## Layout and typography

- Use a consistent 16:9 layout, fonts, margins, palette and title placement.
  Keep titles to one or two lines, normally **34–40 pt**.
- Body text is normally **24–28 pt or larger**. Do not shrink text to make a
  crowded slide fit; remove words or split the content. Avoid dense bullet lists.
- Tables normally use **22–24 pt** text, with few columns and roughly six
  data rows. Keep necessary labels readable; use notes or an appendix for detail.
- Reserve clear space between title, figure, captions and footer. Prevent
  wrapped units, clipped labels, orphaned words and overlaps. Do not place tiny
  footnotes over a figure. Short source attribution can go in notes.
- Keep text, tables and labels editable where practical. Preserve the original
  figure data and plotting source alongside rendered images.

## Plots

### PMFs and umbrella-sampling diagnostics

**Primary relay protocol revision, 26 November:** the PI requested a defensible
replacement for the latest-of-thirteen truncation heuristic. New publications
use a fixed additional 5 ps production discard once/window, retain every
subsequent frame, and display pointwise 95% intervals from 200 within-window
contiguous 5 ps block-bootstrap MBAR refits. Separate heating and 5 ps native
equilibration remain excluded. Compare 0/2.5/10 ps production discards, early
and late halves, and 2.5 versus 5 ps bootstrap blocks. No IID analytical band
is valid for the correlated all-frame estimate. No automatic cutoff or
universal duration certifies equilibration. Bins lacking support in any
bootstrap draw have no displayed interval. These estimates remain provisional
until temporal, orthogonal and initialization checks pass. Use
`analyze_relay_fixed_pmf.py`, `plot_relay_fixed_pmf.py` and the bounded
`watch_relay_fixed_pmf.py` publisher. Keep all prior profiles and their methods
explicitly archived; the historical instructions below describe those earlier
estimators and do not override this revision.

**Continuous-window update, 26 November:** at the PI's request, subsequent
relay PMFs concatenate the uninterrupted production segments for each window
before estimating **one** equilibration discard and correlation factor.
Use `analyze_relay_continuous_pmf.py`; preserve the former blockwise analyses.
For an analysis-policy comparison, show **matched PyMBAR analytical pointwise
95% intervals for both selections**, keeping the coordinate, bins and E–S
reference identical. Supplement these with 200 within-window bootstrap refits
of the whole-window-decorrelated configurations (`window_bootstrap`). These
assume effective independence after thinning. They differ from the archived
completed-5 ps block bootstrap and must be labeled accordingly. Eight windows
at the 20 ps checkpoint retain less than 5 ps after their single discard;
do not disguise these short records as multiple independent 5 ps blocks.
The continuous-window watcher replaces the blockwise figure watcher for
`*_latest` publications; the older source and figures remain archived.

- Use the shared `relay_pmf_uncertainty`, `relay_pmf_block_bootstrap` and
  `relay_pmf_plotting` modules and their `scripts/analyze_relay_pmf_uncertainty.py`
  and `scripts/plot_relay_pmf_diagnostics.py` entry points. The bounded
  `scripts/watch_relay_pmf_figures.py` publishes the same figures at subsequent
  completed 5 ps cumulative checkpoints. Preserve the source report, request,
  source hashes, uncertainty JSON, PNG and PDF. Use dated, horizon-specific
  filenames and stable `*_latest` copies for inspection.
- Within a campaign, keep coordinate definition, histogram edges, reference
  bin, units, colors and labels fixed across checkpoints. Do not smooth the
  primary curve or change bin widths to make it appear converged. Compare
  different campaigns only after reconciling their coordinates and estimators.
- Plot pointwise **95% statistical intervals**. State the method explicitly:
  histogram asymptotic covariance or completed-block bootstrap with MBAR
  refitting. For the relay campaign, use **200 within-window, completed 5 ps
  block-bootstrap refits** as the displayed interval and retain analytical
  errors as a cross-check. Preserve the original per-block decorrelation and
  discard choices when resampling blocks; do not bootstrap fixed MBAR weights.
  Keep the E–S reference bin fixed in every replicate, rather than independently
  moving the zero to each replicate's minimum.
- Retain missing bin/reference support in bootstrap diagnostics. Show bootstrap
  intervals only where all replicates have support; leave unsupported estimates
  or intervals missing. Never hide failed resamples by plotting successful-only
  percentiles without the support fraction. Intervals are pointwise, not a
  simultaneous confidence band or a confidence interval for a selected maximum.
- State that statistical intervals are conditional on sampled configurations
  and estimated decorrelation. Four completed 5 ps blocks/window do not exclude
  longer correlations, hidden chemical branches or preparation bias. Neither
  a smooth curve nor a narrow interval establishes equilibrium or a validated
  activation barrier. Preserve the separate scientific acceptance gates.
- Distinguish three support measures: selected/decorrelated frames per bin;
  **bin-local weight ESS**, `(sum w)^2 / sum(w^2)` within that bin; and number
  of contributing windows. A bin may contain many selected frames but have
  weight ESS near one. Global target-state ESS is another statistic and does
  not describe the precision of every PMF bin.
- Draw neighboring-window comparisons as **unconnected points**, with exact
  window-pair identities in metadata and labels for nonadjacent manifest
  comparisons. A plot line must not imply an unmeasured connection. Distinguish
  empirical 3CV histogram intersection from the directed MBAR overlap matrix;
  their numerical thresholds are different. Compare histogram-overlap trends
  across time on a fixed physical-CV grid; pooled-data-dependent grids are
  valid within a checkpoint but are not a matched-grid trend analysis.
- Use support and overlap together to plan sampling: sparse bins motivate
  targeted extensions of the windows visiting them; low weight ESS despite
  many frames motivates improved sampling/reweighting support and inspection
  of slow orthogonal coordinates. Add bridging windows only for persistent
  distribution gaps. Changing the number of display bins does not create
  missing equilibrium sampling.

### General figure conventions

- Size plots for their actual placement. Target axis labels **24–28 pt** and
  tick/legend labels **20–24 pt after scaling onto the slide**. A large source
  font becomes small if the image is reduced to a thumbnail.
- Use correct symbols and units: **ΔG‡**, **ΔE‡**, **φ**, **α_BD**, **Å**,
  **kcal mol⁻¹**, **K** and **ps**. Use proper subscripts/superscripts or
  mathtext rather than ASCII substitutes in final figures.
- Distinguish physical PES energy, restraint-inclusive energy and PMF. State
  reference state, physical temperature, sampling scope and uncertainty when
  relevant. Do not equate sampling duration with convergence.
- Add concise, preferably editable **plot annotations** when a scientific
  qualification is essential: “TS verification pending”, “Mixed proton-transfer
  channels”, “Sampling not converged”, or “PMF pending”. Place text boxes in
  clear plot space or beside the figure with an unambiguous pointer; preserve
  curves, error bars, ticks and axis labels. Use “unconverged TS” only when the
  TS optimization itself is unconverged: a converged NEB band can still lack
  a verified first-order saddle. Do not rely solely on speaker notes for these
  distinctions. The PI explicitly approved these callouts on 2026-11-25.
- Use consistent, colorblind-compatible group colors across related plots.
  Keep WT/control colors and state labels stable; avoid decorative gradients.
- Reduce grid clutter and unnecessary legend entries. Preserve honest axes
  and show uncertainty where available. Do not conceal missing windows or
  interpolate unsupported PMF segments.
- Retain vector PDF/SVG for plots where practical and a high-resolution PNG
  for slides. Check the exported file, not just the plotting code's font values.

## Molecular images

- Use **PyMOL** or the optional **mol_blender** renderer for publication and
  presentation images. VMD remains useful for trajectory inspection. The PI
  approved a Violet mol_blender installation/trial on 2026-11-25; installation
  and runtime verification are recorded in `MOL_BLENDER_RUNTIME.md`. It is not
  a scientific-pipeline dependency. Obtain approval before adding other software.
- For mol_blender start with `vivid_binding_site` or `northern_nights` (the
  package name for the requested northern-lights style). Reduce stick thickness
  toward the 0.2 Å display convention after checking the renderer's units.
  Preserve orange substrate C, blue N and red O rather than blindly adopting
  preset element colors. Use light-grey H as specified below. If ribbons are shown, color by explicit
  chain ID: A/E green, B/F cyan, C/G magenta and D/H yellow for the mature
  alpha/beta assembly; record any canonical display relabeling. Verify geometry,
  camera and chain selections in an actual rendered preview.
- **Mol_blender backgrounds must always be opaque pure white (`#FFFFFF`).**
  Override preset background colors, gradients and compositor backgrounds.
  Verify the actual PNG background pixels are white; a white setting in a
  JSON file alone is insufficient. The PI made this mandatory on 2026-11-25.
- Against that white background, use **light-grey mol_blender hydrogens**
  (initial target `#B8B8B8`) with soft contact shadows and restrained glare.
  Verify that transferring and polar hydrogens remain visible in the actual
  image. This PI-approved mol_blender override does not change white-H VMD
  defaults. When following PyMOL atom colors, update the hydrogen display
  colors there too so color-follow does not override the grey palette.
- Use a clean white or explicitly requested transparent background, appropriate
  antialiasing and sufficient export resolution. Avoid screenshots of GUI panels.
- **Default active-site camera for every figure and movie unless the PI
  explicitly requests otherwise:** Thr233 below the reactive substrate amide,
  Phe113 upper left and Asn186 upper right, with the Thr233 N-terminus on the
  right as in the PI's 27 November reference image. Define the camera from
  protein landmarks and preserve it across matched states and pose classes.
  Use proper rigid rotations only, never reflection/mirroring. Do not rotate
  each ligand separately to make opposite substrate directions look alike.
  In the approved Boltz pose the substrate amide N and Thr Nα appear on the
  right in this camera; this visual arrangement does not define a universal
  torsion sign. Matched images share scale. Use `set valence, 0` in PyMOL.
- Use orange substrate carbons, cyan protein carbons, white H, blue N and red O
  unless a specific comparison needs a declared alternative. Use few large
  labels, positioned away from atoms and other labels.
- Show the physical connectivity appropriate to each state. Source E–S topology
  bonds must not create long phantom bonds in acylenzyme/intermediate images.
  Record visualization distance criteria; rendered sticks do not determine
  electronic bond order. Never add virtual QM/MM caps as real enzyme atoms.
- Save the renderer/PML/PY script, atom map, source hashes, frame identifiers
  and PyMOL session where practical. Document any rigid display transformation
  or hidden distal substrate portions; never optimize coordinates for appearance.

## Speaker notes

For any Herschlag-phi comparison, state the ordered convention. Our function
uses five atoms: Thr OG1, carbonyl O, carbonyl C, substrate amide N and the
carbon attached to carbonyl C on its acyl side. Carbonyl C is projected onto
the plane defined by O, amide N and acyl-side C. The signed torsion is then
**OG1 → O → projected C → acyl-side C**, evaluated with the unchanged project
`atan2` convention. Thr Nα is not part of phi. In the approved B/F Boltz
reference this is OG1/O50/C55/N57/C63, giving +81.196° before and +88.228°
after relaxation. The approved E–S target is positive; the prior negative-only
selection rule was an interpretation error. Never change measured signs to
match a label, and never infer the sign from screen-left/right alone.

Use notes for definitions, abbreviations, full references, source paths/frame
IDs, methods/settings, controls, uncertainties and interpretation limits. Include
the take-home point and what the figure does and does not establish. For Eyring
tables, state κ=1, elementary-step rate interpretation, and that crude endpoint
conversion does not directly determine k_cat or an experimental activation barrier.
Separate 300 K simulations from the 348.15 K lysate assay and from other purified
enzyme assay conditions. Do not mix their measurements or temperatures.

## Build, review and delivery

1. Preserve the prior deck/build sources in central `clutter` before replacement.
   Keep one clearly identified current deck at its established logical path.
2. Build from reviewed source data; record timestamps, source/asset hashes and
   software environment. Do not install new production dependencies silently.
3. Render **every slide** to PDF/PNG, inspect the full deck/contact sheet and
   examine dense slides at full size. Correct actual rendering failures rather
   than rely only on shape-bound or mock tests. Check unusual glyphs and fonts.
4. Check speaker notes, numerical values, energy labels, source provenance and
   on-slide text size. Record limitations and which visual checks were performed.
5. Deliver the editable PPTX, a PDF preview, and reproducible figure/build inputs.
   Ensure the curated pull script downloads the current deck and policy.


### Sequence-protection figures (PI revision, 27 November 2026)

Use vivid magenta for other protected identities, very light gray for variable
identities, and teal for the F113/K156/N186/T233/D272/D274 mechanistic subset.
The current disjoint counts are 95+6+227; the 101 protected total includes teal.
Label the subset as mechanistic constraints, not six independently proven
Nyl12 catalytic residues. Use a white background and lighting that retains
the light-gray/magenta contrast. PA6 carbon stays orange. This explicit
protection-status view overrides chain coloring for this figure only.
