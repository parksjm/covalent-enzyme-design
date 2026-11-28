# Nyl12 mechanism and redesign briefing

The editable source is `slides.json`. This build contains 29 main slides and 3
appendices. Scientific evidence is current through **25 November 2026**; the
cover and footer build date is generated from the machine's
`America/New_York` date.

The dated outputs for this build are
`Nyl12_PA6_progress_20261125.pptx`, matching PDF and hashed build receipt.
Numbered slide renders and a contact sheet are in `rendered/`. The stable
working-asset directory remains `nyl12_progress_20261123/` so existing source
references continue to resolve.

## Scientific status

- PA6 is nylon-6. The force-field release is checked for a zwitterionic,
four-repeat PA6 L4 model; it does not validate reactive states or arbitrary
polymer lengths.
- Matched wild-type protonation-state simulations found no geometry-occupancy
advantage for the Asp274-COOH state. These fixed-charge structural occupancies
are not activity or pKa populations.
- The direct-transfer umbrella dataset is a preliminary mixed-channel,
unconverged q_acyl marginal. The N-terminal proton-relay campaign has no PMF
yet.
- A coordinate-driven relay path has a physical potential-energy maximum of
about 16.31 kcal/mol. This is ΔE, not ΔG‡. Saddle curvature and downhill
connectivity are unverified.
- The addition-like relay structure came from direct-seeded umbrella MD with
relay-like configurations, then bias removal and QM/MM optimization. It is not
an established stable tetrahedral intermediate.
- A compatibility audit found the A/B physical models and frozen references
agree within restart precision; numerical tolerance checks and shared-coordinate
free-energy reanalysis remain pending.
- Water-mediated deacylation and complete product release remain a separate
follow-up. The 30 November redesign pilot is a planning target; activity is
untested.

Every slide carries detailed scientific notes, definitions, sources and
limitations. `asset_sources.json` maps each used figure to its retained source
and SHA-256. The build receipt records slide, source, builder, figure and output
hashes. `MANIFEST.sha256` covers the current dated deliverables, used figures,
source files and all slide previews.

## Rebuild

From the repository root, run:

```bash
bash docs/presentations/nyl12_progress_20261123/rebuild.sh
```

The script archives prior dated deck outputs and renders in central `clutter`,
uses the existing offline `uv` report environment for `python-pptx`, and
converts the editable PowerPoint to PDF with LibreOffice. It also exports every
PDF page to PNG and creates the contact sheet. The filename date is computed
on each build in `America/New_York`; it is not hardcoded in the builder. The
scientific `as_of` date in `slides.json` changes only when evidence is reviewed.
No production dependency or scientific calculation is added by this workflow.
