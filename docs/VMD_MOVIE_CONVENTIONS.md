# Reaction-path and QM/MM movie conventions

PI-requested defaults, 2026-11-25. Apply to reaction-path, NEB, umbrella and
optimization movies. Ordered umbrella windows and optimized endpoints are not
continuous reactive dynamics; preserve each movie's provenance and ordering.

## Shared implementation

Use `nyl12_covalent.vmd_movie_policy.render_movie_policy` with verified physical
QM/substrate indices and known C–C edges. It embeds self-contained Tcl, so laptop
viewers need no installed helper. Use `mapped_cc_bonds` for retained native edges
and `vmd_chain_groups.display_chain_groups` for validated chain selections.
Reject ambiguous maps. Never insert virtual QM/MM link caps as physical atoms.

## Representations and colors

- Licorice and DynamicBonds bond radius: **0.2 Å**, consistently across QM,
  substrate and displayed context. Proton VDW scale: **0.2**; VDW scale is not
  an absolute bond radius. Preserve each route's declared distance cutoffs.
- Show QM-residue **CA atoms and CA–CB bonds** (Gly has no CB), plus the native
  C–C connections between QM and MM substrate regions. Add missing known C–C
  edges symmetrically, once, in VMD memory. Do not infer arbitrary carbon bonds
  from close contacts or restore intentionally removed reactive C–N/C–O/N–H/O–H.
- Use a connected static selection, avoiding gaps from split representations.
  Substrate **all carbons, including QM carbons, are orange**. Protein carbons
  remain cyan; H white, N blue, O red. Apply the same colors to dynamic bonds.
  Carbon ends of substrate C–H/C–N/C–O bonds are orange; other ends retain their
  element colors. Both ends of substrate C–C bonds are orange.
- Do not overlay orange carbon-only sticks on cyan substrate sticks: that
  leaves inconsistent half-bonds. The policy instead uses registered display
  Type names. **VMD Type coloring uses the first nonnumeric character**, so
  every display type must start with its color bucket: C (cyan protein C),
  Q (orange substrate C), H (white), N (blue), O (red), and the declared
  additional buckets. A common `NYL12_` prefix makes all atoms use N/blue,
  regardless of separately registered full-name colors. The namespace follows
  the bucket, for example `Q_NYL12_MOVIE_SUB_C`. Register the single-character
  bucket colors. Atom names/elements and native files are unchanged. Original Amber
  types remain in `::nyl12_movie_policy::original_types($molid)` and are restored
  with `::nyl12_movie_policy::restore_types $molid`. Do not export display types
  as force-field parameters.
- Keep route-specific distance/angle labels, camera and frame order. Set
  **`display depthcue off`**. Preserve the small transferring-proton sphere;
  a floating proton is acceptable. A displayed bond is not a computed bond order.

## Hidden chain cartoons

Create separate NewCartoons, **off by default**:

| Alpha chain | Color | Corresponding mature beta chain |
|---|---|---|
| A | Green | E, green |
| B | Cyan | F, cyan |
| C | Magenta | G, magenta |
| D | Yellow | H, yellow |

The current Amber files have blank chain IDs. Mature Nyl12 contains four
232-residue alpha chains followed by four 96-residue beta chains, not four
contiguous precursors. Reconstruct display groups only after verifying repeat
sequences, Thr beta starts and peptide breaks. Topology residue positions are:
A 1–232; B 233–464; C 465–696; D 697–928; E 929–1024; F 1025–1120;
G 1121–1216; H 1217–1312. This is not a new canonical numbering convention.
Preserve complete native labels when available; otherwise use explicit mapped
atom selections. The helper does not rewrite native topology chain attributes.

## Verification and maintenance

1. Execute the strict Tcl fixture: C–C repair without duplicates, CA/CB inclusion,
   effective per-atom colors using VMD's first-character Type lookup, radii,
   hidden cartoons, molecule scoping and repeated use. Include a regression
   reproducing the all-blue result of the obsolete common-prefix scheme.
   Register custom colors explicitly and validate every shade against VMD's
   `colorinfo colors` **before mutating the scene**. PyMOL and VMD do not share
   the same color vocabulary: `brown` is invalid in this VMD palette; use the
   validated `tan`. Test the complete published script and real pull filters,
   not only isolated helper procedures.
2. Validate real topology/index maps and portable checksums. Record whether
   an actual VMD image was rendered: headless Tcl tests do not establish visual
   quality. VMD currently runs on the PI's laptop, not Violet.
3. Refresh viewers and download manifests together. Preserve old small display
   files/checksums in `clutter` with source/output hashes and rollback records.
   Never edit frozen scientific inputs or regenerate large trajectories merely
   to restyle them. Keep existing logical movie paths.
4. Publish the latest US movie at **`vmd_progress/current/view_path.vmd`**.
   Record its actual sampled horizon, production-block count and frames/window
   in its manifest, deriving them from accepted source receipts. Keep the
   user-facing path stable when exporting additional completed blocks; archive
   superseded bundles under central `clutter`. Preserve portable redirects from
   previously published viewer paths. PMF and movie scopes must each be explicit.
5. Pull complete curated viewer/topology/trajectory bundles. Exclude native
   Hessian columns, raw full-system MD and archived display backups.

The corrected viewer prints `NYL12 DISPLAY OK version=20261126-typeprefix3 ...
radius=0.2` after successful application. An explicit DISPLAY FAILED message is
an error, not an accepted visual check. The original palette failure and its
repair are retained in `../provenance/vmd_palette_failure_repair_20261125.json`.
The subsequent Type-prefix repair, 21-viewer refresh and frozen replacement
watcher are recorded under `../output_vmd_typeprefix_repair_20261126/`.

Official VMD references: [color commands](https://www.ks.uiuc.edu/Research/vmd/current/ug/node124.html)
and [in-memory atom attributes](https://www.ks.uiuc.edu/Research/vmd/doxygen/AtomSel_8C-source.html).
The first-character Type lookup is explicit in
[MoleculeList::add_color_names](https://www.ks.uiuc.edu/Research/vmd/doxygen/MoleculeList_8C-source.html).
