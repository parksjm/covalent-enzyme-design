# mol_blender runtime on Violet

Installed from the official [mol_blender repository](https://github.com/magnusbauer/mol_blender/tree/f7ace1ebbaedfc4bc805489ce5777e09c6f43795) at commit `f7ace1ebbaedfc4bc805489ce5777e09c6f43795`. The checkout, uv environment, Blender runtime, and configuration are isolated under:

- Checkout: `tools/mol_blender/`
- uv environment: `tools/mol_blender/.venv/`
- `MOL_BLENDER_HOME`: `tools/mol_blender/.runtime/`
- Blender: 5.0.1; Molecular Nodes 4.5.10; NumPy 1.26.4; Biotite 1.6.0
- CCD SHA-256: `dffdb3a19600d51854058daecfbc9bd1b94aab413a04e7afed8d4b7670621d47`
- PyMOL interpreter: `/data/jpc/code/miniforge3/envs/nyl12-pymol/bin/python` (PyMOL 3.1.0)

The checkout requests Python 3.11.15, unavailable on Violet. The pinned renderer requires the Python 3.11 ABI, so setup used the installed managed Python 3.11.13 explicitly; no upstream lock or project `pyproject.toml` was changed. No production dependency was added.

## Setup and rendering

Setup ran in tmux session `mol_blender_setup_20261125_violet_retry01`; its script and full log are under `docs/presentations/nyl12_progress_20261123/assets/mol_blender/attempt_20261125_1545_retry01/`. The final validated render ran in tmux session `mol_blender_final_grayH_20261125_violet` on CPUs 4 and 12 with CUDA unused. Its exact invocation, compositor script, source hashes, and combined output are retained in `docs/presentations/nyl12_progress_20261123/assets/mol_blender/attempt_20261125_1620_retry01/`.

For the retained relay candidate, the style derives from `vivid_binding_site`; use the schema-valid overrides in `attempt_20261125_1620_retry01/relay_vivid_thinsticks.json` and generator `build_style.py`. Binding-site sticks are 0.18 Å (the pinned Molecular Nodes `StyleSticks` API documents Angstroms), rendered at 800×800 and 32 samples. The input contains disconnected local residues, so the style omits a misleading ribbon trace. Hydrogen display hints and palette are both `#B8B8B8`; the stick materials use roughness 0.58 and specular 0.12 with the preset’s area lights for soft shadows.

The white background is set to solid `#FFFFFF` and composited over white with Blender's native Alpha Over node. The final scene sets `film_transparent` only for the molecular foreground, disables output dithering, and saves the compositor in the editable `.blend`. The final PNG has four corner pixels equal to opaque RGBA `(255, 255, 255, 255)`, 67,442 nonwhite foreground pixels in bounding box `(69, 86, 654, 673)`, and zero transparent pixels. `render.sh` enforces the white and nonblank checks.

## Structure and bond provenance

Source PyMOL input and acceptance are retained in `docs/presentations/nyl12_progress_20261123/assets/relay_state_pymol/`. `attempt_20261125_1555_retry02/relay_candidate.pse` was derived by recalling the saved `relay` scene; coordinates agree with `input.json` within 4.71×10⁻⁷ Å and the saved camera round-trips exactly. That view keeps Thr233 below, Phe113 upper left, and Asn186 upper right.

The official PSE exporter’s mmCIF text contains no `_struct_conn` records. Its PDB text retains the 148 explicit relay CONECT edges and 152 per-atom PyMOL color hints. The rendering bundle therefore leaves `cif_text` empty and uses the PDB text. `bond_color_audit.json` records the intended relay changes from the 147-edge base graph: OG1–HG1 is removed; OG1–scissile-C and N–HG1 are present. Native topology and restart files were not edited.

The initial schema-rejected style and the first two nonfinal PNGs remain in their attempt directories with logs. Only the gray-hydrogen white-composited output that passes the pixel check is final: `attempt_20261125_1620_retry01/relay_candidate_white_greyH.png` (SHA-256 `27441b76d47e4701f2b7ad5d7003abd60697f64119882a93a4b1361140d140fc`) and `relay_candidate_white_greyH.blend` (SHA-256 `0e8fb66888c79d64b2e40d69cf57303ae4b6f9d93e2291e576acb6559b39f442`). Its resolved style, runtime manifest, color/bond audits, source hashes, launch scripts, and logs are in the same attempt directory.
