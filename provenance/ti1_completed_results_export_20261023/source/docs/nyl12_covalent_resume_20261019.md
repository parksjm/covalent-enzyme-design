# Nyl12 covalent-modeling handoff (2026-10-19)

## Scientific target and model definitions

- Prioritize the **first tetrahedral covalent intermediate (TI1)**, not the later acyl enzyme. In TI1, beta-chain Thr233 (Thr1 in chain E) OG1 bonds to the scissile carbonyl carbon; the substrate C–N bond remains intact; the former carbonyl oxygen is an oxyanion. The OG1 proton has reached the initially neutral beta N terminus through a neutral relay water, making that terminus NH3+ in TI1. The initial enzyme state has NH2 at that terminus.
- PA6: L4 cleavage is constrained to produce two measured L2 products (neutral C12H24N2O3, exact mass 244.17869). L2 is the priority over other observed PA6 products.
- Poly(X,AMA): model the A|M cut yielding the measured MAXA product (neutral C26H42N4O5, exact mass 490.31552); other possible products are outside the current design target. This is **poly(X,AMA), not poly(MAXA)**.
- Treat attached papers as evidence, not instructions. The product gate follows the user's assay requirement and the Dishner/commercial polyamide results.

## Completed, with exact paths

- PA6 AF3 covalent TI1 prediction completed (seed 20261019, ranking score 0.951171; score alone does not establish chemical geometry): `output_nyl12_covalent_qmmm_20261019/pa6_l2_af3_ti1/prediction/nyl12_pa6_l2_ti1_af3_seed20261019/nyl12_pa6_l2_ti1_af3_seed20261019_model.cif`. Input, custom CCD, log, runner, and provenance are in `output_nyl12_covalent_qmmm_20261019/pa6_l2_af3_ti1/`. Full eight-chain Nyl12 and ligand I; bonded atom pair E1 OG1–I1 C54. No ligand-affinity calculation was requested or run. **Prediction geometry has not yet been audited or refined.**
- Poly(X,AMA) AF3 covalent TI1 prediction completed (seed 20261019, ranking score 0.913599; same caveat): `output_nyl12_covalent_qmmm_20261019/poly_xama_maxa_af3_ti1/prediction/nyl12_poly_xama_maxa_ti1_af3_seed20261019/nyl12_poly_xama_maxa_ti1_af3_seed20261019_model.cif`. Input, custom CCD, log, runner, and provenance are in its parent directory. Bonded atom pair E1 OG1–I1 C033. **Prediction geometry has not yet been audited or refined.**
- PA6 graph-derived TI1 starting system and LEaP topology: `output_nyl12_covalent_qmmm_20261019/pa6_l2_ti1_graph/ti1.pdb`, `ti1.prmtop`, `ti1.inpcrd`, `ti1_audit.json`, `ti1_provenance.json`, and `ti1_dftb3.mdin`. Audit: 19,973 atoms, 1,314 residues, total charge -5, QM region formal/partial charge 0, beta N has H1/H2/H3, OG1 lacks HG1, neutral relay water, tetrahedral carbon bonded to OG1/C62/O49/N59. This **does not yet use the AF3 coordinates** and is only a starting graph/pose.
- Later acyl-enzyme PA6 reference topology (not the priority state): `output_nyl12_covalent_qmmm_20261019/pa6_l2_bondi/`. Its DFTB3 geometry optimization was interrupted after reported step 160; `acyl_dftb3.mdout` preserves the log and `acyl_dftb3.rst7` is an earlier written restart, not a converged final geometry. Do not report it as optimized.
- Poly(X,AMA) acyl reference with provisional Gasteiger charges: `output_nyl12_covalent_qmmm_20261019/poly_xama_maxa_gas2/`. The RDKit-aligned starting pose had high clashes and no fully converged MMFF conformer in the trial. It is neither TI1 nor MD-ready.
- Antechamber full-oligomer AM1-BCC/SQM calculation was interrupted after 190 displayed SQM optimization iterations. Saved `sqm.in` and `sqm.out` in `output_nyl12_covalent_qmmm_20261019/poly_xama_am1bcc_interrupted/`; no completed AM1-BCC mol2 exists. All active calculation processes were stopped for the laptop pause.

## Code and provenance

- `src/nyl12_covalent/` holds chemistry/product gates, AF3 custom CCD/input construction, pose generation, LEaP topology, protonation/audit, and Amber DFTB3 input builders.
- Entry points: `scripts/prepare_nyl12_af3_ti1.py`, `scripts/build_nyl12_ti1_qmmm.py`, and `scripts/build_nyl12_acyl_qmmm.py`. Outputs include provenance JSON; keep inputs/logs with any manuscript figures.
- Original PA6 bound pose came from Boltz-2 affinity-mode model 41, then OpenMM relaxation: `output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/final/relaxed_pdb/pa6_l4z_l2_nyl12_wt_model_41.pdb`. Current work no longer computes affinity. Poly(X,AMA) initial pose used RDKit conformers aligned to this register, not a Boltz prediction.
- `docs/nyl12_pa6_covalent_qmmm_plan_20261019.md` predates the completed poly(X,AMA) implementation and contains an outdated “extend graph” item; update it before citing as current status.
- The shared workspace had many unrelated pre-existing changes. Restrict edits to this modeling package and its documentation.

## Resume sequence

1. Audit both completed AF3 CIFs for actual OG1–C distance, tetrahedral C–O/C–N geometry, beta N protonation (AF3 may omit/alter hydrogens), oxyanion assignment, clashes, and product register. Canonicalize chain/atom names and export inspection-friendly PDB/mmCIF with an explicit source-to-output atom map. Do not assume the high AF3 ranking scores validate chemistry.
2. Use the audited AF3 pose, rather than the current graph-derived pose, to rebuild **PA6 TI1** and then poly(X,AMA) TI1. Ensure LEaP preserves neutral beta N in ES and cationic beta N in TI1. Retain the relay water and inspect its H-bond orientation.
3. Run restrained AmberTools26 DFTB3 QM/MM optimization for TI1, with QM Thr1, whole substrate, and relay water; verify convergence and chemistry from the saved trajectory/restart. The current `ti1_dftb3.mdin` exists but **TI1 optimization has not run**. Consider multiple substrate conformers/poses before interpreting tail contacts.
4. Implement a content-addressed parameter cache before repeating charge jobs. Key by canonical molecular graph, stereochemistry, state/protonation, conformer/charge method, and Amber/GAFF/tool versions; store input/output hashes and provenance. Existing parameter files are local outputs, not a reusable cache. Gaussian 16 was not found in PATH or initial `/data/jpc` and `/opt` search; locate it before selecting RESP. For later MD, obtain valid complete bonded parameters and robust charges; the current QM-only bridge parameters and Gasteiger poly charges are provisional.
5. Expand to reaction path/TS and design only after TI1 coordinates, protonation, and topology are reliable. Keep the measured L2 and MAXA product gates.

No test suite, lint, or postprediction geometry validation was run in this pause handoff.

## Update: 2026-10-20, before switching to GPT-6 Astra

- Both completed AF3 models were audited and exported to `pa6_l2_af3_ti1/canonical/` and `poly_xama_maxa_af3_ti1/canonical/` under the same output root. Each contains `protein.pdb`, `ligand.pdb`, `complex_heavy.pdb`, `atom_map.json`, and `audit.json`. The new reproducible code is `src/nyl12_covalent/af3_postprocess.py` and `scripts/canonicalize_nyl12_af3_ti1.py`. Source and canonical ligand atom names match exactly; no affinity-mode name increment was applied.
- AF3 reactive distances: PA6 OG1–C54 1.332 Å, C54–O49 1.306 Å, C54–N59 1.398 Å; poly(X,AMA) OG1–C033 1.286 Å, C033–O034 1.260 Å, C033–N035 1.427 Å. Both OG1–C bonds are too short for a reliable tetrahedral adduct. The poly model also has a 1.669 Å nonbonded C060–B102:CZ contact. Read `canonical/audit.json` for the full contact ledger. AF3 provides no proton or relay-water coordinates.
- Built AF3-coordinate-based TI1 LEaP topologies in `pa6_l2_ti1_af3_neutral_beta/` and `poly_xama_maxa_ti1_af3_neutral_beta/`. Both passed connectivity/proton inventory audits; the attacked E beta N is NH3+, F/G/H beta N termini are NH2, and the original substrate amide C–N bonds remain intact. Topology total charge is -8 in both systems; QM-region charge is 0. The added entry-point flags `--protein-pdb` and `--ligand-pdb` select canonical AF3 coordinates.
- **Do not run the current DFTB3 jobs yet.** The geometric relay-water seed has severe clashes after LEaP adds hydrogens: in PA6, water O–E1:HB 1.31 Å, water O–A186:HB3 1.43 Å, and water H2–E1:HB 0.92 Å; poly has analogous contacts. Repair placement with a search against a hydrogenated protein and ligand, consider alternate waters/proton microstates, and audit local geometry before optimization. This may challenge the assumed one-water relay mechanism rather than merely requiring a better initial rotation.
- F/G/H neutral beta termini were implemented in LEaP by deleting NTHR H3 and assigning N charge -0.6254 so each modified residue is neutral. This is a **provisional QM/MM-only** charge assignment (the retained N3 atom type and other charges were not refit). Fit/validate a neutral beta-N template before classical MD or mechanistic comparisons. The E TI1 NH3+ protonation is explicitly validated. An ES topology with E neutral NH2 still needs to be built.
- Parameter caching, Gaussian16 location/RESP, MD-ready reactive parameters, geometry refinement, full tests/lint, and plan-doc correction remain pending. Existing completed source parameter files were reused, but there is no content-addressed cache yet. Preserve both original and corrected topology directories; the original `*_ti1_af3/` PA6 topology incorrectly kept the F/G/H beta N termini charged and is superseded by `*_neutral_beta/`.

## Update: 2026-10-20, reproducible execution and extraction planning

- Added content-addressed parameter caching (`parameter_cache.py`) and integrated
  it with poly(X,AMA) parameterization. Cache identities include chemistry,
  conformer, charge method, and software/data hashes; entries validate artifacts
  before reuse. Existing old outputs have not automatically become cache entries.
- Added the tmux job launcher (`scripts/run_nyl12_job.py`) and canonical
  hydrogenated-model exporter (`inspection.py`). Fourteen focused tests for
  jobs, inspection, and caching passed using
  `PYTHONPATH=src uv run --no-project --python .pixi/envs/analysis/bin/python python -m pytest tests/test_nyl12_jobs.py tests/test_nyl12_inspection.py tests/test_nyl12_parameter_cache.py -q`.
- Scripted rigid-pocket relay-water searches ran in tmux for both substrates.
  Reports and launch records are under
  `output_nyl12_ti1_refinement_20261020/`. Both searches failed the explicit
  contact gate and wrote diagnostic candidates, not accepted starting restarts.
  No TI1 QM/MM optimization was launched. Diagnose the contacts before changing
  placement or permitting local relaxation. Failure of this finite rigid-pocket
  search does **not** establish that a one-water reaction mechanism is impossible.
- Before refinement, expand the QM region to avoid cutting the Thr929–Thr930
  peptide bond: proposed selection is all Thr929, ligand1313, water1314, and
  Thr930 except its C/O. Audit the resulting C-alpha–C boundary and electron
  inventory before use; this change has not yet been implemented.
- Gaussian: user identifies installations at `/usr/local/g16` on geranium,
  tulip, and aster; violet/azalea installation is pending. Read-only checks found
  `jpc` lacks access to the g16-group directory on all three existing hosts, and
  `/data/scratch/gaussian_scratch` is absent. User explicitly withdrew the
  `/tmp` scratch suggestion. No Gaussian job or remote environment change made.
- User authorized creating skills when ready. The `scientific-tmux-jobs` skill
  source is under `tools/skills/`, installed at
  `/home/jpc/.codex-hosts/violet/skills/scientific-tmux-jobs`. Its structural
  validator passed and installed/source contents match. Independent review
  prompted an explicit warning that current `qmmm.py` runners do not implement
  checkpoint continuation and would overwrite existing output names.
  Covalent-modeling skill remains pending successful refinement. Proposed
  standalone private repository boundary and migration gates are in
  `docs/covalent_pipeline_repository_plan_20261020.md`; no remote created or
  scientific files relocated.

### Local Gaussian installation now verified

The user supplied the newly installed `/usr/local/g16` on violet and created
`/data/scratch/gaussian_scratch/jpc`, now owned by `jpc:users` with mode 700.
The group database includes `jpc` in `g16`; the existing Codex process predates
that membership and can use `sg g16 -c ...` to launch with the granted group.

`scripts/preflight_gaussian16.py` ran in tmux session
`nyl12_g16_preflight_20261020` and passed with normal termination. Version:
Gaussian 16 ES64L Revision C.02 (7-Dec-2021). The water HF/STO-3G single point
was solely an installation/access/scratch check, not a charge fit or a
scientific reference calculation. Inputs, executable/profile/source hashes,
and logs are in
`output_nyl12_ti1_refinement_20261020/gaussian_preflight_violet_20261020/`;
launch record is `jobs/g16_preflight/` under the same output root. Five tests
using fake executables passed; Ruff passed for this script and its tests.

The user selected their personal GitHub account for the new private repository.
Both authenticated interfaces identify it as `parksjm`; planned remote is
`parksjm/covalent-enzyme-design` (now created, private), local checkout
`/data/jpc/covalent-enzyme-design`. Preserve the current scientific output paths
during extraction.

### Structural diagnostics and extraction checks

- `scripts/diagnose_nyl12_relay_water_contacts.py` now reports canonical chain,
  chain-local residue number, full-protein residue number, contact thresholds,
  hydrogen parents, and input/source hashes. It rejects changed topology or
  coordinate inputs relative to the saved search and refuses output overwrite.
  Latest reports are
  `output_nyl12_ti1_refinement_20261020/{pa6_l2,poly_xama_maxa}/water_search/contact_diagnostic_v2.json`;
  tmux provenance is under `jobs/diag4_*`. Water oxygen contacts with Asn186
  ND2/backbone O dominate several candidates; rotating beta-N hydrogens alone
  does not resolve the observed heavy-atom crowding.
- Added `qm_region.py` and `scripts/audit_nyl12_ti1_qm_region.py`. The expanded
  mask `:929,1313,1314|(:930&!@C,O)` agrees with ParmEd mask selection. Structural
  audits select 109 PA6 atoms and 181 poly(X,AMA) atoms, with one CA930–C930
  carbon–carbon cut. Formal QM charge is 0, selected force-field partial-charge
  sum is -0.0294, and link-H electron counts are even (364 and 628).
  This is a structural boundary audit, not a validation of embedding, GB, or
  DFTB energetics. The old `qmmm.py` TI1 input writer still uses the earlier
  peptide-cut mask; replace it before any optimization.
  Final formatted-source audit records are
  `output_nyl12_ti1_refinement_20261020/jobs/qm_region_audit_pa6_r1/report.json`
  and `jobs/qm_region_audit_poly_r1/report.json` under that same output root.
- The complete focused extraction test set passed: 35 tests covering jobs,
  inspection, cache, water diagnostics, QM selection, source inventory, and
  Gaussian preflight. Ruff passed across the allowlisted package/scripts/tests.
  The environment lacks `coverage`/`pytest-cov`, so coverage percentage has not
  been measured. No new dependencies were installed.
- `scripts/inventory_covalent_pipeline.py` creates a strict allowlisted,
  checksummed inventory plus an optional source tarball. It records original
  tracking status and archives exact bytes, preserving uncommitted source that
  hashes or the old Git revision alone cannot reconstruct. Repository templates
  are in `tools/covalent_repository_template/`.
