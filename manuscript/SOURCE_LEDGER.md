# Manuscript source and claim ledger — 2026-10-20

Paths below are source records, not interchangeable datasets. File hashes and
external repository revisions are in `../provenance/manuscript_sources_20261020.json`.
The original manuscript files were read and remain unchanged.

| Key | Source | Reuse and limits |
|---|---|---|
| S1 | `/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/manuscript_modeling_section_20261011/METHODS_SECTION_DRAFT.md` | Mature assembly, product/register framing and historical modeling methods. Its Boltz, GAFF/AM1-BCC, MD and Rosetta protocols are historical and must not be silently substituted for the new covalent workflow. |
| S2 | `/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/NYL12_PA6_PA66_FINAL_INTERPRETATION_AND_NEXT_PLAN_20260928.md` | Bound-start MD and score/activity limitations; invalidated PA66 histories. No quantitative activity conclusions imported into the new Introduction. |
| S3 | `/data/jpc/enzyme-substrate-modeling/manuscript_nyl12_wt_vs_mutants/NYL12_WT_MUTANT_PROTOCOL_NARRATIVE_AND_AF3_RESTRAINT_TEST_20261009.md` | Historical AF3 restraint work and interpretation. Its evolving dated sections must be read chronologically. |
| S4 | `/data/jpc/enzyme-substrate-modeling/output_nyl12_af3_contact_pa6_l4z_l2_four_site_static_20261009/` | Production RUN_NOTES, input/canonical manifests, and pre/post-relax geometry summaries corroborate 90 models, 360 sites and the implemented AF3 capability. Historical single-site atom naming differs from current TI1. |
| S5 | `../docs/standalone_handoff_20261020.md`, `../docs/ti1_attempt_protocol_20261020.md`, `../provenance/resume_20261020.json` | Current audited topology/preflight status; original water clashes still unresolved. Latest dated checkpoint supersedes prior job-state statements. |
| S6 | `../output_nyl12_resume_20261020/confidence_audit_v2/report.json` | Literal MAXA|MAXA graph/products and ligand-local confidence with original AF3 hashes. Initial failed precision-check attempt retained separately. |
| S7 | `../output_nyl12_resume_20261020/active_site_baseline/pa6_l2.json` and `poly_xama_maxa.json` | Descriptive unrefined geometry; includes hashes of current metric implementation and original projected-phi math. |
| S8 | `../output_nyl12_resume_20261020/neutral_thr_nme_wb97xd/` | Prepared Gaussian inputs, source snapshots, job status and live calculation logs. Launch is verified; acceptance remains pending. |
| S9 | `/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/sander/dlfind_module.F90`, `sander.F90`, and adjacent `libdlfind/src/dlf_formstep.F90` | Native DL-FIND capability, numerical Hessian and mode-selection evidence. Source inspection and library resolution are verified; runtime saddle-search validation is pending. |

## Newly supplied information and verification

- **User instruction, 2026-10-20:** the contact-restraint fork exists at
  `/data/jpc/code/alphafold3-contact-restraints`, has tests and successful
  production use, and remains available even though independent ES prediction
  should normally omit restraints. Verified fork HEAD:
  `ef5dc91d65b4845f63e95ea03316f1613c6dfb5e`, clean working tree.
- **Production evidence:** S4 pre-relax summary reports 90 models, 360 sites,
  250 Tier A/B sites and 35 all-four-pass models. Post-relax summary reports
  193 Tier A/B sites and 11 all-four-pass models. Both carry original AF3
  all-restraints-satisfied diagnostics for 90 models; the latter is not a
  fresh post-relax contact-restraint measurement. These are historical campaign
  outputs, not new validation of enzyme turnover.
- **User instruction, 2026-10-20:** an unrestrained substrate may choose any
  of four equivalent sites. Verified code maps AE, BF, CG and DH. Future
  independent ES scoring must be symmetry-aware; the current AE-only TI1
  metric CLI is insufficient for that purpose. Record all sites and avoid
  treating within-model sites as independent replicates.
- **Atom-name caveat:** historical restraint analysis uses C12/O4/N3/C11;
  current PA6 TI1 uses C54/O49/N59/C62. The earlier WT RUN_NOTES describes the
  graph-checked carbonyl crosswalk. Preserve an explicit mapping for every
  state; do not copy atom selectors by assumption.
- **User preference, 2026-10-20:** wB97XD/def2TZVP geometry, no extra dispersion;
  HF/6-31G(d) ESP for planned RESP. DFTB for initial paths; calibrated chemistry
  and Hessian/endpoint validation are required for TS claims.
- **User instruction, 2026-10-20:** maintain manuscript-grade records and start
  Introduction/Methods using existing Nyl12 text; compare newer saddle
  optimizers, and commit/push at the evening stopping point. This draft and
  ledger implement the writing/record-keeping portion; optimizer integration
  and TS searches remain future work.

## Publication readiness gaps

No accepted TI1/TS1/ES geometry, full reaction profile, completed RESP fit,
validated new FF, RFD3 redesign campaign or new sequence prediction campaign
can yet be described as completed. Before including historical quantitative
activity comparisons, recover their original assay definitions and raw tables.
Before release, complete the bibliography and archive accepted input/output
artifacts with a repository revision and retrievable storage location.
