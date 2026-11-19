# Prepared-system MD runner

The delivered executable datasets and configs are WT only; the legacy mutant
panel is historical evidence, not a prepared mutant MD release.

This tool runs the repository's existing protocol. It does not mutate proteins,
assign protonation, fit charges, rebuild a topology, or launch anything during
`plan`. Supply the repository checkout explicitly with `--repo` and use its
existing Python environment with NumPy and ParmEd. No package installation is
performed. An existing CUDA-capable Amber installation and tmux are required
for `launch`.

Choose `config_wt_asp272_ash274.json` (primary),
`config_wt_asp272_asp274.json` (baseline), or `config_wt_ash272_asp274.json`
(contextual control). Set the explicit environment variables `AMBERHOME` and
`PMEMD_CUDA`, or replace those fields with actual installation paths. Missing
variables fail. Review output paths, independent replica seeds and GPU indices
before planning. The supplied configs already have distinct seed ranges. Paths
for the dataset manifest and output resolve relative to that config file. Retain
a separate config/output per dataset or mutant. Each replica's six heating,
equilibration and production seeds must be distinct across replicas.

```bash
export AMBERHOME=/path/to/amber
export PMEMD_CUDA=/path/to/amber/bin/pmemd.cuda
python tools/run_prepared_md.py plan --repo /path/to/covalent-enzyme-design --config tools/config_wt_asp272_ash274.json
python tools/run_prepared_md.py launch --repo /path/to/covalent-enzyme-design --config tools/config_wt_asp272_ash274.json
```

Run these from the extracted handoff directory, using the pinned Python
environment described in `../SOFTWARE_AND_REPRODUCIBILITY.md` (relative to this
README's directory). The `--repo` path can be its supplied `project_code/`
snapshot or an explicitly reviewed checkout. `plan` creates source/input hashes,
environment metadata, exact stage-input previews and per-replica plans.
`launch` verifies those frozen hashes, checks each explicit NVIDIA GPU and
starts one recorded tmux job per replica. There is no CPU fallback. Every
stage uses the specified `pmemd.cuda`; engine or acceptance failures stop that
replica. Job command, environment, log and status are under
`runs/<dataset>/<replica>/job/`; native stage files and acceptance records are
under its `dynamics/` directory. No source plan or run is overwritten.

## Second wave: reach four independent primary and baseline replicas

The original configs remain the two-replica starter wave. For at least three,
preferably four independent replicas, use these additional configs unchanged:

- `config_wt_asp272_ash274_wave2.json`: primary replicas 03/04.
- `config_wt_asp272_asp274_wave2.json`: baseline replicas 03/04.

**Wait until both wave-1 replicas using GPUs 0/1 finish before launching the
corresponding second wave.** The second wave explicitly reuses GPUs 0/1; it
must not overlap any other wave using those devices. Nothing launches or queues
a second wave automatically. Check job/status and native acceptance records;
process completion is separate from equilibration or sampling convergence.

Each second-wave config has a distinct `_wave2` output directory and tmux
session prefix, so wave-1 plans/results remain untouched. Replicas 03/04 have
new seed ranges; all six stage seeds are distinct across all supplied dataset
and wave configs. They retain exactly the same 100 ns (10 + 90 ns) protocol.
Use the same explicit engine environment and plan/launch commands, substituting
the desired second-wave config, for example:

```bash
python tools/run_prepared_md.py plan --repo /path/to/covalent-enzyme-design --config tools/config_wt_asp272_ash274_wave2.json
python tools/run_prepared_md.py launch --repo /path/to/covalent-enzyme-design --config tools/config_wt_asp272_ash274_wave2.json
```

## Exact protocol boundary

Starts must already be accepted **unrestrained MM-minimized** systems. The
pinned `md_start_receipt.json` must certify the exact provided topology and
coordinate hashes. Raw solvated coordinates, another replica's minimized
receipt, or an unreviewed protonation/parameter change are rejected. MD
preparation acceptance does not imply a stationary point or sampling convergence.
Retain the supplied original minimization input and acceptance evidence.

The unchanged repository protocol is 100 ps NVT heating, 100 ps NPT density,
100 ps weak NPT, 200 ps unrestrained NPT equilibration, then 10 ns + 90 ns
production, using 2 fs timesteps at 300 K. Heavy-solute restraints are
5, 2 and 0.5 kcal/mol/Å² in the first three stages, referenced to the same
accepted initial minimum; equilibration and both production segments have
`ntr=0`. Production also has `nmropt=0`. The runner imports
`pilot_md_protocol.py`, `run_model30_md_replica.py` and
`continue_model30_md.py`; it does not maintain a second MD implementation.

## New mutants and microstates require explicit preparation

A changed sequence or microstate requires a separately prepared, validated
system and dataset entry. Current execution supports only ASP272/ASP274,
ASP272/ASH274 and ASH272/ASP274 with neutral Thr NH2. Other protonation
patterns require explicit method/tool review; changing a JSON label is not
sufficient. Use the repository's authoritative protein and PA6
parameter release, retaining its charge method and topology conventions.
Rebuild changed residues with an explicit source-chain/full-sequence/Amber
mapping, review termini/proton inventory, recount charge and ions, then run the
supplied unrestrained MM preparation and issue a new matching receipt. Copying
or relabeling a WT topology is not mutant preparation. The existing
`prepare_flipped_model30_ash272.py` is WT Asp272/Asp274 staging logic, not a
mutant builder; it requires dry-charge/ion and minimized-start verification.

The execution gate checks Phe113, Lys156, Asn186, Thr233, Asp272 and Asp274 in
**all four** AE/BF/CG/DH sites using 24 explicit manifest records. Beta-chain
PDB residues 1/40/42 correspond to full-sequence 233/272/274; these numbering
systems are not interchangeable. Neutral Thr NH2 and declared Asp/ASH proton
inventories are checked. For conservative design protection also retain the
maturation boundary at 232, Thr234 and Gly273 in every symmetry copy. Neither
these protections nor successful MD certify a productive catalytic pathway.
LigandMPNN design and mutation selection are separate reviewed steps.
