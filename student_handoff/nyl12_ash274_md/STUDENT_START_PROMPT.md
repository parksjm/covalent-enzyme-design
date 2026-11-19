# Prompt for the student's first agent session

Read `AGENTS.md`, `START_HERE.md`, `MEMORY.md`, `STUDY_PLAN.md` and
`STUDY_SUMMARY.md` in this handoff. Our task is the **noncovalent E–S MD
comparison**, not a new reaction-path campaign. Use the bundled `project_code/`
snapshot or a verified checkout of the private project repository.

1. Verify `CONTENTS.sha256`, the PA6 release, delivered WT systems, neutral
   mature Thr233 termini and all four Asp-pair states. Reproduce existing WT
   occupancy counts from cached tables. Inventory available CPUs/GPUs and
   Amber software before selecting any simulation hardware.
2. Create the locked analysis environment described in
   `SOFTWARE_AND_REPRODUCIBILITY.md`. Run the prepared-WT `plan` command and
   inspect its topology/charge/mapping checks and full MD input previews.
   Do not silently substitute force fields, water models or protonation states.
3. Present the exact matched WT ASH274 versus ASP274 schedule, hardware,
   replica seeds and output paths to the student/PI before launching the
   campaign. Use at least three, preferably four, independent 100 ns replicas
   per condition, at 300 K. The primary comparison keeps Asp272 deprotonated
   and catalytic Thr233 N termini neutral in all four beta chains. Optional
   ASH272 is a separately labeled control. Do not share busy GPUs or launch
   every example configuration simultaneously.
4. Add and validate four-site-aware catalytic occupancy alongside the existing
   site-AE metric before design ranking. Preserve signed phi, native Thr
   geometry, substrate integrity, binding/site identity, Asp hydrogen bonding,
   histograms/time series and replica-level uncertainty.
5. Rebuild the eight historical mutant pose references with the same current
   RESP/GAFF2 PA6 release and matched protonation controls. Those supplied
   mutant PDBs are not ASH274-prepared systems. Retain genotype/model-quality
   metadata. Validate each new topology/start pair and counterions before MD.
6. Only if the matched controls support a reproducible benefit, propose focused
   LigandMPNN candidates with the catalytic/protected-residue and symmetry
   maps in `STUDY_PLAN.md`. Treat occupancy as a structural proxy, not proof
   of enhanced catalysis.

Keep an append-only local action journal with intent/outcome, source hashes,
explicit software environments and rollback guidance. Use persistent scripted
tmux jobs and update concise manuscript Methods/Results and a task dashboard.
Do not delete old data, publish unpublished structures, invite collaborators,
modify workstation configuration or launch QM/MM jobs without task-specific
authorization. Ask when chemistry or a required input is ambiguous.
