# Accepted-prefix main57 continuation to100ps

PI-authorized continuation of the same57 accepted60ps histories, eight additional5ps blocks (native stages15–22), with57 single-thread CPU workers onAzalea. No extra equilibration; accepted restart coordinates, velocities, time andbox are preserved byte-for-byte. Same topology, Hamiltonian, QM98, movable4165, total104922, SCC controls, restraints, temperature protocol andnativeCPU Sander. New disjoint per-block ig values follow the existing continuation policy; this is not a claim of bitwise thermostat RNG continuity.

`native_advancement_receipt.json` proves57nativeprocesses and first SCC-completedsteps40–47; this establishes actualprogress, not futureblockacceptance orsamplingconvergence. Historical57 last5ps nativewalltime3.997–4.028h gives approximately4h to65ps and32.2h to100ps (roughly29Nov11UTC). Update estimates fromnative runtime_summary intervalrates.

Native source andinputs are pinned byplan.json/plan.sha256 andsource_snapshot. Copiedruntime differs only bymaximum12→20 blockguard. Existingaccepted-prefixvalidation andfair minimum-block queue areunchanged. All OMP/MKL/OPENBLAS/NUMEXPR threads1; CUDAhidden. First preparationattempt preserved after terminatingonly verified duplicate SCC-replayPython; no nativeMD was stopped. Secondpreparation verifies existing exit/SCC/checkpoint receipts andallactualrestartvelocities beforelaunch.

Each block requires nativeexit0, all10000integrationsteps SCCconverged, finite energies/temperatures, finalT225–375K, unchanged frozenoutercoordinates andfinite validcheckpoint. Derivedanalysis must independently pass physicalHamiltonian/nativebias/CV parity. Processcompletion is not samplingconvergence.

Futurederivedprocessing: ../analysis/main100_continuation_20261127/continuation_registration.json. Persistentmain-onlycompact sync300s, fixed65–100ps checkpointmonitor (same200bootstrapreplicates,widths5/2.5ps, one5psdiscard perwholelineage), andseparatecurrent_ext100pool publisher waitingforqualified65ps mainfit. Existing60ps analysis/mainmovie arepreserved. Pooledpoint estimates do not inheritconfidence intervals orconvergence claims.
