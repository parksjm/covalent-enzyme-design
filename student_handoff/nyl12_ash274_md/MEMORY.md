# Scientific memory — student MD arm

- Goal: test whether ASP272/ASH274 improves persistent near-reactive PA6 L4
  binding in WT and prior mutants; only then explore sequence-design ranking.
- The hypothesis is unproven. Two WT ASH274 replicas yielded 1,677 and zero
  Tier A frames; ASP/ASP gave 37 and 238, each out of 10,000 production frames.
- ASH272/ASP274 gave 3,886 and 1,724. Keep it as context/optional control;
  do not silently substitute ASH272 for the requested ASH274 state.
- Tier A plus positive phi/native-Thr/intact-substrate counts are only
  5+0 for ASH274, 0+21 for ASP/ASP, and 7+0 for ASH272. These differ from the
  wider Tier A/B joint counts. Higher Tier A alone is not complete competence.
- All six matched WT trajectories used 300 K, the same flipped Model 30 and
  modular PA6 RESP/GAFF2 release. The old WT+eight-mutant panel used legacy
  AM1-BCC/GAFF treatment and is a historical comparator, not a matched control.
- Prior mutants: T53Q, F55N, K59D, K59E, F101E, F101H, P105K and V111C.
- Microstates: ASP272/ASP274 versus ASP272/ASH274 at **all four copies**;
  keep mature Thr233 N termini neutral NH2, Lys156 in its prepared state,
  substrate termini/charge and all other settings fixed. Rebuild ions/charges.
- Use ff19SB/TIP3P and the exact released P6N–P6I–P6I–P6C zwitterionic L4;
  no TI1 or acyl-enzyme classical parameters are needed for this MD arm.
- ASH274 is a model hypothesis. Fixed-charge MD cannot move its proton or
  determine a pKa. Negoro's Asp protonation assignment is circumstantial;
  maturation and substrate hydrolysis must remain distinct.
- Positive signed projected Herschlag phi denotes the flipped local pose
  under the frozen current atom convention. Never replace it by absolute phi
  when assigning orientation; retain both signed and historical Tier counts.
- Four mature alpha/beta pairs are AE/BF/CG/DH. Precursor Thr233 is beta-local
  residue 1 and Amber929 in site AE. Never reuse raw atom indices across mutants.
- Protect F113/K156/N186/T233/D272/D274 in every symmetry copy for design;
  conservatively preserve cleavage junction 232–233, Thr234 and Gly273 too.
- MD is a near-reactive-pose screen, not a rate/barrier predictor. No stable
  TI1 or verified realistic TS has been established by the current QM/MM work.
- Record every action, exact source/version/hash, microstate, seeds, restraints,
  process completion and scientific acceptance. Never silently fix failed input.

## pKa documentation addition — 2026-11-19

Read `docs/propka_asp_pair_pilot_20261116.md` and the manuscript pKa sections
in the current project checkout. Corrected mature NylC, deposited uncleaved
mutant controls and mature Nyl12 were screened with PROPKA; none determines
a unique experimental Asp-pair/Thr-N microstate. MCCE results are rejected
or input-preflight failures, not accepted independent pKas. Compact tables
and rejection receipts are in `evidence/pka/`. The original data archives
remain immutable and predate this addition; use current Git documents/evidence.
