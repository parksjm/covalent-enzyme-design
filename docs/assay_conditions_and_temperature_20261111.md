# Assay conditions and simulation interpretation — 2026-11-11

User-provided experimental conditions (crude lysate assays): Nyl12 variants
expressed in E. coli; PA6 or PA66 polymer powder, 5 mg per well in 96-well PCR
plates; enzyme-containing lysates in 200 mM Tris, reported pH 8.0; incubation
75 °C for 6 h. Reaction volume and buffer pH adjustment/measurement temperature
are not established here. User explicitly answered “Not sure” about the latter.
Do not infer a 200 microliter reaction volume from the buffer concentration.

The current matched neutral/charged model30 MD used 300 K, not 348.15 K
(75 °C). These simulations characterize soluble oligomer binding and selected
protonation microstates, not heterogeneous polymer-powder turnover in lysate.
Temperature affects conformational sampling, solvent and protonation equilibria;
crude-lysate conversion also depends on enzyme abundance/stability and substrate
accessibility. Low conversion does not directly establish an elementary barrier.
PES minimization has no thermostat temperature; its starting conformations and
fixed environment nevertheless inherit the preparation protocol. PMFs are
temperature-specific and must not be equated with experimental barriers.

Tris has a substantial negative pKa temperature coefficient (approximately
−0.028 per °C near ordinary reference conditions). A simple linear illustration
would lower pH by about 1.4 units on heating from 25 to 75 °C at fixed acid/base
ratio. This is NOT a measured assay pH or a validated high-temperature correction:
composition, activity, and temperature dependence require appropriate calibration.
Do not silently assign catalytic protonation from reported pH 8.0.
Reference: Brooke, Movahed and Bothner, Universal buffers for use in biochemistry
and biophysical experiments, AIMS Biophysics (2015):
https://pmc.ncbi.nlm.nih.gov/articles/PMC8956001/

Plan: retain the current 300 K campaign as a matched baseline. After finding a
credible reaction route, consider a matched 348.15 K equilibration/sampling check
with explicit protonation assumptions and the same force field. No hot MD or
buffer/protonation changes were authorized or launched by this note.

Current neutral route decision: use DFTB3/MM with the accepted 12 Å movable
region (outer atoms fixed), first drive HG1 from OG1 to the neutral Thr233
N-terminus and release the reaction bias. Endpoint N–OG1=2.679 Å and N–HG1=2.038 Å
support testing this route; N–leaving-N=3.786 Å does not prove later acid donation.
Asp272/274 OG1–OD1 distances of 5.235/4.875 Å do not currently support a direct
one-step Asp relay. Water-assisted and Asp-assisted routes remain hypotheses;
explicit QM waters and changed QM regions require separate mapped comparisons.
A coupled asynchronous proton-transfer/attack path is a bounded alternative if
the proton-transferred state is not retained. NEB follows suitable endpoints;
neither a driven peak nor a geometry gate establishes a true transition state.
