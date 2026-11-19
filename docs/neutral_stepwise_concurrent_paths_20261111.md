# Neutral model30: stepwise and concurrent paths

User authorized both branches on 2026-11-11. Keep them under
qmmm/pa6_flipped_model30/05_neutral_reaction_path.

02_nucleophilic_attack starts from the reaction-unrestrained, converged
01_direct_pt/02_unrestrained_endpoint/optimized.ncrst. Drive only OG1–scissile C
in approximately0.10 Å increments toward1.50 Å, then release the reactionbias.
No extra N–H restraints.

03_concurrent_pt_attack starts from inputs/accepted_es.ncrst, the original
optimized E-S state. Drive xi=d(OG1,HG1)-d(N,HG1)-d(OG1,C) with increments
no greater than0.15 Å, bounded finaltarget0.0 Å and possible early release
when TI1-like geometry is reached. This scalar coordinate allows asynchronous
component motion but is degenerate: monitor individual distances and all
proton identities. It is a diagnostic PES scan, not an approved US coordinate.
The earlier recommendation for separate coordinates in production free-energy
sampling remains in force. No inference of chemistry from xi alone.

Both branches retain DFTB3/MM,98realQM atoms+8linkH,charge−2,singlet, and the
same4165movable atoms selected for the original12 Å optimization. Use the
verified compact-residue mask, not the failed long atommask. Outer atoms remain
fixed. Track HG1/H1/H2 relative to ThrN andOG1,scissileC-N,C-O,andOG1-C.
Unintended chemistry is a diagnostic outcome, not a successfulTI1.
Release all reactionbiases at endpoints; forceconvergence/chemical retention
do not establish Hessian minimum character or a true transitionstate.

Preserve original data; retain input/source hashes and per-window energy,
distance and VMD outputs. Compare both pathways to the same originalES energy
when assembling a combined profile; label any per-branch zero explicitly.

## Direct four-center acylation branch (additional user authorization)

04_direct_acylenzyme starts from the SAME originalES. Unlike03, HG1acceptor
is substrateamideN19929, notThrN13905. Drive
xi=d(OG1,HG1)-d(HG1,substrateN)-d(OG1,C)+d(C,substrateN).
Initialxi=-3.273755808 Å; boundedfinal+2.5 gives39windows at0.148044 Å.
Singlecompositeconstraint allowsasynchronousbondchanges; four-centerroute
is the testedhypothesis, not establishedTSgeometry. NoextraN-Hrestraints.
ThrN H1/H2 and originalsubstrateHN19937 mustremain intact. HG1may be shared
betweenOG1 andsubstrateN duringtransfer; transfer toThrN is differentchemistry
and stops thisbranch. AllowC-Ncleavage andintermediateC-Oelongation.
Earlyacyllikecandidate: OG1-C<=1.60,C-N>=2.40,C-O<=1.30,substrateN-HG1<=1.25,
OG1-HG1>=1.50 Å, withprotoninventoriesretained. Removeallreactionbias for
endpointoptimization; fixedouteratoms remain as inotherpaths.
NoHessianminimum orTSvalidation inferred fromgeometry/forceconvergence.

VMDdiagnostics now measureHG1-substrateN; updatedexporter includesdynamic
HG1-substrateN cutoff1.40 Å andscissileC-N cutoff1.90 Å, removingstatic
topologybonds thatcould hidecleavage. Originalcalculationinputs unchanged.

## Third path: direct four-member acyl-enzyme route

Launched 2026-11-11 19:54 UTC in tmux
`model30_direct_acylenzyme_20261111`, using one CPU. Retained inputs, immutable
source snapshot, protocol and launcher records are in
`04_direct_acylenzyme/`; changing state is in `work/status.json` and
`job/status.json`. This path starts from the original accepted E-S restart,
not either existing driven branch. The verified compact-residue input retains
4165 movable atoms, 98 real QM atoms plus 8 link H, charge -2, DFTB3/MM.

Drive only xi=d(OG1,HG1)-d(HG1,substrateN)-d(OG1,C)+d(C,substrateN),
with one generalized restraint (weights 1,-1,-1,+1; force constant 100).
Mapping is OG1=13916, HG1=13917, C=19916, O=19917, substrateN=19929;
Thr terminal N=13905 is not the intended acceptor. Substrate original HN=19937
must remain bound; Thr H1/H2=13906/13907 must remain NH2. Initial distances
are OG1-HG1=0.981202, HG1-substrateN=2.466322, OG1-C=3.148239,
C-substrateN=1.359603, C-O=1.258474 A. Initial xi=-3.273755808 A;
39 targets advance by 0.148045 A toward +2.5 A, with early endpoint release.

Native convergence and mapped chemistry are checked after every window.
Intermediate carbonyl elongation is allowed (C-O 1.10–1.55 A), as is
scissile cleavage (C-N 1.20–4.50 A). All other QM H-parent bonds must remain
0.75–1.30 A; the terminal NH2 and original substrate NH additionally satisfy
0.85–1.25 A. HG1 may be shared between OG1 and substrate N with nearest
bond 0.75–1.70 A. ThrN-HG1 shorter than 1.30 A stops this branch as unintended
chemistry. These checks apply no additional atom restraints.

An ester/NH2 candidate requires OG1-C<=1.60, C-N>=2.40, C-O<=1.30,
HG1-substrateN<=1.25, and HG1-OG1>=1.50 A plus retained proton inventories.
Release all reaction restraints and optimize for up to 10000 baseline cycles;
record convergence and chemical retention separately. A retained candidate
is not a validated minimum or TS. The scalar scan is degenerate, diagnostic
PES work, not a PMF or validated barrier. Per-window exports retain individual
reaction distances, potential energies, HG1-substrateN labels/dynamic bonds,
and dynamic scissile C-N bonds. The reused historical internal folder names
contain `ti1`; the protocol and direct chemistry gates define the acyl endpoint.

Verification: Ruff passed; 13 focused chemistry/exporter tests passed and the
three chemistry tests passed again after the reviewed terminal-N gate update.
Coverage percentage has not been measured. Runtime inspection confirmed
Sander evaluating window 1 with the exact four-pair restraint and retained
input hashes. No other calculation was stopped or changed.

## Direct-route stop clarified from raw records

The original04path stopped at restrained asynchronous_21, not after an
unrestrainedendpointoptimization. Window20 was converged and retained
HG1-substrateN=1.213516 Å,OG1-HG1=1.459490 Å,OG1-C=2.086188 Å,C-N=1.533745 Å.
Window21 stillhad thefour-pairreactionrestraint (target−0.1648103729 Å),
butHG1-ThrN became1.052814 Å andHG1-substrateN1.954540 Å. It was included
in completed_windows forpublication, but chemistry_retained failed; it must
not be treated as anaccepted direct-routewindow. The endpointrelease directory
is empty. Thirty-ninewindows wereplanned;21triggeredthestop.

Userrequestscontinuation. Proposedlastgoodrestart iswindow20 (not21 orany
release). Rootrequestedchoice: addflatbottomHG1-substrateN restraint free
0.95–1.25 Å,k50 (outer0.8/1.5 Å) toretainintendeddirectroute, orallowthe
Nterminaldetour underoriginalcoordinate. Nochoiceassumed; noextensionlaunched
atthisnote. Anynewbiasdefinesapiecewisemodifiedscan andmustbelabeled, with
biasenergiesseparate andallreactionbiasremovedforendpointretentiontest.

## User-approved endpoint-first pivot

User allowed either the additional proton guard or building an acylenzyme
from window20, and expressed concern that a guard would be artificial. Lead
selected endpoint-first construction. The prepared guarded-continuation code
was NOT launched; original21 and alloriginalscanrecords remain unchanged.

Build an atom-preserving ester plus neutral leaving amine fromwindow20:
Thr retainsH1/H2, HG1 joins substrateN withoriginalHN19937, totalcharge−2.
Use scripted internalcoordinates and measured clash checks, preserving exact
atomorder and the original4165movable set. Optimize with DFTB3/MM without
reaction or protonholding restraints; outeratoms remainfixed as in the
matchedE-S calculation. Constructedcoordinates are a hypothesis, not an
optimizedstate. Reuse the E-S topology only for QM/MM with reactivebonds
insideQM; it is not a classicalacylenzyme forcefield/topology.

If chemistryretains and forceconverges, use accepted E-S/acylenzyme endpoints
for a pathconnection (NEBpreferred). Otherwise preserve the negative result
and reassess, rather than label or force a stableacylenzyme. Endpointretention
does not prove Hessianminimum character; a pathmaximum is not a trueTS.
