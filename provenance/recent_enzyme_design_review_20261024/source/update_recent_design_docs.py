"""Add verified reference identifiers and a bounded Nyl12 design strategy."""
from pathlib import Path

ROOT=Path('/data/jpc/covalent-enzyme-design')
refs=[
(23,'Wei et al.','Fundamental reaction pathway for peptide metabolism by proteasome: insights from first-principles quantum mechanical/molecular mechanical free energy calculations','J. Phys. Chem. B','2013','10.1021/jp405337v','https://pubmed.ncbi.nlm.nih.gov/24111489/'),
(24,'Chilov, Sidorova and Svedas','Quantum chemical studies of the catalytic mechanism of N-terminal nucleophile hydrolase','Biochemistry (Moscow)','2007','10.1134/S0006297907050057','https://pubmed.ncbi.nlm.nih.gov/17573703/'),
(25,'Negoro et al.','X-ray crystallographic and mutational analysis of the NylC precursor: catalytic mechanism of autocleavage and substrate hydrolysis of nylon hydrolase','FEBS Journal','2023','10.1111/febs.16755','https://febs.onlinelibrary.wiley.com/doi/10.1111/febs.16755'),
(26,'','A unified mechanism for proteolysis and autocatalytic activation in the 20S proteasome','Nature Communications','2016','10.1038/ncomms10900','https://pmc.ncbi.nlm.nih.gov/articles/PMC4792962/'),
(27,'Lauko et al.','Computational design of serine hydrolases','Science','2025','10.1126/science.adu2454','https://pmc.ncbi.nlm.nih.gov/articles/PMC12288761/'),
(28,'Kim et al.','Computational design of metallohydrolases','Nature','2026','10.1038/s41586-025-09746-w','https://www.nature.com/articles/s41586-025-09746-w'),
(29,'Choi et al.','De novo design of cysteine proteases','Nature','2026','10.1038/s41586-026-11143-w','https://www.nature.com/articles/s41586-026-11143-w'),
(30,'Chen et al.','Computational design of metalloproteases','Nature','2026','10.1038/s41586-026-11167-2','https://www.nature.com/articles/s41586-026-11167-2'),
(31,'Butcher et al.','Atomic functional site scaffolding in biopolymers with RFdiffusion3','Science','2026','10.1126/science.aeh2982','https://pubmed.ncbi.nlm.nih.gov/42879533/'),
(32,'Wang et al.','Computational loop reconstruction based design of efficient PET hydrolases','Communications Biology','2025','10.1038/s42003-025-08364-6','https://pmc.ncbi.nlm.nih.gov/articles/PMC12174350/'),
(33,'Listov et al.','Complete computational design of high-efficiency Kemp elimination enzymes','Nature','2025','10.1038/s41586-025-09136-2','https://www.nature.com/articles/s41586-025-09136-2'),
]
p=ROOT/'manuscript/REFERENCES.md'
text=p.read_text()
assert '**R23.' not in text
text += '\n## Mechanistic alternatives and recent enzyme design — added 2026-10-24\n\n'
for key,authors,title,journal,year,doi,url in refs:
    text += f'- **R{key}.** {authors+". " if authors else ""}*{title}*. {journal} ({year}).\n  https://doi.org/{doi}\n  Primary access record: {url}\n'
text += '''
R23/R24 support alternative relay/direct-transfer models in other enzymes,
not an assigned Nyl12 mechanism. R25 separates precursor autocleavage from
mature substrate hydrolysis. R26 discusses competing proteasome general-base
models. Full author lists and publication formatting remain editorial tasks.
R27 is the final Science article, not the 2024 ChemNet preprint. R28 appeared
online 2025-12-03 and in the January 2026 issue. R29/R30 appeared online
2026-10-14; R31's primary PubMed record gives online 2026-10-15. R30 detailed
motif provenance was checked in the authors' preprint, not final SI. Final
R31 full Methods/SI were inaccessible; final abstract and official code/docs
are distinguished from the older preprint. Do not transplant preprint metrics.
Official implementation sources: https://github.com/hojaec/RFD2_MI_cysteine_protease_pipeline
and https://rosettacommons.github.io/foundry/models/rfd3/tutorials/advanced_enzyme_design_tutorial.html .
Additional preserved relay structural reference (full bibliographic details
pending): *Crystal structures of penicillin acylase enzyme-substrate complexes:
structural insights into the catalytic mechanism* (2001), publisher record
https://www.sciencedirect.com/science/article/pii/S0022283601950430 .
'''
text=text.replace('  Add the final RFD3 primary-paper citation after reconciling the project\'s\n  stored Science 2026 citation with the exact release used for design.','  See R31 for the verified final 2026 primary-paper record. Any production\n  design must additionally pin its actual software release/commit.')
p.write_text(text)
bib=ROOT/'manuscript/recent_design_and_relay_references.bib'
assert not bib.exists()
bib.write_text('% Curated DOI records, accessed2026-10-24. Complete author lists before submission.\n'+ '\n'.join(f'@article{{R{k},\n  title = {{{title}}},\n  journal = {{{journal}}},\n  year = {{{year}}},\n  doi = {{{doi}}},\n  url = {{{url}}},\n  note = {{Working reference; verified source/access limits in REFERENCES.md}}\n}}\n' for k,_,title,journal,year,doi,url in refs))
p=ROOT/'manuscript/INTRODUCTION_DRAFT.md'
p.write_text(p.read_text()+'''
Mechanistic models for N-terminal nucleophile hydrolases do not establish a
universal requirement for a water bridge: water-assisted and direct-transfer
models have been reported in other family members, and a carboxylate-containing
network has been proposed for NylC [R23–R26]. We therefore treat relay geometry
and protonation as hypotheses requiring structural and chemical validation.
An unsuccessful placement search in one fixed conformation does not exclude
transient solvation or establish an alternative mechanism.

Recent generative enzyme-design studies motivate complementing local chemical
models with structural assessment across reaction states [R27–R31]. These
approaches do not convert a predicted intermediate or transition-state analogue
into a verified electronic saddle point. Native-fold and polymer-hydrolase
engineering provide additional precedents for focused redesign and conformer-
aware evaluation [R32,R33]. For Nyl12, this distinction is essential: local
amide chemistry must be organized while long substrate segments retain access,
conformational diversity and product-specific cleavage registers. We propose
separating economical catalytic-state screening from detailed reference and
shortlist reaction-path calculations, while retaining the maturation and
oligomeric architecture of the existing enzyme.
''')
p=ROOT/'manuscript/METHODS_DRAFT.md'
p.write_text(p.read_text()+'''
## Relay-water feasibility and staged design scope — evidence and proposal, 2026-10-24

A deterministic front-side oxygen search was performed in the original full
alpha/beta protein/substrate context. The user-specified oxygen defined the
front hemisphere in a camera-independent frame perpendicular to the terminal
N-to-OG1 axis. At 0.025 A distance and 0.5-degree azimuth sampling, none of
225,000 positions at 2.6–3.2 A from both atoms passed the declared steric
screen; heavy-only and opposite-side diagnostics also failed. These finite
fixed-conformation screens concern spatial compatibility, not activation
energies or mechanistic exclusion. Thresholds, atom maps, source hashes,
sensitivity analyses and failed placements are retained in the source ledger.

The subsequent design strategy is proposed, not executed. We will prioritize
a small set of chemically justified local-state models and explicit alternative
protonation/relay hypotheses, using the native structure and mapped catalytic
residues as constraints. Exploratory motif screening may use an intermediate
or transition-state-like template labelled with its actual provenance. A claim
of TS1 optimization or a quantitative barrier comparison will still require
the saddle, reaction-mode and endpoint tests specified above.

Multistate structural evaluation will cover independent noncovalent ES, TI1,
acyl enzyme and deacylation-compatible states, with matched WT controls and
apo/product access checks [R27–R33]. An initial focused LigandMPNN arm will
preserve WT mature lengths and explicit alpha/beta symmetry ties; local RFD3
backbone sampling is an optional arm rather than a replacement of the native
fold. Repeated structural predictions will be interpreted as model ensembles,
not equilibrium populations. For large polyamides, local QM truncation will
be separated from full oligomer specificity testing: distal conformers,
cleavage registers, chain ends and open-pocket geometry will be varied.
Neither activated-ester activity in a published design nor a predicted
polymer-bound pose establishes nylon turnover, L2 selectivity or processivity.
''')
p=ROOT/'manuscript/SOURCE_LEDGER.md'
p.write_text(p.read_text()+'''
## Relay references and latest design-method review — 2026-10-24

Primary relay source summaries: ../docs/water_relay_mechanism_evidence_20261024.md.
Finite packing evidence and exact screening limits:
../docs/front_water_feasibility_20261024.md. Recent-paper comparisons and
project-specific proposals: ../docs/recent_enzyme_design_lessons_20261024.md.
Stable reference keys R23–R33 and DOI-only BibTeX records are now retained;
complete author lists remain pending editorial formatting. No new QM or design
campaign was launched by this review. Do not turn proposals into executed methods.

Retrieval distinctions: Science2025 final PMC12288761 vs preprintPMC11384011;
Nature metallohydrolase onlineDec2025 vs issue2026; Oct2026 cysteine/protease
final abstract/extended data and official code vs metalloprotease detailed
motif provenance verified only in primary preprint. Final RFD3 bibliographic
record/abstract and official tutorial verified; final full Methods/SI could
not be retrieved. RFD3 hydrogen stripping means heavy-atom donor/acceptor
conditioning is not a calculation of proton locations or reaction barriers.
Primary Nature source for Listov2025 was indexed/searchable; direct full-page
retrieval redirected to authentication. PET2025 full Methods were available.

Approval-review timeouts occurred before two root/agent read commands; no
scientific input changed. A competing-literature subagent terminated with a
model content-access error; its partial suggestions were not treated as verified
method facts. Root independently checked benign native-fold/PET primary sources.
''')
