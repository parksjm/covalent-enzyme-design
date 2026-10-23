# Focused PA6 redesign preparation

The current inputs and manifest are in `protected_inputs/`. They prepare six
LigandMPNN jobs (18 proposed sequences before deduplication) but launch none.
The original `backbones/` and root `preparation_manifest.json` record the first
preparation snapshot; `pilot_inputs/` is also superseded. Both old manifests
are marked `superseded_incomplete_protection_do_not_execute` and are retained
for provenance. The corrected mask excludes the historical maturation core
Glu229/Asn232/Thr233/Asp272/Asp274, Phe113/Asn186, and alpha 222-232, in all
copies. The current config reports mutation count rather than imposing the
provisional four-position post-design cap in the initial snapshot.

Protocol: `manuscript_nyl12_wt_vs_mutants/NYL12_FOCUSED_REDESIGN_WORKFLOW_20261016.md`.
Review backbone diversity, uncertain loops, and the candidate design mask
before executing commands. No covalent parameterization or MD has been run.
Tests and preparation run in tmux via
`scripts/launch_nyl12_focused_redesign_preflight_20261016.sh`; durable log:
`logs/preflight.log`.
