# QM/MM work

Current path: [PA6 / Asp272-assisted relay](pa6_asp272/README.md).

Start with `pa6_asp272/plots/distances.png` or
`pa6_asp272/vmd/view_path.vmd`. The exact executed controller is exposed as
`pa6_asp272/scripts/run_path.py`; it is a provenance copy for inspection.

`01_og1_to_asp272` and `02_asp272_to_nterm` name the elementary-step hypotheses.
Within each, `window_00`, `window_01`, ... follow calculation order. Window
inputs and geometry records are retained, and `window_order.json` maps VMD
frames to those directories and the unchanged original scientific run.

The curated path folder is generated and excluded from Git; its builder,
configuration, tests and this guide are versioned. On the laptop, the same
layout is inside `nyl12-inspection-20261021/qmmm/` after the ordinary pull.
