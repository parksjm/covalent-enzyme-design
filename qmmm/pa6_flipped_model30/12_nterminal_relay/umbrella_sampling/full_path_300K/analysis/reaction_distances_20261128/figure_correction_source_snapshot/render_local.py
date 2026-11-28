"""Regenerate only local figure from verified retained-row artifact."""

import csv
from pathlib import Path
from plot_distances import plot_local

WORK = Path(__file__).resolve().parent
with (WORK / "retained_rows.csv").open() as handle:
    rows = list(csv.DictReader(handle))
for row in rows:
    for key in (
        "production_time_ps",
        "r_og1_hg1",
        "r_nleave_hg1",
        "r_nalpha_hg1",
        "r_c_og1",
        "r_c_n",
        "q_del_A",
    ):
        row[key] = float(row[key])
plot_local(WORK, rows)
