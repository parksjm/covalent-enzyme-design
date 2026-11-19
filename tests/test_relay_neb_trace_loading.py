"""Validate source trace model/exterior/end-label gates with synthetic restarts."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import relay_neb_trace as trace


@pytest.mark.parametrize(
    "bad",
    [
        None,
        "links",
        "active",
        "water",
        "shape",
        "box",
        "frozen",
        "mdin",
        "label",
        "endpoint",
    ],
)
def test_trace_model_and_frame_validation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bad: str | None
) -> None:
    size = 20400
    active = np.zeros(size, dtype=bool)
    active[:4657] = True
    active[[13915, 13916, 19932, 20316, 20317]] = True
    if bad == "active":
        active[0] = False
    source = np.zeros((size, 3))
    paths = [
        tmp_path / "03_reverse_path/unbiased_ES_release/optimized.ncrst",
        tmp_path / "02_acyl/02_unbiased_release/optimized.ncrst",
    ]
    if bad == "endpoint":
        paths[-1] = paths[-1].with_name("accepted_00000.ncrst")
    for path in paths:
        path.parent.mkdir(parents=True)
        path.write_text("synthetic accepted coordinates")
        (path.parent / "one_cycle.mdin").write_text("fixed native input")
    if bad == "mdin":
        (paths[-1].parent / "one_cycle.mdin").write_text("changed method")
    top = tmp_path / "topology.prmtop"
    top.write_text("synthetic topology identity")
    case = tmp_path / "case.json"
    case.write_text("synthetic fixed case identity")
    atoms = [SimpleNamespace()] * size
    atoms[20316] = SimpleNamespace(
        residue=SimpleNamespace(idx=1575 if bad != "water" else 1, name="WAT")
    )
    bonds = [
        SimpleNamespace(
            atom1=SimpleNamespace(idx=13915), atom2=SimpleNamespace(idx=6000)
        )
        for _ in range(8 if bad != "links" else 7)
    ]
    box = np.array([40, 40, 40, 90, 90, 90])
    model = SimpleNamespace(
        qm=np.array([13915]),
        active=active,
        coordinates=source,
        topology=SimpleNamespace(atoms=atoms, bonds=bonds),
        restart=SimpleNamespace(box=box),
        atom_map=dict(water_o=20317, proton=13917, thr_og=13916, leaving_n=19933),
        config={"topology": str(top), "mdin": "fixed native input"},
    )
    monkeypatch.setattr(trace, "load_case", lambda _: model)
    monkeypatch.setattr(trace, "stage_trace_paths", lambda _: paths)

    def restart(path: str) -> SimpleNamespace:
        index = paths.index(Path(path))
        xyz = source.copy()
        xyz[0, 0] = index + 1
        if bad == "shape":
            xyz = xyz[:-1]
        elif bad == "frozen":
            xyz[6000, 0] = 1
        return SimpleNamespace(coordinates=xyz, box=box if bad != "box" else None)

    monkeypatch.setattr(trace.Rst7, "open", restart)

    def audit(_: object, xyz: np.ndarray) -> dict:
        parent = 13916 if xyz[0, 0] == 1 else 19933
        if bad == "label":
            parent = 1
        return {"hydrogen_parents_1based": {"20318": parent, "13917": 20317}}

    monkeypatch.setattr(trace, "proton_audit", audit)
    if bad:
        with pytest.raises(ValueError):
            trace.load_trace(case)
    else:
        result = trace.load_trace(case)
        assert result.active_coordinates.shape == (2, 4662, 3)
        assert result.paths[-1] == paths[-1]
        assert len(result.sha256[str(paths[-1])]) == 64
