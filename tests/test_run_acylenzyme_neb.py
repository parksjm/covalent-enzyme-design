"""Verify chemical gates and NEB/restraint input separation."""

from pathlib import Path
import pytest
from scripts.run_acylenzyme_neb import acyl_retained, input_text, restraints

BASE = "&cntrl nmropt=0,ntr=0,maxcyc=10000,ntwx=25,/\n&dlfind maxcycle=10000,maxstep=0.5,verbosity=1,/\n"


def test_neb() -> None:
    text = input_text(BASE, Path("ref"), images=[Path("product")])
    assert "neb='FREE'" in text and "neb_nimage=17" in text
    assert "coords2(1)='product'" in text and "refcrd='ref'" in text
    assert "DISANG" not in text and "nmropt=0" in text
    assert "maxcycle=1000" in text


def test_release() -> None:
    assert "DISANG" in input_text(BASE, Path("ref"), restrained=True)
    assert "DISANG" not in input_text(BASE, Path("ref"))
    assert "iat=19929,13917" in restraints()


@pytest.mark.parametrize("images", [[], [Path("p")] * 13])
def test_invalid_images(images: list[Path]) -> None:
    with pytest.raises(ValueError):
        input_text(BASE, Path("ref"), images=images)


def test_chemistry() -> None:
    g = dict(
        OG1_C=1.43,
        C_O=1.23,
        C_Nleave=2.8,
        Nleave_HG1=1.03,
        N_H1=1.02,
        N_H2=1.02,
        N_HG1=2.0,
        OG1_HG1=2.0,
    )
    assert acyl_retained(g)
    assert not acyl_retained({**g, "OG1_C": 2.9})
    assert not acyl_retained({**g, "N_HG1": 1.03})
    assert not acyl_retained({**g, "C_Nleave": 1.4})


def test_queue_and_process_failures(tmp_path: Path) -> None:
    import json
    import numpy as np
    from parmed.amber import Rst7
    from scripts.run_acylenzyme_neb import run, run_engine, geometry

    xyz = np.zeros((19930, 3))
    xyz[13915] = [0, 0, 0]
    xyz[19915] = [1.43, 0, 0]
    xyz[19916] = [1.43, 1.23, 0]
    xyz[19928] = [4.23, 0, 0]
    xyz[13916] = [5.26, 0, 0]
    xyz[13904] = [0, 4, 0]
    xyz[13905] = [1.02, 4, 0]
    xyz[13906] = [0, 5.02, 0]
    restart = Rst7(natom=len(xyz))
    restart.coordinates = xyz
    seed = tmp_path / "seed.ncrst"
    restart.write(str(seed), netcdf=True)
    engine = tmp_path / "fake_sander"
    engine.write_text(
        '#!/usr/bin/env python3\nimport sys,shutil\na=sys.argv\nshutil.copyfile(a[a.index("-c")+1],a[a.index("-r")+1])\nopen(a[a.index("-o")+1],"w").write("Optimisation converged\\n")\n'
    )
    engine.chmod(0o700)
    base = tmp_path / "base.mdin"
    base.write_text(BASE)
    config = dict(
        base_input=str(base),
        reactant=str(seed),
        seed=str(seed),
        sander=str(engine),
        topology=str(tmp_path / "topology"),
        intermediates=[],
    )
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(config))
    run(cfg)
    assert json.loads((tmp_path / "status.json").read_text())["phase"] == "finished"
    assert (tmp_path / "03_neb/process.json").exists()
    assert geometry(seed)["Nleave_HG1"] == pytest.approx(1.03)
    engine.write_text(
        '#!/usr/bin/env python3\nimport sys\na=sys.argv\nopen(a[a.index("-o")+1],"w").write("SANDER BOMB")\nsys.exit(1)\n'
    )
    with pytest.raises(RuntimeError, match="Sander failed"):
        run_engine(config, tmp_path / "failure", seed, BASE)
    bad = tmp_path / "bad"
    bad.mkdir()
    badcfg = bad / "config.json"
    badcfg.write_text(json.dumps(config))
    with pytest.raises(RuntimeError):
        run(badcfg)
    assert json.loads((bad / "status.json").read_text())["phase"] == "stopped"


def test_invalid_source_bias() -> None:
    with pytest.raises(ValueError):
        input_text(BASE.replace("nmropt=0", "nmropt=1"), Path("ref"))
