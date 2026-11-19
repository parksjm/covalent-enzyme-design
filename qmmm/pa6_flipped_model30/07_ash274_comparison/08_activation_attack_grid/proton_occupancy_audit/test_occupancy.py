"""Focused tests for proton occupancy invariance."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import numpy as np

SCRIPT = Path(__file__).with_name("analyze_occupancy.py")
sys.modules.setdefault("parmed", ModuleType("parmed"))
matplotlib = ModuleType("matplotlib")
matplotlib.__path__ = []
sys.modules.setdefault("matplotlib", matplotlib)
sys.modules.setdefault("matplotlib.pyplot", ModuleType("matplotlib.pyplot"))
SPEC = importlib.util.spec_from_file_location("occupancy_analysis", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Could not load occupancy analysis module")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_coordination_is_invariant_to_hydrogen_permutation() -> None:
    """Verify atom-order changes preserve smooth coordination values."""
    hydrogens = np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [0.0, 2.0, 0.0]])
    heavy_atoms = np.array([[0.1, 0.0, 0.0], [3.0, 0.0, 0.0]])
    reference = MODULE.occupancy(hydrogens, heavy_atoms)
    permutation = MODULE.occupancy(hydrogens[[2, 0, 1]], heavy_atoms)
    assert np.allclose(reference, permutation, atol=1e-14)
    MODULE.check_permutation_invariance()
