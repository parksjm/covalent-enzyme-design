"""Guarded ctypes access to the installed DL-FIND L-BFGS implementation.

The algorithm is in AmberTools/libdlfind/src/dlf_lbfgs.F90 (Nocedal 1980;
Liu and Nocedal 1989). This module supplies no optimizer or line search.
Callers own units, convergence checks, and displacement limits. The ABI is
Fortran default INTEGER by reference and explicit-shaped REAL(8) arrays.
"""

from __future__ import annotations

import ctypes
import hashlib
import threading
from pathlib import Path
from types import TracebackType

import numpy as np
from numpy.typing import ArrayLike, NDArray

_NATIVE_STATE = threading.Lock()
_DOUBLE_PTR = ctypes.POINTER(ctypes.c_double)


class DLFindLBFGS:
    """Own one exclusive native optimizer, closed explicitly or by context.

    Args:
        n_variables: Positive number of scalar coordinates.
        library_path: Explicit path to the installed DL-FIND shared library.
        source_path: Source implementing the inspected native ABI.
        memory: Positive number of L-BFGS history pairs.
    """

    def __init__(
        self,
        n_variables: int,
        library_path: Path,
        source_path: Path,
        memory: int = 5,
    ) -> None:
        self._closed = True
        for value in (n_variables, memory):
            if type(value) is not int or not 0 < value <= 2**31 - 1:
                raise ValueError("Dimensions and memory must be positive int32 values")
        self._n = n_variables
        self._owner = threading.get_ident()
        self._previous: tuple[NDArray[np.float64], NDArray[np.float64]] | None = None
        library_path, source_path = library_path.resolve(), source_path.resolve()
        self.provenance = {
            "library_path": str(library_path),
            "library_sha256": hashlib.sha256(library_path.read_bytes()).hexdigest(),
            "source_path": str(source_path),
            "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "abi": "Fortran int32-by-reference, explicit-shaped float64 arrays",
        }
        self._lib = ctypes.CDLL(str(library_path))
        signatures = {
            "init": [ctypes.POINTER(ctypes.c_int)] * 2,
            "step": [_DOUBLE_PTR] * 3,
            "restart": [],
            "destroy": [],
        }
        for name, arguments in signatures.items():
            function = getattr(self._lib, f"dlf_lbfgs_{name}_")
            function.argtypes, function.restype = arguments, None
        if not _NATIVE_STATE.acquire(blocking=False):
            raise RuntimeError("A DL-FIND L-BFGS instance already owns global state")
        try:
            nvar, nmem = ctypes.c_int(n_variables), ctypes.c_int(memory)
            self._lib.dlf_lbfgs_init_(ctypes.byref(nvar), ctypes.byref(nmem))
        except BaseException:
            _NATIVE_STATE.release()
            raise
        self._closed = False

    def _check_open(self) -> None:
        if self._closed:
            raise RuntimeError("DL-FIND L-BFGS instance is closed")
        if threading.get_ident() != self._owner:
            raise RuntimeError("DL-FIND L-BFGS must run on its owning thread")

    def step(self, coordinates: ArrayLike, gradient: ArrayLike) -> NDArray[np.float64]:
        """Return the raw native displacement without evaluating energy.

        Args:
            coordinates: Finite flat coordinate vector in caller-selected units.
            gradient: Finite flat energy gradient in corresponding units.

        Returns:
            Independent flat displacement array; caller must apply its step cap.

        Raises:
            ValueError: Inputs are invalid, stationary, or have unsafe curvature.
            RuntimeError: Instance is closed, foreign-thread, or output is invalid.
        """
        self._check_open()
        arrays = [
            np.array(a, dtype=np.float64, order="C", copy=True)
            for a in (coordinates, gradient)
        ]
        x, g = arrays
        if any(a.shape != (self._n,) or not np.isfinite(a).all() for a in arrays):
            raise ValueError("Coordinates and gradient must be finite flat vectors")
        norm = float(np.linalg.norm(g))
        if not np.isfinite(norm) or norm == 0.0:
            raise ValueError("Gradient norm must be finite and nonzero")
        if self._previous is not None:
            old_x, old_g = self._previous
            curvature = float(np.dot(x - old_x, g - old_g))
            if not np.isfinite(curvature) or curvature <= 0.0:
                raise ValueError("Nonpositive secant curvature; explicitly restart")
        result = np.empty(self._n, dtype=np.float64)
        self._lib.dlf_lbfgs_step_(
            *(a.ctypes.data_as(_DOUBLE_PTR) for a in (x, g, result))
        )
        self._previous = (x, g)
        if not np.isfinite(result).all():
            raise RuntimeError("DL-FIND returned a nonfinite step; explicitly restart")
        return result

    def restart(self) -> None:
        """Explicitly clear native history before requesting a fresh step."""
        self._check_open()
        self._lib.dlf_lbfgs_restart_()
        self._previous = None

    def close(self) -> None:
        """Release native arrays and global ownership; repeated calls are safe."""
        if not self._closed:
            self._check_open()
            self._lib.dlf_lbfgs_destroy_()
            self._closed = True
            _NATIVE_STATE.release()

    def __enter__(self) -> DLFindLBFGS:
        self._check_open()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()
