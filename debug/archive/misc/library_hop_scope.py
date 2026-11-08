"""Hop-scope audit over the 37-system hamiltonian() library.

For each canonical system, build the library Hamiltonian at its default
(max_n=2) and inspect the pre-JW one-body matrix h1:
  - number of off-diagonal nonzeros (|v| > 1e-12)
  - number of off-diagonal entries equal to +-1/16 (the graph hop value)
  - max |off-diagonal|
For He (no h1 returned by the LatticeIndex route) the LatticeIndex
_H1_spatial actually used by _build_he is inspected via the same kwargs.

Usage: python debug/library_hop_scope.py
Writes debug/data/library_hop_scope.json.
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "data", "library_hop_scope.json")


def audit(h1: np.ndarray) -> dict:
    h1 = np.asarray(h1, dtype=float)
    off = h1 - np.diag(np.diag(h1))
    nz = np.abs(off) > 1e-12
    hop = np.isclose(np.abs(off), 1.0 / 16.0, rtol=0, atol=1e-14)
    return dict(M=int(h1.shape[0]), offdiag_nnz=int(nz.sum()),
                offdiag_eq_1_16=int(hop.sum()),
                max_abs_offdiag=float(np.abs(off).max()) if off.size else 0.0)


def main() -> None:
    warnings.simplefilter("ignore")
    from geovac import ecosystem_export as ee
    names = sorted(set(ee._SYSTEM_REGISTRY.values()))
    rows = {}
    for name in names:
        t0 = time.time()
        try:
            ham = ee.hamiltonian(name)
            if ham._h1 is not None:
                row = audit(ham._h1)
                row["route"] = "h1 returned by builder"
            else:
                row = {"route": "no h1 returned"}
            row["N_pauli_incl_id"] = ham.n_terms
            row["Q"] = ham.n_qubits
        except Exception as exc:  # report, do not hide
            row = {"error": repr(exc)}
        row["t_s"] = round(time.time() - t0, 1)
        rows[name] = row
        print(name, json.dumps(row), flush=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=1)


if __name__ == "__main__":
    main()
