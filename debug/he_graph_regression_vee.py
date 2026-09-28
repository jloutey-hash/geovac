"""He graph-native regression (2026-09-25), stage 1: cache the V_ee config matrix.

V_ee enters the graph-native FCI at k = Z and is exactly linear in k
(two_electron_integral(..., k) = k * value at k=1), so the singlet M_L=0 V_ee
config matrix is computed ONCE per n_max at k = 1 (Z=1, one-body zeroed) and
reused for every Z:   H(Z) = Z^2 * D1 + Z * V1 + hop   (hop Z-independent).

Stage 2 (he_graph_regression.py) validates this recomposition bit-for-bit
against build_graph_native_fci and build_fci_matrix before using it.

Writes debug/data/he_graph_regression_V1_n{n}.npz
"""
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np
import geovac.casimir_ci as C

_orig = C._build_graph_h1


def _zero_h1(Z, n_max):
    h1, orb = _orig(Z, n_max)
    return np.zeros_like(h1), orb


def main() -> None:
    nmaxes = [int(a) for a in sys.argv[1:]] or [2, 3, 4, 5, 6, 7]
    for n in nmaxes:
        out = f"debug/data/he_graph_regression_V1_n{n}.npz"
        if os.path.exists(out):
            print(f"n_max={n}: cached", flush=True)
            continue
        t0 = time.time()
        C._build_graph_h1 = _zero_h1
        try:
            V1 = C.build_graph_native_fci(Z=1, n_max=n)
        finally:
            C._build_graph_h1 = _orig
        np.savez(out, V1=V1)
        print(f"n_max={n}: dim={V1.shape[0]}  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
