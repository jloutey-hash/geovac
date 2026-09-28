"""Independent check (PM, 2026-09-25): does the graph-native He CI's 0.2% come from the
(1/16)*A one-body hop, which is not part of He's physical Hamiltonian in the k=Z
hydrogenic basis (where T + V_ne is exactly diagonal, -Z^2/2n^2)?

Reports, per n_max: E with hop (the published method), E with hop zeroed (the physical
Hamiltonian in the same basis), and <psi_hop | H_phys | psi_hop>.
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import geovac.casimir_ci as C

E_EXACT = -2.903724377  # He nonrelativistic, infinite nuclear mass

_orig = C._build_graph_h1


def _diag_only(Z: int, n_max: int):
    h1, orb = _orig(Z, n_max)
    return np.diag(np.diag(h1)), orb


def ground(H: np.ndarray):
    w, v = np.linalg.eigh(H)
    return w[0], v[:, 0]


for n_max in (2, 3, 4, 5):
    H_hop = C.build_graph_native_fci(Z=2, n_max=n_max)
    C._build_graph_h1 = _diag_only
    try:
        H_phys = C.build_graph_native_fci(Z=2, n_max=n_max)
    finally:
        C._build_graph_h1 = _orig
    E_hop, psi = ground(H_hop)
    E_phys, _ = ground(H_phys)
    E_mix = float(psi @ H_phys @ psi)
    pct = lambda e: 100.0 * (e - E_EXACT) / abs(E_EXACT)
    print(f"n_max={n_max}  dim={H_hop.shape[0]:5d}  "
          f"E_hop={E_hop:.6f} ({pct(E_hop):+.3f}%)  "
          f"E_phys={E_phys:.6f} ({pct(E_phys):+.3f}%)  "
          f"<psi_hop|H_phys|psi_hop>={E_mix:.6f} ({pct(E_mix):+.3f}%)  "
          f"hop contribution={1e3*(E_hop-E_phys):+.1f} mHa", flush=True)
