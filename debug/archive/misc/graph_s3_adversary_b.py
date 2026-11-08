"""Expectation of the TRUE (diagonal-h1) He Hamiltonian in the graph-native ground state.
If <psi_on|H_true|psi_on> is far above E_on, the graph off-diagonal 'correlation' is
an energy-lowering term absent from the physical Hamiltonian (not a variational bound)."""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "1"
import sys
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
import geovac.casimir_ci as cc

E_exact = -2.903724377
orig = cc._build_graph_h1
def diag_only(Z, n_max):
    h, orb = orig(Z, n_max)
    return np.diag(np.diag(h)), orb
for nm in (3, 4, 5):
    cc._build_graph_h1 = orig
    H = cc.build_graph_native_fci(2, nm)
    w, V = np.linalg.eigh(H)
    psi = V[:, 0]
    cc._build_graph_h1 = diag_only
    H0 = cc.build_graph_native_fci(2, nm)
    cc._build_graph_h1 = orig
    et = psi @ H0 @ psi
    print(f"n_max={nm} E_on={w[0]:.6f} <psi_on|H_true|psi_on>={et:.6f} "
          f"({100*(et-E_exact)/abs(E_exact):+.4f}%)  offdiag contribution={w[0]-et:+.6f}")
