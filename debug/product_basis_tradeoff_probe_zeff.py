"""Addendum to product_basis_tradeoff_probe.py: best effective-charge screened Lambda(Zeff,1)
for the lowest sigma valence state of h_eff (Hartree + core projector), and cond(W)."""
import os, sys
import numpy as np
from scipy.linalg import eigh
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from prolate_separation_constant import OneE
from product_basis_tradeoff_probe import build_eff, wframe
out = []
def log(s):
    print(s, flush=True); out.append(s)
Nx, Ne, al = 30, 40, 1.2
for R in (2.6, 3.015, 3.45):
    d = build_eff(R, 0, Nx, Ne, al, True, True)
    log(f"R={R}  cond(W)={np.linalg.cond(d['W']):.1e}")
    for Z in (1.0, 1.25, 1.5, 2.0, 2.5, 3.0):
        N = OneE(R, Z, 1.0, 0, Nx, Ne, al).N
        L, (Hh, Nh) = wframe(d["W"], d["H"], N)
        E, V = eigh(Hh)
        s = []
        for i in range(3):
            v = V[:, i]; lam = v @ Nh @ v
            ev = np.linalg.eigvalsh(Nh)
            s.append(f"leak={np.linalg.norm(Nh @ v - lam * v):.2f}")
        # label purity: largest |<Lambda-eigvec|v0>|^2
        w, U = eigh(Nh); p = (U.T @ V[:, 0]) ** 2
        log(f"   Zeff={Z:4.2f}  " + "  ".join(s) + f"   max purity(v0)={p.max():.3f}")
with open(os.path.join(HERE, "data", "product_basis_tradeoff_probe_zeff.log"), "w") as f:
    f.write("\n".join(out))
