"""Addendum 2: overlap of the lowest sigma valence state of h_eff with the JOINT (E,Lambda)
eigenstates of the separable h(Zeff,1) -- the label-purity that matters (Lambda-only
eigvectors can be near-degenerate). Also the Lambda(Zeff,1) spacing near the state."""
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
for m in (0, 1):
    for R in (2.6, 3.015, 3.45):
        d = build_eff(R, m, Nx, Ne, al, True, m == 0)
        for Z in (1.0, 1.5, 2.0, 2.5):
            op = OneE(R, Z, 1.0, m, Nx, Ne, al)
            L, (Hh, Hs, Nh) = wframe(d["W"], d["H"], op.Hm, op.N)
            E, V = eigh(Hh); Es, Vs = eigh(Hs)
            lamj = np.einsum('ij,ik,kj->j', Vs[:, :12], Nh, Vs[:, :12])
            row = []
            for i in range(2):
                p = (Vs[:, :60].T @ V[:, i]) ** 2
                k = int(p.argmax())
                row.append(f"st{i}: E={E[i]:+.4f} best={p.max():.3f} (joint #{k}, lam={lamj[k] if k < 12 else float('nan'):+.3f}) top3={np.sort(p)[-3:].sum():.3f}")
            log(f"m={m} R={R} Zeff={Z:3.1f}  " + " | ".join(row))
with open(os.path.join(HERE, "data", "product_basis_tradeoff_probe_zeff2.log"), "w") as f:
    f.write("\n".join(out))
