"""Does the e-e split generalize to L > 0 multipole channels?

Each multipole kernel is K_L = r_<^L / r_>^{L+1} = phi(min) psi(max) with
phi = r^L, psi = r^{-(L+1)}.  The min/max identity gives, exactly and per L:

    K_L = 1/2 [ phi(x)psi(y) + phi(y)psi(x) ]  -  1/2 | phi(x)psi(y) - phi(y)psi(x) |
          [ separable, smooth ]                   [ W_L >= 0, all the kink ]

so the tensor head is  1/2 [ A_ik B_jl + B_ik A_jl ]  with  A = <r^L>, B = <r^-(L+1)>
--- ALWAYS rank <= 2 in the pair matricization, at every L.

What is L=0-SPECIAL (measured here): only at L=0 is the head pure labels-x-metric
(A = S, B = (k/n) diag).  At L >= 1, A and B are banded matrices (Sturmian recursion:
r is tridiagonal, so <r^L> is (2L+1)-banded) --- structured/skeleton-flavored but not
label-diagonal.

Measured per L in the radial family l = L (the lowest channel that carries the
multipole): split identity, head rank, A/B band structure, W_L spectrum decay.
"""
import io
import json
import os
import sys

import numpy as np
from scipy.special import eval_genlaguerre
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from geovac import transcorrelated_sturmian as TC  # noqa: E402


def R_nl(n, l, r, k):
    """L2-normalized Coulomb-Sturmian radial function (shared k)."""
    x = 2 * k * r
    norm = np.sqrt((2 * k) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * x ** l * np.exp(-x / 2) * eval_genlaguerre(n - l - 1, 2 * l + 1, x)


out = {}
k, Ng, nfun = 2.0, 900, 5
r, wr = TC.make_grid(k, Ng=Ng)
Wt = r * r * wr

print("=" * 84)
print("per-L split: identity | head rank | A,B band structure | W_L spectrum decay")
print("=" * 84)
print(f"{'L':>2}{'split dev':>12}{'head sv3/sv1':>14}{'A band-resid':>14}"
      f"{'B band-resid':>14}{'W tail m=4':>12}{'m=8':>8}")
for L in range(0, 4):
    l = L
    ns_list = list(range(l + 1, l + 1 + nfun))
    R = {n: R_nl(n, l, r, k) for n in ns_list}
    # radial one-body matrices for the head
    A = np.array([[np.sum(R[m] * R[n] * r ** L * Wt) for n in ns_list] for m in ns_list])
    B = np.array([[np.sum(R[m] * R[n] * r ** (-(L + 1)) * Wt) for n in ns_list]
                  for m in ns_list])
    # kernels, exact pointwise
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    lo, hi = np.minimum(R1, R2), np.maximum(R1, R2)
    K_full = lo ** L / hi ** (L + 1)
    t1 = R1 ** L / R2 ** (L + 1)
    t2 = R2 ** L / R1 ** (L + 1)
    K_sep = 0.5 * (t1 + t2)
    K_W = 0.5 * np.abs(t1 - t2)
    kernel_dev = float(np.abs(K_full - (K_sep - K_W)).max() / np.abs(K_full).max())

    D = {(i, j): R[ns_list[i]] * R[ns_list[j]] * Wt
         for i in range(nfun) for j in range(nfun)}

    def eri(K):
        g = np.zeros((nfun,) * 4)
        for i in range(nfun):
            for j in range(nfun):
                for kk in range(nfun):
                    for ll in range(nfun):
                        g[i, j, kk, ll] = D[(i, kk)] @ K @ D[(j, ll)]
        return g

    g, gs, gw = eri(K_full), eri(K_sep), eri(K_W)
    split_dev = float(np.abs(g - (gs - gw)).max() / np.abs(g).max())
    # head = 1/2 (A x B + B x A): verify + rank
    head = 0.5 * (np.einsum("ik,jl->ijkl", A, B) + np.einsum("ik,jl->ijkl", B, A))
    head_dev = float(np.abs(gs - head).max() / np.abs(gs).max())
    Pm = gs.transpose(0, 2, 1, 3).reshape(nfun * nfun, nfun * nfun)
    sv = np.linalg.svd(Pm, compute_uv=False)
    # band structure: residual outside |m-n| <= band
    def band_resid(M, band):
        mask = np.abs(np.subtract.outer(range(nfun), range(nfun))) > band
        return float(np.linalg.norm(M[mask]) / np.linalg.norm(M))
    a_res = band_resid(A, L if L > 0 else 0)          # <r^L> expected (2L+1)-banded
    b_res = band_resid(B, nfun)                        # B not claimed banded; report full
    b_res_diag = band_resid(B, 0)
    # W spectrum
    Wm = gw.transpose(0, 2, 1, 3).reshape(nfun * nfun, nfun * nfun)
    lam = np.abs(np.linalg.eigvalsh(Wm))
    lam = np.sort(lam)[::-1]
    tot = np.sum(lam ** 2)
    tail4 = float(np.sqrt(np.sum(lam[4:] ** 2) / tot))
    tail8 = float(np.sqrt(np.sum(lam[8:] ** 2) / tot))
    print(f"{L:>2}{split_dev:>12.1e}{sv[2]/sv[0]:>14.1e}{a_res:>14.1e}"
          f"{b_res_diag:>14.1e}{tail4:>12.4f}{tail8:>8.4f}")
    out[f"L={L}"] = dict(kernel_dev=kernel_dev, split_dev=split_dev,
                         head_dev=head_dev, head_sv=[float(x) for x in sv[:4]],
                         A_band_resid=a_res, B_diag_resid=b_res_diag,
                         W_tail4=tail4, W_tail8=tail8)

print()
print("A-band residual uses band = L (i.e. <r^L> claimed (2L+1)-banded);")
print("B column shows the residual off the DIAGONAL (only L=0 should be ~0).")
print("head verification max dev per L:",
      ", ".join(f"L={L}: {out[f'L={L}']['head_dev']:.1e}" for L in range(4)))

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/ee_split_L_probe.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/ee_split_L_probe.json")
