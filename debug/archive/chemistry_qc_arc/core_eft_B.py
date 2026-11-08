"""CE step B: LiH+ with the exact frozen core against the effective operators.

Record: docs/predictions/2026-10-20_core_eft.md (committed 97901bbd before computing).

Exchange is applied exactly through scaled running sums (the kernel r_<^l / r_>^(l+1) is
semiseparable), as auxiliary sparse unknowns in an augmented shift-invert system.  Orthogonality to
the core is a level shift Lambda |c><c| (one auxiliary unknown).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).parent))
import molfns_solver as S  # noqa: E402
import core_eft_A as CA  # noqa: E402

OUT = Path(__file__).parent / "data" / "core_eft_B.json"


def exchange_aug(g: S.Grid, wc: np.ndarray, ls: List[int], lK: int, lam: float):
    """Blocks of the augmented system for X = -K + lam B|c><c|B (t-major w ordering).
    Returns (Cw, Rw, D): w-rows += Cw @ aux;  Rw @ w + D @ aux = 0."""
    n, h, nc = g.n, g.h, len(ls)
    r, e2t = g.r, np.exp(2 * g.t)
    f = e2t * wc
    Crow, Ccol, Cval = [], [], []
    Rrow, Rcol, Rval = [], [], []
    Drow, Dcol, Dval = [], [], []
    na = 0
    for c, l in enumerate(ls):
        if l > lK:
            continue
        coef = 1.0 / (2 * l + 1)
        P0, Q0 = na, na + n
        na += 2 * n
        idx = np.arange(n)
        wi = idx * nc + c
        # w-rows:  -K w = -coef f_i (P_i + Q_i)
        Crow += list(wi) * 2; Ccol += list(P0 + idx) + list(Q0 + idx)
        Cval += list(-coef * f) * 2
        # P_i - e^{-h(l+1)} P_{i-1} - (h f_i / r_i) w_i = 0
        Drow += list(P0 + idx); Dcol += list(P0 + idx); Dval += [1.0] * n
        Drow += list(P0 + idx[1:]); Dcol += list(P0 + idx[:-1]); Dval += [-np.exp(-h * (l + 1))] * (n - 1)
        Rrow += list(P0 + idx); Rcol += list(wi); Rval += list(-h * f / r)
        # Q_i - e^{-h l} Q_{i+1} - h e^{-h l} (f_{i+1}/r_{i+1}) w_{i+1} = 0
        Drow += list(Q0 + idx); Dcol += list(Q0 + idx); Dval += [1.0] * n
        Drow += list(Q0 + idx[:-1]); Dcol += list(Q0 + idx[1:]); Dval += [-np.exp(-h * l)] * (n - 1)
        Rrow += list(Q0 + idx[:-1]); Rcol += list(wi[1:]); Rval += list(-h * np.exp(-h * l) * f[1:] / r[1:])
    # level shift: z - h sum e^{2t} wc w_s = 0 ; w_s rows += lam e^{2t} wc z
    if 0 in ls:
        c = ls.index(0)
        z = na; na += 1
        wi = np.arange(n) * nc + c
        Crow += list(wi); Ccol += [z] * n; Cval += list(lam * e2t * wc)
        Rrow += [z] * n; Rcol += list(wi); Rval += list(-h * e2t * wc)
        Drow += [z]; Dcol += [z]; Dval += [1.0]
    N = n * nc
    Cw = sps.csr_matrix((Cval, (Crow, Ccol)), shape=(N, na))
    Rw = sps.csr_matrix((Rval, (Rrow, Rcol)), shape=(na, N))
    D = sps.csr_matrix((Dval, (Drow, Dcol)), shape=(na, na))
    return Cw, Rw, D


def eigs_aug(sysm: S.System, Cw, Rw, D, sigma: float, k: int):
    """Lowest-k near sigma of (A + X) w = E B w with X = -Cw D^{-1} Rw (exact Schur complement)."""
    A, B = sysm.matrices()
    N = A.shape[0]
    if Cw is None:
        vals, vecs = spla.eigsh(A, k=k, M=B, sigma=sigma, which="LM")
    else:
        K = sps.bmat([[A - sigma * B, Cw], [Rw, D]], format="csc")
        lu = spla.splu(K)
        na = D.shape[0]

        def opinv(v):
            return lu.solve(np.concatenate([v, np.zeros(na)]))[:N]

        def amul(v):
            aux = spla.spsolve(D.tocsc(), -(Rw @ v))
            return A @ v + Cw @ aux
        OPinv = spla.LinearOperator((N, N), matvec=opinv, dtype=float)
        Aop = spla.LinearOperator((N, N), matvec=amul, dtype=float)
        vals, vecs = spla.eigsh(Aop, k=k, M=B, sigma=sigma, which="LM", OPinv=OPinv)
    o = np.argsort(vals)
    return vals[o], vecs[:, o]


def eft_vchan(g: S.Grid, fits: Dict[str, list], K: int, L: int) -> Dict[int, np.ndarray]:
    out = {}
    for l in range(L):
        c = fits[f"l{l}_K{K}"]
        if isinstance(c, str):
            raise ValueError(f"l={l} K={K} not constructible")
        out[l] = sum(ck * gg for ck, gg in zip(c, CA.gauss_terms(g.r, K)))
    return out


def setup(h: float = 0.01, t0: float = -16.0, rmax: float = 70.0, anchor: Optional[float] = None):
    g = S.Grid(t0=t0, t1=np.log(rmax), h=h, anchor=anchor)
    wc, Ec, Jc = CA.li_core(g)
    # reference atomic levels on this grid (channels decouple in the atom)
    levels = {}
    for l in (0, 1, 2):
        at = S.System(grid=g, lmax=l, Zd=3.0, vsph=2 * Jc)
        at.channels = lambda l=l: [l]                         # type: ignore[method-assign]
        Cw, Rw, D = exchange_aug(g, wc, [l], 8, CA.LAMBDA)
        vals, _ = eigs_aug(at, Cw, Rw, D, sigma=-0.25, k=5)
        levels[l] = [float(x) for x in vals if x < 0][:5]
    fits = {}
    for l in (0, 1, 2):
        c = None
        for K in range(3):
            c_init = None if c is None else np.concatenate([c, [0.0]])
            c = CA.fit_channel(g, l, np.array(levels[l]), K, c_init)
            fits[f"l{l}_K{K}"] = c.tolist()
    return g, wc, Jc, levels, fits


def molecule(g, wc, Jc, fits, R: float, lmax: int, sigma: float, k: int = 2):
    out = {}
    ref = S.System(grid=g, lmax=lmax, Zd=3.0, others=[(1.0, R)], vsph=2 * Jc)
    Cw, Rw, D = exchange_aug(g, wc, ref.channels(), 8, CA.LAMBDA)
    vals, _ = eigs_aug(ref, Cw, Rw, D, sigma=sigma, k=k)
    out["ref"] = vals.tolist()
    for K in range(3):
        for L in (1, 2, 3):
            eff = S.System(grid=g, lmax=lmax, Zd=1.0, others=[(1.0, R)], vchan=eft_vchan(g, fits, K, L))
            v, _ = eigs_aug(eff, None, None, None, sigma=sigma, k=k)
            out[f"K{K}_L{L}"] = v.tolist()
    # Paper 17 PK form on l = 0, A refit to the 2s reference level; local screening 2J, no exchange
    out["pk"] = pk_molecule(g, Jc, R, lmax, sigma, k, fits_levels=None)
    return out


def pk_A(g: S.Grid, Jc: np.ndarray, target2s: float) -> float:
    vpk = np.exp(-7.0 * g.r**2) / g.r**2
    A = 6.93
    for _ in range(60):
        vals, vecs, sysm = CA.channel_levels(g, 0, 2 * Jc + A * vpk, 1, Zd=3.0, sigma=target2s - 0.02)
        res = vals[0] - target2s
        if abs(res) < 1e-13:
            return A
        A -= res / sysm.inner(vecs[:, 0], vpk * vecs[:, 0])
    raise RuntimeError("PK refit failed")


PK_A: Dict[str, float] = {}


def pk_molecule(g, Jc, R, lmax, sigma, k, fits_levels=None):
    A = PK_A["A"]
    vpk = A * np.exp(-7.0 * g.r**2) / g.r**2
    sysm = S.System(grid=g, lmax=lmax, Zd=3.0, others=[(1.0, R)], vsph=2 * Jc, vchan={0: vpk})
    v, _ = eigs_aug(sysm, None, None, None, sigma=sigma, k=k)
    return v.tolist()


if __name__ == "__main__":
    lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    Rs = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [2.0, 3.0, 4.0, 6.0, 10.0]
    rec: Dict = {"lmax": lmax}
    for R in Rs:
        t = time.time()
        g, wc, Jc, levels, fits = setup(anchor=R)
        PK_A["A"] = pk_A(g, Jc, levels[0][0])
        sig = -1.5
        res = molecule(g, wc, Jc, fits, R, lmax, sigma=sig)
        res["pkA"] = PK_A["A"]
        rec[f"R{R}"] = res
        e = np.array(res["ref"])
        line = " ".join(f"{key}:{(np.array(v)[0]-e[0]):+.2e}" for key, v in res.items()
                        if key.startswith("K") or key == "pk")
        print(f"R={R} lmax={lmax} ref={e} ({time.time()-t:.0f}s)\n   X-state errors: {line}", flush=True)
        OUT.with_name(f"core_eft_B_l{lmax}.json").write_text(json.dumps(rec, indent=1))
