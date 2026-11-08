"""CE step A: Li+ frozen core, reference valence levels, effective-operator fits.

Record: docs/predictions/2026-10-20_core_eft.md (committed 97901bbd before computing).

Reference valence operator (exact frozen core): F = T - 3/r + 2J - K on the space orthogonal to the
core orbital (level shift Lambda |c><c|).  Effective operator at order (K, L):
  V_eff = -1/r + sum_{l<L} P_l sum_{k<=K} c_lk (r/a)^(2k) exp(-(r/a)^2),  a = 0.6.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sps

sys.path.insert(0, str(Path(__file__).parent))
import molfns_solver as S  # noqa: E402

OUT = Path(__file__).parent / "data" / "core_eft_A.json"
A_RANGE = 0.6
LAMBDA = 1e7


def li_core(g: S.Grid, iters: int = 60) -> Tuple[np.ndarray, float, np.ndarray]:
    """Restricted HF 1s^2 of Li+: F = h + J (2J - K with K = J). Returns w_c (normalized), E_total, J(r)."""
    sysm = S.System(grid=g, lmax=0, Zd=3.0)
    J = np.zeros(g.n)
    E_old = 0.0
    for _ in range(iters):
        sysm.vsph = J
        vals, vecs = S.ground_state(sysm, sigma=-3.0, k=1)
        w = vecs[:, 0]
        J = hartree(g, w)
        # E = 2<h> + <J>,  <h> = eps - <J>
        eps = float(vals[0])
        EJ = sysm.inner(w, J * w)
        E = 2 * (eps - EJ) + EJ
        if abs(E - E_old) < 1e-13:
            break
        E_old = E
    return w, E, J


def hartree(g: S.Grid, w: np.ndarray) -> np.ndarray:
    """J(r) of one electron with radial u = e^{t/2} w: (1/r) int_0^r u^2 + int_r^inf u^2 / r'."""
    r = g.r
    u2dr = np.exp(2 * g.t) * w**2 * g.h               # u^2 dr = e^{t} w^2 e^{t} dt
    inner = np.cumsum(u2dr) - 0.5 * u2dr
    outer = np.cumsum((u2dr / r)[::-1])[::-1] - 0.5 * u2dr / r
    return inner / r + outer


def exchange_blocks(g: S.Grid, wc: np.ndarray, ls: List[int], lmaxK: int, nc: int,
                    cut: float = 1e-13) -> sps.csr_matrix:
    """Per-channel exchange, t-major indexing: M_ij = (1/(2l+1)) e^{2t_i} wc_i K_l(r_i, r_j) e^{2t_j} wc_j h."""
    r, e2t = g.r, np.exp(2 * g.t)
    f = e2t * wc
    keep = np.where(np.abs(f) > cut * np.abs(f).max())[0]
    rk = r[keep]
    rows, cols, vals = [], [], []
    for c, l in enumerate(ls):
        if l > lmaxK:
            continue
        rl, rg = np.minimum.outer(rk, rk), np.maximum.outer(rk, rk)
        Kl = rl**l / rg ** (l + 1)
        M = (1.0 / (2 * l + 1)) * np.outer(f[keep], f[keep]) * Kl * g.h
        ii = keep * nc + c
        I, Jx = np.meshgrid(ii, ii, indexing="ij")
        rows.append(I.ravel()); cols.append(Jx.ravel()); vals.append(M.ravel())
    N = g.n * nc
    if not rows:
        return sps.csr_matrix((N, N))
    return sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))


def level_shift(g: S.Grid, wc: np.ndarray, ls: List[int], nc: int, lam: float = LAMBDA) -> sps.csr_matrix:
    """Lambda B|c><c|B in the s channel."""
    e2t = np.exp(2 * g.t)
    f = e2t * wc * np.sqrt(g.h)
    keep = np.where(np.abs(f) > 1e-15 * np.abs(f).max())[0]
    c = ls.index(0)
    ii = keep * nc + c
    I, Jx = np.meshgrid(ii, ii, indexing="ij")
    M = lam * np.outer(f[keep], f[keep])
    N = g.n * nc
    return sps.csr_matrix((M.ravel(), (I.ravel(), Jx.ravel())), shape=(N, N))


def reference_system(g: S.Grid, wc: np.ndarray, Jc: np.ndarray, lmax: int, others=(), lmaxK: int = 8):
    sysm = S.System(grid=g, lmax=lmax, Zd=3.0, others=list(others), vsph=2 * Jc)
    ls = sysm.channels()
    nc = len(ls)
    extra = -exchange_blocks(g, wc, ls, lmaxK, nc) + level_shift(g, wc, ls, nc)
    return sysm, extra


def gauss_terms(r: np.ndarray, K: int) -> List[np.ndarray]:
    x = (r / A_RANGE) ** 2
    return [x**k * np.exp(-x) for k in range(K + 1)]


def channel_levels(g: S.Grid, l: int, vextra: np.ndarray, nlev: int, Zd: float = 1.0, sigma: float = -0.6):
    """Levels of -Zd/r + vextra in one channel l (single-channel system with only that l)."""
    sysm = S.System(grid=g, lmax=l, Zd=Zd, vsph=None, parity=None)
    sysm.lmax = l
    sysm.channels = lambda: [l]                              # type: ignore[method-assign]
    sysm.vchan = {l: vextra}
    vals, vecs = S.ground_state(sysm, sigma=sigma, k=nlev)
    return vals, vecs, sysm


def _levels_jac(g: S.Grid, l: int, c: np.ndarray, K: int, sig: float, nlev: int):
    gk = gauss_terms(g.r, K)
    v = sum(ck * gg for ck, gg in zip(c, gk))
    vals, vecs, sysm = channel_levels(g, l, v, nlev, sigma=sig)
    Jm = np.array([[sysm.inner(vecs[:, n], gg * vecs[:, n]) for gg in gk] for n in range(K + 1)])
    return vals, Jm


def fit_channel(g: S.Grid, l: int, targets: np.ndarray, K: int, c0: Optional[np.ndarray] = None,
                tol: float = 1e-13, steps: int = 40) -> np.ndarray:
    """c (K+1 values) so the K+1 lowest levels of -1/r + sum c_k g_k equal targets[:K+1].
    Branch rule (disclosed in the disposition): K = 0 by Newton from c = 1 (s) or 0 (l > 0); K >= 1 by
    continuation from the K-1 solution with c_K = 0, moving the targets linearly from the levels
    there to the reference levels in `steps` steps, Newton at each step."""
    sig = min(targets[0] - 0.05, -0.6) if l == 0 else targets[0] - 0.05
    if c0 is None:
        c = np.array([1.0 if l == 0 else 0.0] + [0.0] * K)
        L0 = None
    else:
        c = c0.copy()
        L0 = _levels_jac(g, l, c, K, sig, K + 1)[0][: K + 1]
    taus = [1.0] if L0 is None else np.linspace(0, 1, steps + 1)[1:]
    for tau in taus:
        tgt = targets[: K + 1] if L0 is None else (1 - tau) * L0 + tau * targets[: K + 1]
        for it in range(60):
            vals, Jm = _levels_jac(g, l, c, K, sig, K + 1)
            res = vals[: K + 1] - tgt
            if np.max(np.abs(res)) < tol:
                break
            c = c + np.linalg.solve(Jm, -res)
        else:
            raise RuntimeError(f"fit l={l} K={K} stuck at tau={tau}: res={res}")
    return c


if __name__ == "__main__":
    rec: Dict = {}
    g = S.Grid(t0=-22.0, t1=np.log(220.0), h=0.005)
    wc, Ecore, Jc = li_core(g)
    print(f"C1: Li+ HF E = {Ecore:.10f}  (ref -7.2364152, diff {Ecore + 7.2364152:.2e})", flush=True)
    rec["C1"] = {"E": Ecore, "diff": Ecore + 7.2364152}
    # reference atomic levels, channel by channel (the atom is spherical: channels decouple)
    levels: Dict[int, List[float]] = {}
    for l in (0, 1, 2, 3):
        sysm, extra = reference_system(g, wc, Jc, lmax=l)
        sysm.channels = lambda l=l: [l]                     # type: ignore[method-assign]
        sysm2, extra2 = sysm, -exchange_blocks(g, wc, [l], 8, 1) + (level_shift(g, wc, [l], 1) if l == 0 else 0)
        vals, _ = S.ground_state(sysm2, sigma=-0.25, k=5, extra=extra2)
        levels[l] = [float(x) for x in vals if x < 0][:5]
        print(f"ref l={l}: {levels[l]}", flush=True)
    rec["levels"] = levels
    # C4: pure -1/r (hydrogenic) 2s level against the reference 2s
    rec["C4"] = {"hydrogenic_2s": -0.125, "ref_2s": levels[0][0], "diff": -0.125 - levels[0][0]}
    print(f"C4: hydrogenic n=2 minus ref 2s = {-0.125 - levels[0][0]:.4f}", flush=True)
    # EFT fits, K = 0..3, l = 0..2
    fits: Dict[str, List[float]] = {}
    for l in (0, 1, 2):
        c = None
        for K in range(4):
            c_init = None if c is None else np.concatenate([c, [0.0]])
            try:
                c = fit_channel(g, l, np.array(levels[l]), K, c_init)
            except (RuntimeError, np.linalg.LinAlgError) as e:
                fits[f"l{l}_K{K}"] = f"NOT CONSTRUCTIBLE: {e}"
                print(f"fit l={l} K={K}: NOT CONSTRUCTIBLE ({e})", flush=True)
                break
            fits[f"l{l}_K{K}"] = c.tolist()
            v = sum(ck * gg for ck, gg in zip(c, gauss_terms(g.r, K)))
            vals, _, _ = channel_levels(g, l, v, K + 2, sigma=levels[l][0] - 0.05)
            print(f"fit l={l} K={K}: c={np.round(c, 6).tolist()}  matched-res="
                  f"{np.max(np.abs(vals[:K+1]-levels[l][:K+1])):.1e}  next level eff={vals[K+1]:.8f} "
                  f"ref={levels[l][K+1]:.8f}", flush=True)
    rec["fits"] = fits
    OUT.write_text(json.dumps(rec, indent=1))
