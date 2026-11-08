"""Probe B -- does the ANALYTIC momentum factorization beat the standard block-encoding
1-norm?  (the "analytic THC" question)

Pre-registered killer test.  The standard measure (v4.94.0,
``geovac/sturmian_molecular_lambda.py``) is

    lambda_std = sum_pq |h_pq| + sum_pqrs |(pq|rs)|        [Loewdin-orthonormal basis]

The factorized side uses the momentum representation (Paper 59 ``sec:f12``):

    (pq|rs) = (1/pi) int_0^inf dk int_{-1}^{1} dmu  rho~_pq^*(k,mu) rho~_rs(k,mu)
            = sum_nu omega_nu [ C^nu_pq C^nu_rs + S^nu_pq S^nu_rs ]

with C = Re rho~, S = Im rho~ (real SYMMETRIC n_orb x n_orb matrices at each node)
and omega_nu = (2/pi) w_k w_mu > 0 on the half-space mu in [0,1] (the mu<0 half is
the complex conjugate).  Each node contributes a SQUARED Hermitian one-body
operator -- exactly the single-factorization (SF) / double-factorization (DF)
block-encoding structure, with a DIAGONAL core.

THREE 1-norm conventions, all directly comparable to sum_pqrs |(pq|rs)|:

  lambda_abs   = sum_nu omega_nu ( sum_pq |rho~_pq| )^2         [protocol convention]
  lambda_absCS = sum_nu omega_nu [ (sum_pq|C_pq|)^2 + (sum_pq|S_pq|)^2 ]
  lambda_spec  = sum_nu omega_nu [ ||C||_*^2 + ||S||_*^2 ]      [Berry/vBurg DF std]

||.||_* = nuclear norm (sum |eigenvalues|) = the honest cost after the Givens
rotation that diagonalizes each one-body leaf.

ORDERING THEOREM (proof in debug/probeB_findings.md; checked numerically here):
    lambda_spec <= lambda_absCS <= lambda_abs   and   lambda_absCS >= lambda_std
so the ONLY convention under which the factorization can win is lambda_spec.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
from numpy.polynomial.legendre import leggauss

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from debug.probeB_thc_density import (SOrb, _dF, feynman_grid,  # noqa: E402
                                      lowdin_X, rho_tilde, system)


# ================================================================= generic rho~
def rho_tilde_terms(termsA, aA: float, zA: float, NA: float,
                    termsB, aB: float, zB: float, NB: float,
                    k, mu, tgrid, chunk: int = 32) -> np.ndarray:
    """rho~ for GENERIC radial term lists (powers m >= -1 allowed).

    int r_A^{m1} r_B^{m2} e^{-a r_A - b r_B} e^{-i k.r} d^3r
        = (-1)^{m1+m2} d_a^{m1+1} d_b^{m2+1} G,   G = FT of the Yukawa product.

    m = -1 is the un-differentiated Yukawa, which is what the r^{-1} term of the
    kinetic operator -1/2 grad^2 needs.
    """
    t, wt = tgrid
    k = np.atleast_1d(np.asarray(k, float))
    mu = np.atleast_1d(np.asarray(mu, float))
    D = abs(zA - zB)
    Pz = (1.0 - t) * zA + t * zB
    out = np.zeros((k.size, mu.size), dtype=complex)
    for i0 in range(0, k.size, chunk):
        ks = np.maximum(k[i0:i0 + chunk], 1e-12)
        G = np.zeros((t.size, ks.size))
        for (m1, c1) in termsA:
            for (m2, c2) in termsB:
                val = _dF(m1 + 1, m2 + 1)(aA, aB, t[:, None], ks[None, :], D)
                G = G + ((-1) ** (m1 + m2)) * c1 * c2 * np.asarray(val, float)
        G = G * wt[:, None]
        if D < 1e-14:
            rad = G.sum(axis=0)
            out[i0:i0 + chunk, :] = rad[:, None] * np.exp(-1j * np.outer(ks, mu) * zA)
        else:
            ph = np.exp(-1j * (ks[None, :, None] * mu[None, None, :])
                        * Pz[:, None, None])
            out[i0:i0 + chunk, :] = np.einsum("tk,tkm->km", G, ph)
    return 2.0 * math.pi * NA * NB * out


def laplacian_terms(o: SOrb):
    """Radial term list of  -1/2 grad^2 phi_o  (same N and rate a as phi_o).

    grad^2 [r^m e^{-a r}] = [ m(m+1) r^{m-2} - 2a(m+1) r^{m-1} + a^2 r^m ] e^{-a r}
    """
    acc: Dict[int, float] = {}
    for m, c in o.terms:
        for pw, co in ((m - 2, -0.5 * m * (m + 1)),
                       (m - 1, o.a * (m + 1)),
                       (m, -0.5 * o.a ** 2)):
            if co != 0.0:
                acc[pw] = acc.get(pw, 0.0) + c * co
    return tuple((pw, co) for pw, co in sorted(acc.items()) if abs(co) > 1e-300)


# ============================================================ one-electron block
def _gl_panels(lo: float, hi: float, n_pan: int, n_g: int, power: float = 1.0):
    x, wx = leggauss(n_g)
    edges = lo + (hi - lo) * np.linspace(0.0, 1.0, n_pan + 1) ** power
    pts, wts = [], []
    for a_, b_ in zip(edges[:-1], edges[1:]):
        m, h = 0.5 * (a_ + b_), 0.5 * (b_ - a_)
        pts.append(m + h * x)
        wts.append(h * wx)
    return np.concatenate(pts), np.concatenate(wts)


def nuclear_pair(p: SOrb, q: SOrb, zC: float, Z: float) -> float:
    """<p| -Z/|r-C| |q>, C on the z axis.  Exact 1-D / 2-D real-space quadrature."""
    if abs(p.z - q.z) < 1e-14:                                  # one-centre density
        d = abs(zC - p.z)
        rate = p.a + q.a
        if d < 1e-14:
            r, w = _gl_panels(0.0, 80.0 / rate, 80, 24)
            return -Z * 4 * math.pi * float(np.sum(w * p.radial(r) * q.radial(r) * r))
        r1, w1 = _gl_panels(0.0, d, 40, 24)
        r2, w2 = _gl_panels(d, d + 80.0 / rate, 80, 24)
        inner = float(np.sum(w1 * p.radial(r1) * q.radial(r1) * r1 ** 2)) / d
        outer = float(np.sum(w2 * p.radial(r2) * q.radial(r2) * r2))
        return -Z * 4 * math.pi * (inner + outer)
    c = 0.5 * abs(p.z - q.z)
    hi_z, lo_z = (p, q) if p.z > q.z else (q, p)
    rate = hi_z.a + lo_z.a
    xi, wxi = _gl_panels(1.0, 1.0 + 80.0 / (c * rate), 80, 24, power=1.5)
    eta, weta = _gl_panels(-1.0, 1.0, 8, 40)
    XI, ETA = xi[:, None], eta[None, :]
    f = hi_z.radial(c * (XI - ETA)) * lo_z.radial(c * (XI + ETA))
    W = wxi[:, None] * weta[None, :]
    if abs(zC - hi_z.z) < 1e-12:                                # nucleus at focus F1
        val = 2 * math.pi * c ** 2 * float(np.sum(W * (XI + ETA) * f))
    elif abs(zC - lo_z.z) < 1e-12:                              # nucleus at focus F2
        val = 2 * math.pi * c ** 2 * float(np.sum(W * (XI - ETA) * f))
    else:
        raise NotImplementedError("nucleus off both foci")
    return -Z * val


def one_electron(orbs: Sequence[SOrb], nuclei, tgrid):
    """(S, h) in the RAW hydrogenic basis; h = T + V_ne."""
    n = len(orbs)
    S = np.zeros((n, n))
    T = np.zeros((n, n))
    V = np.zeros((n, n))
    k0, mu0 = np.array([1e-12]), np.array([0.0])
    for i in range(n):
        for j in range(n):
            if j >= i:
                S[i, j] = S[j, i] = rho_tilde(orbs[i], orbs[j], k0, mu0,
                                              tgrid)[0, 0].real
                V[i, j] = V[j, i] = sum(nuclear_pair(orbs[i], orbs[j], zc, zz)
                                        for zc, zz in nuclei)
            T[i, j] = rho_tilde_terms(orbs[i].terms, orbs[i].a, orbs[i].z,
                                      orbs[i].norm, laplacian_terms(orbs[j]),
                                      orbs[j].a, orbs[j].z, orbs[j].norm,
                                      k0, mu0, tgrid)[0, 0].real
    T = 0.5 * (T + T.T)
    return S, T + V


# ==================================================================== k/mu grid
def build_grid(Kmax: float, n_pan: int, n_g: int, dP_max: float,
               pad: int = 12, n_mu_min: int = 6, power: float = 1.6):
    """Nodes (k_i, mu_j) with a k-ADAPTIVE mu order.

    The mu-integrand of conj(rho~_pq) rho~_rs carries the net phase
    e^{-i k mu (P_rs - P_pq)}, |dP| <= dP_max = max centre separation, so the
    Gauss-Legendre order must satisfy n_mu >~ k dP/pi + pad (calibrated: GL on
    [0,1] integrates cos(a mu) to 1e-13 at n = a/pi + 12).  A product grid wastes
    nodes at small k; a QROM over a flat node list needs no product structure, so
    the adaptive count is the honest M.

    Uses rho~(k,-mu) = conj(rho~(k,mu)): integrate mu in [0,1], weight doubled.
    """
    k, wk = _gl_panels(0.0, Kmax, n_pan, n_g, power=power)
    orders = np.array([6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80,
                       96, 112, 128, 160, 192, 224, 256, 320, 384, 448, 512, 640,
                       768, 896, 1024])
    need = np.maximum(n_mu_min, np.ceil(k * dP_max / math.pi) + pad)
    idx = np.clip(np.searchsorted(orders, need), 0, len(orders) - 1)
    bands = []
    for oi in np.unique(idx):
        sel = idx == oi
        m = int(orders[oi])
        x, wx = leggauss(m)
        bands.append((k[sel], wk[sel], 0.5 * (x + 1.0), 0.5 * wx))
    return bands


def node_count(bands) -> int:
    return int(sum(len(kb) * len(mb) for kb, _, mb, _ in bands))


def kof_bands(bands) -> np.ndarray:
    return np.concatenate([np.repeat(kb, len(mb)) for kb, _, mb, _ in bands])


# ============================================================ densities on grid
def grid_densities(orbs: Sequence[SOrb], bands, tgrid, X: np.ndarray):
    """(omega (M,), C (M,n,n), S (M,n,n)) in the Loewdin basis."""
    n = len(orbs)
    om_all, C_all, S_all = [], [], []
    for (kb, wkb, mub, wmub) in bands:
        raw = np.zeros((n, n, kb.size, mub.size), dtype=complex)
        for i in range(n):
            for j in range(i, n):
                v = rho_tilde(orbs[i], orbs[j], kb, mub, tgrid)
                raw[i, j] = v
                raw[j, i] = v
        rot = np.einsum("pi,ijkm,jq->pqkm", X, raw, X, optimize=True)
        rot = np.transpose(rot, (2, 3, 0, 1)).reshape(-1, n, n)
        om = ((2.0 / math.pi) * wkb[:, None] * wmub[None, :]).reshape(-1)
        om_all.append(om)
        C_all.append(rot.real.copy())
        S_all.append(rot.imag.copy())
    return (np.concatenate(om_all), np.concatenate(C_all, axis=0),
            np.concatenate(S_all, axis=0))


def eri_from_nodes(om, C, S) -> np.ndarray:
    return (np.einsum("n,npq,nrs->pqrs", om, C, C, optimize=True)
            + np.einsum("n,npq,nrs->pqrs", om, S, S, optimize=True))


def lambdas(om, C, S, n_elec: int = 2) -> Dict[str, float]:
    """The three head-to-head conventions plus the identity-shifted variant.

    IDENTITY SHIFT (the small-k cure).  As k -> 0, rho~ -> S = I in an
    orthonormal basis, so the leaf becomes the number operator N and its square
    is a c-number on a fixed particle-number sector -- the small-k 1-norm mass is
    spurious.  Splitting  Chat = c I + Ctilde,  c = tr C / n,  tr Ctilde = 0:

        Chat^2 = c^2 N^2 (constant)  +  2 c N Ctilde (one-body)  +  Ctilde^2

    so the operator-carrying 1-norm is 2 n_elec |c| ||Ctilde||_* + ||Ctilde||_*^2.
    """
    n = C.shape[1]
    absC = np.abs(C).sum(axis=(1, 2))
    absS = np.abs(S).sum(axis=(1, 2))
    absR = np.abs(C + 1j * S).sum(axis=(1, 2))
    nucC = np.abs(np.linalg.eigvalsh(C)).sum(axis=1)
    nucS = np.abs(np.linalg.eigvalsh(S)).sum(axis=1)
    cC = np.trace(C, axis1=1, axis2=2) / n
    cS = np.trace(S, axis1=1, axis2=2) / n
    eye = np.eye(n)[None, :, :]
    ntC = np.abs(np.linalg.eigvalsh(C - cC[:, None, None] * eye)).sum(axis=1)
    ntS = np.abs(np.linalg.eigvalsh(S - cS[:, None, None] * eye)).sum(axis=1)
    shift = (2 * n_elec * (np.abs(cC) * ntC + np.abs(cS) * ntS)
             + ntC ** 2 + ntS ** 2)
    return {"abs": float(np.sum(om * absR ** 2)),
            "absCS": float(np.sum(om * (absC ** 2 + absS ** 2))),
            "spec": float(np.sum(om * (nucC ** 2 + nucS ** 2))),
            "shift": float(np.sum(om * shift))}


# ======================================================================== main
#  label : (Kmax, n_k_panels, n_g, mu pad)
GRIDS = {
    "g0": (8.0, 10, 6, 1),
    "g1": (14.0, 16, 6, 2),
    "g2": (22.0, 24, 8, 3),
    "g3": (34.0, 34, 8, 5),
    "g4": (55.0, 55, 10, 8),
    "g5": (90.0, 90, 10, 10),
    "ref": (220.0, 190, 12, 14),
}
ORDER = ("g0", "g1", "g2", "g3", "g4", "g5", "ref")


NELEC = {"H2": 2, "LiH": 4, "H2_4o": 2, "LiH_5o": 4}


def run(sysname: str, tgrid, out: List[str]):
    orbs, nuclei = system(sysname)
    n_elec = NELEC[sysname]
    n = len(orbs)
    dP = max(max(abs(o1.z - o2.z) for o1 in orbs for o2 in orbs), 1e-6)

    S, h_raw = one_electron(orbs, nuclei, tgrid)
    X = lowdin_X(S)
    hm = X.T @ h_raw @ X
    lam1 = float(np.abs(hm).sum())
    ev = np.linalg.eigvalsh(S)

    out.append("")
    out.append(f"### {sysname}   n_orb = {n}   centres {[o.z for o in orbs]}   "
               f"dP_max = {dP:.4f}")
    out.append(f"    S eigenvalues {np.array2string(ev, precision=6)}   "
               f"cond(S) = {ev[-1]/ev[0]:.4g}")
    out.append(f"    lambda_1body = sum_pq|h_pq| (Loewdin) = {lam1:.6f}")

    results = {}
    for label in ORDER:
        Kmax, n_pan, n_g, pad = GRIDS[label]
        t0 = time.time()
        bands = build_grid(Kmax, n_pan, n_g, dP, pad=pad)
        om, C, Si = grid_densities(orbs, bands, tgrid, X)
        eri = eri_from_nodes(om, C, Si)
        lam = lambdas(om, C, Si, n_elec)
        lam["std"] = float(np.abs(eri).sum())
        results[label] = dict(M=node_count(bands), Kmax=Kmax, pad=pad, eri=eri,
                              lam=lam, om=om, C=C, S=Si, kof=kof_bands(bands),
                              secs=time.time() - t0)
    ref = results["ref"]

    out.append("")
    out.append("    grid convergence of the reconstructed ERI tensor and of the 1-norms")
    out.append("    grid  Kmax  pad   M nodes  log2 M   max|d(pq|rs)|    lam2_std    "
               "lam2_spec   lam2_shift   lam2_absCS    lam2_abs   secs")
    for label in ORDER:
        r = results[label]
        err = float(np.abs(r["eri"] - ref["eri"]).max())
        out.append(f"    {label:>4} {r['Kmax']:>5.0f} {r['pad']:>4}  {r['M']:>8d}  "
                   f"{math.log2(max(r['M'], 1)):>6.2f}   {err:>13.3e}  "
                   f"{r['lam']['std']:>11.5f} {r['lam']['spec']:>11.5f} "
                   f"{r['lam']['shift']:>11.5f} "
                   f"{r['lam']['absCS']:>12.5f} {r['lam']['abs']:>11.5f} "
                   f"{r['secs']:>6.1f}")

    out.append("")
    out.append("    M(eps): smallest tested grid whose max ERI deviation <= eps")
    for eps in (1e-4, 1e-6, 1e-8):
        pick = None
        for label in ORDER:
            if float(np.abs(results[label]["eri"] - ref["eri"]).max()) <= eps:
                pick = label
                break
        if pick is None:
            out.append(f"      eps={eps:.0e}: no tested grid reaches it")
        else:
            r = results[pick]
            out.append(f"      eps={eps:.0e}: grid {pick}  M = {r['M']}  "
                       f"(log2 M = {math.log2(r['M']):.2f} -> "
                       f"{math.ceil(math.log2(r['M']))} k/mu ancillas, +1 C/S bit)"
                       f"   lam2_spec = {r['lam']['spec']:.5f}"
                       f"   lam2_std = {r['lam']['std']:.5f}")

    L = ref["lam"]
    ok = (L["spec"] <= L["absCS"] <= L["abs"]) and (L["absCS"] >= L["std"])
    out.append("")
    out.append(f"    ordering check: spec {L['spec']:.5f} <= absCS {L['absCS']:.5f} "
               f"<= abs {L['abs']:.5f};  absCS >= std {L['std']:.5f}  -> "
               f"{'OK' if ok else 'VIOLATED'}")

    om, C, Si, kof = ref["om"], ref["C"], ref["S"], ref["kof"]
    absC = np.abs(C).sum(axis=(1, 2))
    absS = np.abs(Si).sum(axis=(1, 2))
    nucC = np.abs(np.linalg.eigvalsh(C)).sum(axis=1)
    nucS = np.abs(np.linalg.eigvalsh(Si)).sum(axis=1)
    per_spec = om * (nucC ** 2 + nucS ** 2)
    per_absCS = om * (absC ** 2 + absS ** 2)
    out.append("")
    out.append("    small-k audit -- where the factorized 1-norm mass sits, and how much")
    out.append("    ERI value the k < k_c band carries (i.e. what a classical handling")
    out.append("    of that band would have to reproduce)")
    out.append("    k_c    frac lam2_spec   frac lam2_absCS   lam2_spec(k>k_c)   "
               "max|(pq|rs)| from k<k_c")
    for kc in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        m = kof <= kc
        drop = eri_from_nodes(om[m], C[m], Si[m])
        out.append(f"    {kc:>5.2f}   {per_spec[m].sum()/per_spec.sum():>13.4f}   "
                   f"{per_absCS[m].sum()/per_absCS.sum():>15.4f}   "
                   f"{per_spec[~m].sum():>16.5f}   {np.abs(drop).max():>18.4e}")
    return dict(n=n, lam1=lam1, ref=ref, results=results, cond=ev[-1] / ev[0])


if __name__ == "__main__":
    which = sys.argv[1:] or ["H2", "LiH", "H2_4o"]
    tgrid = feynman_grid(y_max=26.0, panel=0.30, n_g=12)
    out: List[str] = []
    out.append("Probe B -- analytic-momentum (SF/THC-style) vs standard block-encoding "
               "1-norm")
    out.append(f"Feynman t-grid: {tgrid[0].size} nodes (two-centre pairs only; "
               "one-centre densities use the closed-form rational FT)")
    summary = {}
    shown = 0
    for sname in which:
        summary[sname] = run(sname, tgrid, out)
        print("\n".join(out[shown:]))
        shown = len(out)
        sys.stdout.flush()

    out.append("")
    out.append("=" * 104)
    out.append("HEAD-TO-HEAD (converged 'ref' grid, Loewdin-orthonormal basis; "
               "lambda = lambda_1body + lambda_2body)")
    out.append("system   n_orb   lam_std      lam_THC(spec)  lam_THC(absCS)  "
               "lam_THC(abs)  lam_THC(shift)   spec/std  shift/std")
    ns, ls, lt = [], [], []
    for sname, d in summary.items():
        l1 = d["lam1"]
        L = d["ref"]["lam"]
        tot = l1 + L["std"]
        out.append(f"{sname:>7}  {d['n']:>5}   {tot:>10.5f}   {l1+L['spec']:>12.5f}  "
                   f"{l1+L['absCS']:>13.5f}  {l1+L['abs']:>12.5f}   "
                   f"{l1+L['shift']:>13.5f}   "
                   f"{(l1+L['spec'])/tot:>7.3f}  {(l1+L['shift'])/tot:>9.3f}")
        ns.append(d["n"])
        ls.append(tot)
        lt.append(l1 + L["spec"])
    if len(ns) >= 3:
        p_std = float(np.polyfit(np.log(ns), np.log(ls), 1)[0])
        p_thc = float(np.polyfit(np.log(ns), np.log(lt), 1)[0])
        out.append("")
        out.append(f"n_orb trend (3 points at n_orb = {ns}, INDICATIVE ONLY): "
                   f"lam_std ~ n^{p_std:.2f},  lam_THC(spec) ~ n^{p_thc:.2f}")
    txt = "\n".join(out)
    Path("debug/data").mkdir(parents=True, exist_ok=True)
    Path("debug/data/probeB_thc_run.txt").write_text(txt, encoding="utf-8")
    print("\n".join(out[shown:]))
    print("\n[written] debug/data/probeB_thc_run.txt")
