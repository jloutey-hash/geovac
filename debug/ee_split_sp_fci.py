"""The l>0 payoff curve: does the e-e split's rank compression survive angular coupling?

Builds a Gaunt-coupled s+p one-centre FCI on shared-k Coulomb-Sturmians, splits EVERY
multipole channel by the min/max identity, truncates each channel's W_L by rank, and
measures the FCI error.  s-only measured 1 mHa at W-rank 3-4 (basis-independent); the
question is whether that survives when L = 0, 1, 2 channels are all active.

  <ab|1/r12|cd> = sum_{LM} (4pi/(2L+1)) gA(a,LM,c) gB(b,LM,d) R^L(ac|bd)
  R^L kernel  r_<^L / r_>^{L+1}  =  K_L,sep - K_L,W     (exact, per channel)

Angular factors from the TRACKED geovac.xtc_angular_sparsity (exact wigner3j).
Complex harmonics are rotated to REAL harmonics so the tensor is real -- the size of
the discarded imaginary part is itself a gate.

VALIDATION LADDER (all must pass before any payoff number is quoted)
  G1  s-only sector reproduces the independent s-only engine (transcorrelated_sturmian)
  G2  assembled tensor is real (imag ~ 0) and obeys the 8-fold real-orbital symmetry
  G3  He FCI: s+p is BELOW s-only (more basis) and ABOVE the exact -2.9037243770
  G4  gamma-free sanity: the split is exact at tensor level, per L and in total
"""
import io
import json
import os
import sys
from math import factorial, pi, sqrt

import numpy as np
from scipy.linalg import eigh
from scipy.special import eval_genlaguerre

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from geovac import transcorrelated_sturmian as TC          # noqa: E402
from geovac.xtc_angular_sparsity import gA, gB             # noqa: E402

FOURPI = 4.0 * pi


def R_nl(n, l, r, k):
    x = 2 * k * r
    norm = sqrt((2 * k) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * x ** l * np.exp(-x / 2) * eval_genlaguerre(n - l - 1, 2 * l + 1, x)


def dR_nl(n, l, r, k):
    """ANALYTIC dR/dr (dL^a_m/dx = -L^{a+1}_{m-1}); finite differences cost 4e-4 in h1."""
    x = 2 * k * r
    norm = sqrt((2 * k) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    m = n - l - 1
    L = eval_genlaguerre(m, 2 * l + 1, x)
    dL = -eval_genlaguerre(m - 1, 2 * l + 2, x) if m >= 1 else np.zeros_like(x)
    e = np.exp(-x / 2)
    xl = x ** l
    xlm1 = x ** (l - 1) if l >= 1 else np.zeros_like(x)
    dRdx = norm * e * (l * xlm1 * L - 0.5 * xl * L + xl * dL)
    return 2 * k * dRdx


def real_harmonic_transform(orbs):
    """U[real, complex] mapping complex-Y orbitals to real-Y orbitals (block per (n,l))."""
    n_orb = len(orbs)
    U = np.zeros((n_orb, n_orb), dtype=complex)
    for a, (na, la, ma) in enumerate(orbs):
        for b, (nb, lb, mb) in enumerate(orbs):
            if (na, la) != (nb, lb):
                continue
            if ma == 0:
                U[a, b] = 1.0 if mb == 0 else 0.0
            elif ma > 0:                       # cosine-type
                if mb == ma:
                    U[a, b] = ((-1) ** ma) / sqrt(2)
                elif mb == -ma:
                    U[a, b] = 1.0 / sqrt(2)
            else:                              # sine-type, ma < 0
                mm = -ma
                if mb == mm:
                    U[a, b] = -1j * ((-1) ** mm) / sqrt(2)
                elif mb == -mm:
                    U[a, b] = 1j / sqrt(2)
    return U


def build(ns, npp, k, Z, Ng=700, Lmax=2):
    """Return orbitals, S, h1, and dict of (g, g_sep, per-L radial W blocks)."""
    r, wr = TC.make_grid(k, Ng=Ng)
    W2 = r * r * wr
    orbs = [(n, 0, 0) for n in range(1, ns + 1)]
    orbs += [(n, 1, m) for n in range(2, 2 + npp) for m in (-1, 0, 1)]
    n_orb = len(orbs)
    rad = {(n, l): R_nl(n, l, r, k) for (n, l, m) in orbs}
    drad = {(n, l): dR_nl(n, l, r, k) for (n, l, m) in orbs}

    # one-body: S and h1 = T + V  (gradient-form T with the centrifugal term)
    S = np.zeros((n_orb, n_orb))
    h1 = np.zeros((n_orb, n_orb))
    for a, (na, la, ma) in enumerate(orbs):
        for b, (nb, lb, mb) in enumerate(orbs):
            if (la, ma) != (lb, mb):
                continue
            Ra, Rb = rad[(na, la)], rad[(nb, lb)]
            dRa, dRb = drad[(na, la)], drad[(nb, lb)]
            S[a, b] = np.sum(Ra * Rb * W2)
            T = 0.5 * np.sum(dRa * dRb * W2) \
                + 0.5 * la * (la + 1) * np.sum(Ra * Rb * wr)
            V = -Z * np.sum(Ra * Rb * r * wr)
            h1[a, b] = T + V

    # radial pair densities and per-L Slater integrals (full + split)
    R1g, R2g = np.meshgrid(r, r, indexing="ij")
    lo, hi = np.minimum(R1g, R2g), np.maximum(R1g, R2g)
    radpairs = sorted({(min((na, la), (nb, lb)), max((na, la), (nb, lb)))
                       for (na, la, _) in orbs for (nb, lb, _) in orbs})
    rp_index = {p: i for i, p in enumerate(radpairs)}
    P = np.array([rad[p0] * rad[p1] * W2 for (p0, p1) in radpairs])

    RL, RL_sep, RL_W = {}, {}, {}
    for L in range(Lmax + 1):
        KL = lo ** L / hi ** (L + 1)
        t1 = R1g ** L / R2g ** (L + 1)
        t2 = R2g ** L / R1g ** (L + 1)
        Ksep, Kw = 0.5 * (t1 + t2), 0.5 * np.abs(t1 - t2)
        RL[L] = P @ KL @ P.T
        RL_sep[L] = P @ Ksep @ P.T
        RL_W[L] = P @ Kw @ P.T
    return orbs, rad, S, h1, rp_index, RL, RL_sep, RL_W


def assemble(orbs, rp_index, RLd, Lmax=2):
    """<ab|cd> from per-L radial blocks; returns the REAL-harmonic tensor + imag size."""
    n_orb = len(orbs)
    g = np.zeros((n_orb,) * 4, dtype=complex)
    ang = {}
    for L in range(Lmax + 1):
        pref = FOURPI / (2 * L + 1)
        for M in range(-L, L + 1):
            A = np.array([[gA(la, ma, L, M, lc, mc) for (nc, lc, mc) in orbs]
                          for (na, la, ma) in orbs])
            B = np.array([[gB(lb, mb, L, M, ld, md) for (nd, ld, md) in orbs]
                          for (nb, lb, mb) in orbs])
            ang[(L, M)] = (pref, A, B)
    for L in range(Lmax + 1):
        RLmat = RLd[L]
        for M in range(-L, L + 1):
            pref, A, B = ang[(L, M)]
            for a, (na, la, ma) in enumerate(orbs):
                for c, (nc, lc, mc) in enumerate(orbs):
                    if abs(A[a, c]) < 1e-14:
                        continue
                    iac = rp_index[(min((na, la), (nc, lc)), max((na, la), (nc, lc)))]
                    for b, (nb, lb, mb) in enumerate(orbs):
                        for d, (nd, ld, md) in enumerate(orbs):
                            if abs(B[b, d]) < 1e-14:
                                continue
                            ibd = rp_index[(min((nb, lb), (nd, ld)),
                                            max((nb, lb), (nd, ld)))]
                            g[a, b, c, d] += pref * A[a, c] * B[b, d] * RLmat[iac, ibd]
    U = real_harmonic_transform(orbs)
    # U[a, p]: real orbital a in terms of complex p.  Contract the FIRST index of U
    # with the tensor's complex index -- "pa,..." would transpose the transform and
    # silently corrupt the four-index tensor (caught by G2: imag = sym residual = 0.146).
    gr = np.einsum("ap,bq,cr,ds,pqrs->abcd", U.conj(), U.conj(), U, U, g,
                   optimize=True)
    imag = float(np.abs(gr.imag).max())
    return np.real(gr), imag, U


def one_body_real(S, h1, U):
    Sr = np.einsum("ap,bq,pq->ab", U.conj(), U, S, optimize=True)
    hr = np.einsum("ap,bq,pq->ab", U.conj(), U, h1, optimize=True)
    return np.real(Sr), np.real(hr), max(float(np.abs(Sr.imag).max()),
                                         float(np.abs(hr.imag).max()))


def fci(S, h1, g, n_elec):
    n_orb = S.shape[0]
    X = TC.lowdin(S)
    h1o = TC.transform_1(h1, X)
    go = TC.transform_2(g, X)
    nso = 2 * n_orb
    dets, didx = TC.make_dets(nso, n_elec)
    H = TC.build_H(dets, didx, TC.h_spin(h1o, nso),
                   TC.asym_from_phys(go, nso), nso)
    return float(eigh(H, eigvals_only=True)[0]), len(dets)


if __name__ == "__main__":
    out = {}
    k, Z, Lmax = 2.0, 2.0, 2
    ns, npp = 3, 1

    print("=" * 78)
    print("G1 -- s-only sector vs the independent s-only engine")
    print("=" * 78)
    orbs0, rad0, S0, h10, rpi0, RL0, RLs0, RLw0 = build(ns, 0, k, Z, Lmax=0)
    g0, im0, U0 = assemble(orbs0, rpi0, RL0, Lmax=0)
    r, wr = TC.make_grid(k, Ng=700)
    Sref, href, Rtab, W2 = TC.build_one_body(ns, r, wr, k, Z)
    Km = TC.build_kernels(r, 0.7, nx=64)
    gref, _, _ = TC.two_body(ns, Rtab, W2, Km)
    print(f"  S   dev {np.abs(S0 - Sref).max():.2e}")
    print(f"  h1  dev {np.abs(h10 - href).max():.2e}")
    print(f"  ERI dev {np.abs(g0 - gref).max():.2e}   (imag discarded {im0:.1e})")
    e_new, nd0 = fci(S0, h10, g0, 2)
    e_ref, _ = fci(Sref, href, gref, 2)
    print(f"  He FCI: new {e_new:.10f}   ref {e_ref:.10f}   dev {abs(e_new-e_ref):.2e}")
    out["G1"] = dict(S=float(np.abs(S0-Sref).max()), h1=float(np.abs(h10-href).max()),
                     eri=float(np.abs(g0-gref).max()), fci_dev=abs(e_new-e_ref))

    print()
    print("=" * 78)
    print("G2/G3/G4 -- s+p tensor: reality, symmetry, variational sanity, split exactness")
    print("=" * 78)
    orbs, rad, S, h1, rpi, RL, RLs, RLw = build(ns, npp, k, Z, Lmax=Lmax)
    g, imag, U = assemble(orbs, rpi, RL, Lmax=Lmax)
    Sr, hr, im1 = one_body_real(S, h1, U)
    gs, _, _ = assemble(orbs, rpi, RLs, Lmax=Lmax)
    gw, _, _ = assemble(orbs, rpi, RLw, Lmax=Lmax)
    sym = max(np.abs(g - g.transpose(2, 1, 0, 3)).max(),
              np.abs(g - g.transpose(0, 3, 2, 1)).max(),
              np.abs(g - g.transpose(1, 0, 3, 2)).max())
    print(f"  orbitals: {len(orbs)} spatial ({ns}s + {npp*3}p)")
    print(f"  max |imag| discarded: ERI {imag:.2e}, one-body {im1:.2e}")
    print(f"  8-fold permutational symmetry residual: {sym:.2e}")
    print(f"  split exactness |g - (gsep - gW)|: {np.abs(g - (gs - gw)).max():.2e}")
    e_sp, ndet = fci(Sr, hr, g, 2)
    print(f"  He FCI s+p = {e_sp:.10f}  ({ndet} dets);  s-only = {e_new:.10f}"
          f";  exact = -2.9037243770")
    print(f"  below s-only: {e_sp < e_new};  above exact: {e_sp > -2.9037243770}")
    out["G2G3G4"] = dict(imag=imag, sym=float(sym),
                         split=float(np.abs(g - (gs - gw)).max()),
                         e_sp=e_sp, e_s=e_new, ndet=ndet)

    print()
    print("=" * 78)
    print("THE PAYOFF CURVE -- truncate every W_L at rank m, assemble, FCI")
    print("=" * 78)
    e_full = e_sp
    print(f"{'rank m':>8}{'E(s+p)':>16}{'E - E_full (mHa)':>20}")
    rows = []
    for m in sorted({0, 1, 2, 3, 4, 6, 8, min(12, len(rpi)), len(rpi)}):
        RLw_t = {}
        for L in range(Lmax + 1):
            Wl = RLw[L]
            lam, V = np.linalg.eigh(Wl)
            o = np.argsort(-np.abs(lam))[:m]
            RLw_t[L] = (V[:, o] * lam[o]) @ V[:, o].T if m > 0 else np.zeros_like(Wl)
        gw_t, _, _ = assemble(orbs, rpi, RLw_t, Lmax=Lmax)
        e_t, _ = fci(Sr, hr, gs - gw_t, 2)
        rows.append((m, e_t, 1000 * (e_t - e_full)))
        print(f"{m:>8}{e_t:>16.10f}{1000*(e_t-e_full):>20.4f}")
    out["curve"] = rows
    out["n_radpairs"] = len(rpi)

    os.makedirs("debug/data", exist_ok=True)
    with io.open("debug/data/ee_split_sp_fci.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=float)
    print("\nwrote debug/data/ee_split_sp_fci.json")
