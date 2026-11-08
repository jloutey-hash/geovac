"""Q3 of the STO -> shared-k Sturmian evaluator scoping: H2 and LiH head-to-head (debug only).

Three evaluators of the SAME s-type STO basis:
  (G)  noci_engine production: McMurchie-Davidson over Gaussian fits of the STO shapes;
  (X)  exact one-electron integrals by the Mulliken auxiliary route (A_m(alpha), B_m(beta) in
       prolate spheroidal coordinates, mpmath, no quadrature) applied DIRECTLY to the STOs;
  (S)  the same exact Mulliken route applied to the shared-k Sturmian EXPANSION of each STO
       (Avery coefficients, truncated at N terms) -- the candidate replacement evaluator.
ERIs: exact closed forms for H2 1s (Slater/Ruedenberg/Sugiura); for LiH the four-index tensor
comes from the Gaussian route at increasing fit depth (the exchange class (AB|AB) has no
exact engine with per-orbital exponents in-repo).

Writes debug/data/sto_sturmian_eval_headtohead.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from functools import lru_cache
from math import comb

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from geovac import noci_engine as E  # noqa: E402
from sto_sturmian_eval_convergence import coeffs as st_coeffs  # noqa: E402

mp.mp.dps = 60
OUT = os.path.join(HERE, "data", "sto_sturmian_eval_headtohead.json")


# ------------------------------------------------------------- Mulliken auxiliaries

@lru_cache(maxsize=None)
def A_aux(m: int, alpha):
    alpha = mp.mpf(alpha)
    return mp.exp(-alpha) * mp.fsum(mp.factorial(m) / (mp.factorial(k) * alpha ** (m - k + 1))
                                    for k in range(m + 1))


@lru_cache(maxsize=None)
def B_aux(m: int, beta):
    beta = mp.mpf(beta)
    if beta == 0:
        return mp.mpf(2) / (m + 1) if m % 2 == 0 else mp.mpf(0)
    tot, k, term_mag = mp.mpf(0), 0, mp.mpf(1)
    while True:
        if (m + k) % 2 == 0:
            t = (-beta) ** k / mp.factorial(k) * mp.mpf(2) / (m + k + 1)
            tot += t
            if k > abs(beta) * 3 + 20 and abs(t) < mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        k += 1
    return tot


def two_center(p: int, a, q: int, b, R):
    """int r_a^p r_b^q e^{-a r_a - b r_b} d^3r, p, q >= -1 (A, B separated by R)."""
    a, b, R = mp.mpf(a), mp.mpf(b), mp.mpf(R)
    al, be = (a + b) * R / 2, (a - b) * R / 2
    P, Q = p + 1, q + 1
    tot = mp.mpf(0)
    for i in range(P + 1):
        for j in range(Q + 1):
            tot += comb(P, i) * comb(Q, j) * (-1) ** j * A_aux(P + Q - i - j, al) * B_aux(i + j, be)
    return 2 * mp.pi * (R / 2) ** (3 + p + q) * tot


def one_center(p: int, a):
    """int r^p e^{-a r} d^3r."""
    return 4 * mp.pi * mp.factorial(p + 2) / mp.mpf(a) ** (p + 3)


# ------------------------------------------------------------- poly-exp radial functions
# A radial function on a centre is a list of (power, coeff, rate) terms: sum c r^p e^{-rate r},
# times Y00 = 1/sqrt(4 pi).

def sto_terms(n: int, zeta):
    zeta = mp.mpf(zeta)
    N = mp.sqrt((2 * zeta) ** (2 * n + 1) / mp.factorial(2 * n))
    return [(n - 1, N, zeta)]


def lap_terms(terms):
    """-1/2 nabla^2 acting on sum c r^p e^{-z r} Y00:  -1/2 [f'' + 2 f'/r]."""
    out = {}

    def add(p, c, z):
        out[(p, z)] = out.get((p, z), 0) + c
    for p, c, z in terms:
        # f = r^p e^{-zr}: f'' + 2f'/r = p(p+1) r^{p-2} - 2 z (p+1) r^{p-1} + z^2 r^p
        if p * (p + 1) != 0:
            add(p - 2, -c * p * (p + 1) / 2, z)
        add(p - 1, c * z * (p + 1), z)
        add(p, -c * z * z / 2, z)
    return [(p, c, z) for (p, z), c in sorted(out.items())]


def pair_integral(ti, ci, tj, cj, R, w_power=None, w_center=None):
    """<f_i| W |f_j> with f on centres ci, cj (0/1, B at distance R); W = 1 or 1/r_{w_center}.
    Angular factor Y00^2 = 1/(4 pi)."""
    tot = mp.mpf(0)
    for (p1, c1, z1) in ti:
        for (p2, c2, z2) in tj:
            if ci == cj:
                pw = p1 + p2
                if w_power is not None and w_center == ci:
                    tot += c1 * c2 * one_center(pw - 1, z1 + z2)
                elif w_power is not None:          # 1/r to the other centre: two-centre, b=0
                    tot += c1 * c2 * two_center(pw, z1 + z2, -1, 0, R)
                else:
                    tot += c1 * c2 * one_center(pw, z1 + z2)
            else:
                # put i on A=ci; j on the other centre
                pa, pb = p1, p2
                if w_power is not None:
                    if w_center == ci:
                        pa -= 1
                    else:
                        pb -= 1
                if pa < -1 or pb < -1:
                    raise ValueError("power below -1")
                tot += c1 * c2 * two_center(pa, z1, pb, z2, R)
    return tot / (4 * mp.pi)


def h_matrix(funcs, R, Zs):
    """funcs: list of (terms, centre). Returns S, H as mp matrices (centres at 0 and R)."""
    m = len(funcs)
    S = mp.matrix(m, m)
    H = mp.matrix(m, m)
    for i in range(m):
        for j in range(i, m):
            ti, ci = funcs[i]
            tj, cj = funcs[j]
            s = pair_integral(ti, ci, tj, cj, R)
            t = pair_integral(ti, ci, lap_terms(tj), cj, R)
            v = mp.mpf(0)
            for C, Z in enumerate(Zs):
                v -= Z * pair_integral(ti, ci, tj, cj, R, w_power=-1, w_center=C)
            S[i, j] = S[j, i] = s
            H[i, j] = H[j, i] = t + v
    return S, H


def sturmian_terms(n_sto: int, zeta, k, nterms: int):
    """Shared-k Sturmian expansion of the NORMALISED STO r^{n-1}e^{-zeta r}, truncated."""
    k = mp.mpf(k)
    zeta = mp.mpf(zeta)
    rho = zeta / k
    c = st_coeffs(n_sto, 0, rho, nmax=nterms)     # expansion of s^{n-1} e^{-rho s}, s = k r
    Nsto = mp.sqrt((2 * zeta) ** (2 * n_sto + 1) / mp.factorial(2 * n_sto))
    # r^{n-1} e^{-zeta r} = k^{-(n-1)} s^{n-1} e^{-rho s};  S^{(1)}(s) = k^{-3/2} S^{(k)}(r)
    pref = Nsto / k ** (n_sto - 1) / k ** mp.mpf(1.5)
    terms = {}
    for idx, n in enumerate(range(1, nterms + 1)):
        # L2-normalised Sturmian at k, l=0:  N (2kr)^0 e^{-kr} L_{n-1}^1(2kr)
        Nst = mp.sqrt((2 * k) ** 3 / (2 * mp.mpf(n) * n))   # (n-l-1)!/(n+l)! = 1/n for l=0
        for j in range(n):
            lag = (-1) ** j * comb(n, n - 1 - j) / mp.factorial(j)
            terms[j] = terms.get(j, 0) + pref * c[idx] * Nst * lag * (2 * k) ** j
    return [(p, cf, k) for p, cf in sorted(terms.items())]


# ------------------------------------------------------------- energies

def energies_from(s, h, g, n_elec, vnn, dets=None):
    s, h, g = np.asarray(s, float), np.asarray(h, float), np.asarray(g, float)
    x = E.lowdin_orbitals(s)
    ht, gt = E.transform_integrals(x, h, g)
    out = {"FCI": float(E.fci_ground(ht, gt, n_elec) + vnn)}
    if dets:
        for lab, d in dets.items():
            e, _ = E.noci_ground_gensc(d, s, h, g)
            out[lab] = float(e + vnn)
    return out


def tomat(M):
    return np.array([[float(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])


# ------------------------------------------------------------- H2 exact ERIs (zeta-scaled)

def h2_exact_eri(zeta, R):
    z, R = mp.mpf(zeta), mp.mpf(R)
    rho = z * R
    aaaa = 5 * z / 8
    aabb = (1 / R) * (1 - (1 + 11 * rho / 8 + 3 * rho ** 2 / 4 + rho ** 3 / 6) * mp.exp(-2 * rho))
    aaab = z * (mp.exp(-rho) * (rho + mp.mpf(1) / 8 + 5 / (16 * rho))
                + mp.exp(-3 * rho) * (-mp.mpf(1) / 8 - 5 / (16 * rho)))
    S = mp.exp(-rho) * (1 + rho + rho ** 2 / 3)
    Sp = mp.exp(rho) * (1 - rho + rho ** 2 / 3)
    abab = (z / 5) * (-mp.exp(-2 * rho) * (-mp.mpf(25) / 8 + 23 * rho / 4 + 3 * rho ** 2 + rho ** 3 / 3)
                      + (6 / rho) * (S ** 2 * (mp.euler + mp.log(rho)) + Sp ** 2 * mp.ei(-4 * rho)
                                     - 2 * S * Sp * mp.ei(-2 * rho)))
    g = np.zeros((2, 2, 2, 2))
    vals = {}
    for p in range(2):
        for q in range(2):
            for r in range(2):
                for t in range(2):
                    idx = (p, q, r, t)
                    if p == q and r == t:
                        v = aaaa if p == r else aabb
                    elif p != q and r != t:
                        v = abab
                    else:
                        v = aaab
                    g[idx] = float(v)
    vals.update(aaaa=float(aaaa), aabb=float(aabb), aaab=float(aaab), abab=float(abab))
    return g, vals


def gaussian_set(shapes_kind, specs, nuclei_pos, Zs, ng):
    """specs: list of (centre_idx, n, zeta). shapes by fit depth ng, or 'prod' for STO-6G/3G."""
    orbs = []
    for ci, n, z in specs:
        if ng == "prod":
            a, d = E.STO6G_1S if n == 1 else E.STO3G_2S
        else:
            a, d = shapes_kind[(n, ng)]
        orbs.append(E.BasisFn(nuclei_pos[ci], (0, 0, 0), a * z ** 2, d))
    return E.integral_set_md(orbs, [(nuclei_pos[i], Zs[i]) for i in range(len(Zs))])


def main() -> None:
    res = {}
    t0 = time.time()
    depths = [4, 6, 8, 10, 12]
    shapes = {}
    for n in (1, 2):
        for ng in depths:
            a, d, q = E.fit_sto_shape(0, n, n_gauss=ng)
            shapes[(n, ng)] = (a, d)
            res.setdefault("fit_quality", {})[f"{n}s_{ng}"] = q
    print("fits done", time.time() - t0)

    # ---------------- H2, R = 1.4, zeta = 1 -----------------------------------------
    R = 1.4
    pos = [np.zeros(3), np.array([0, 0, R])]
    specs = [(0, 1, 1.0), (1, 1, 1.0)]
    funcs = [(sto_terms(1, 1.0), 0), (sto_terms(1, 1.0), 1)]
    Sx, Hx = h_matrix(funcs, R, [1.0, 1.0])
    gx, gvals = h2_exact_eri(1.0, R)
    Sx, Hx = tomat(Sx), tomat(Hx)
    h2 = {"exact": energies_from(Sx, Hx, gx, 2, 1 / R), "eri_exact": gvals}
    rows = []
    for ng in depths + ["prod"]:
        s, h, g = gaussian_set(shapes, specs, pos, [1.0, 1.0], ng)
        eG = energies_from(s, h, g, 2, 1 / R)["FCI"]
        e_hx = energies_from(Sx, Hx, g, 2, 1 / R)["FCI"]     # exact h, Gaussian g
        e_gx = energies_from(s, h, gx, 2, 1 / R)["FCI"]      # Gaussian h, exact g
        rows.append({"ng": ng, "E_gauss": eG, "dE": eG - h2["exact"]["FCI"],
                     "dE_from_h": e_gx - h2["exact"]["FCI"], "dE_from_g": e_hx - h2["exact"]["FCI"],
                     "max_dS": float(np.max(abs(s - Sx))), "max_dh": float(np.max(abs(h - Hx))),
                     "max_dg": float(np.max(abs(g - gx)))})
        print("H2", rows[-1])
    h2["gauss"] = rows
    res["H2"] = h2

    # ---------------- LiH, R = 3.25, Li 1s 2.69 / 2s 0.65, H 1s 1.0 ------------------
    R = 3.25
    pos = [np.zeros(3), np.array([0, 0, R])]
    Zs = [3.0, 1.0]
    specs = [(0, 1, 2.69), (0, 2, 0.65), (1, 1, 1.0)]
    UP, DN = 0, 1
    core = [(0, UP), (0, DN)]
    dets = {"NOCI2": [core + [(1, UP), (2, DN)], core + [(2, UP), (1, DN)]],
            "NOCI3": [core + [(1, UP), (2, DN)], core + [(2, UP), (1, DN)], core + [(2, UP), (2, DN)]]}
    funcs = [(sto_terms(n, z), ci) for ci, n, z in specs]
    t1 = time.time()
    Sx, Hx = h_matrix(funcs, R, Zs)
    Sx, Hx = tomat(Sx), tomat(Hx)
    print("LiH exact h", time.time() - t1)
    # Gaussian sequences
    gauss = {}
    for ng in depths + ["prod"]:
        s, h, g = gaussian_set(shapes, specs, pos, Zs, ng)
        gauss[ng] = (s, h, g)
    g_ref = gauss[12][2]
    lih = {"exact_h": {"S": Sx.tolist(), "H": Hx.tolist()},
           "E_exact_h_g12": energies_from(Sx, Hx, g_ref, 4, 3 / R, dets)}
    rows = []
    for ng, (s, h, g) in gauss.items():
        eG = energies_from(s, h, g, 4, 3 / R, dets)
        e_hx = energies_from(Sx, Hx, g, 4, 3 / R, dets)
        rows.append({"ng": ng, "E": eG, "E_exact_h": e_hx,
                     "max_dS": float(np.max(abs(s - Sx))), "max_dh": float(np.max(abs(h - Hx))),
                     "max_dg_vs_g12": float(np.max(abs(g - g_ref)))})
        print("LiH", ng, rows[-1]["max_dS"], rows[-1]["max_dh"], rows[-1]["max_dg_vs_g12"],
              {k: v - lih["E_exact_h_g12"][k] for k, v in eG.items()})
    lih["gauss"] = rows

    # ---------------- Sturmian-expanded evaluator (S), global k --------------------------
    stu = []
    for k in (1.261, 1.0):
        for nt in (4, 8, 12, 16, 20, 26, 32):
            t1 = time.time()
            fs = [(sturmian_terms(n, z, k, nt), ci) for ci, n, z in specs]
            Ss, Hs = h_matrix(fs, R, Zs)
            Ss, Hs = tomat(Ss), tomat(Hs)
            e = energies_from(Ss, Hs, g_ref, 4, 3 / R, dets)
            row = {"k": k, "nterms": nt, "max_dS": float(np.max(abs(Ss - Sx))),
                   "max_dh": float(np.max(abs(Hs - Hx))),
                   "dE_FCI": e["FCI"] - lih["E_exact_h_g12"]["FCI"],
                   "dE_NOCI3": e["NOCI3"] - lih["E_exact_h_g12"]["NOCI3"], "sec": time.time() - t1}
            stu.append(row)
            print("Sturm", row)
    lih["sturmian_expanded"] = stu
    res["LiH"] = lih
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print("wrote", OUT, time.time() - t0)


if __name__ == "__main__":
    main()
