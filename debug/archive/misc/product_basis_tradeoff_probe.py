"""Product-basis trade-off probe (scoping sprint, 2026-09-26). float64, ~1 min.

Uses the stage-A product basis of debug/prolate_separation_constant.py:
  f_a(xi) = (xi^2-1)^{m/2} e^{-alpha(xi-1)} L_a^{(m)}(2 alpha (xi-1)),  g_b(eta) = P_l^m(eta)
(a FULL tensor product, so Atkinson exactness holds for any separable one-body operator).

P1  Core resolution in a pure product basis.  The Li 1s  e^{-zeta r_A} = e^{-(zR/2)(xi-1)} e^{-(zR/2)(1+eta)}
    is ONE product function, but its xi-rate zR/2 and its eta-exponential are outside a
    single-alpha Laguerre x Legendre span.  Measures (i) the lowest eigenvalue of the bare
    one-centre h(Z_A=3, Z_B=0) (exact -4.5 Ha, R-independent) and (ii) the fraction of the
    zeta=2.6875 core norm outside the span, vs R, basis size and alpha.  The R-spread at
    fixed basis is the product-basis analogue of increment-7's R-dependent core resolution.

P2  Screened separation constant.  A valence one-body operator for LiH:
      h_eff = T - 1/r_A - 1/r_B            (separable: Lambda(1,1) commutes exactly)
            + dV(r_A)                       (Hartree residual of Li 1s^2: -2 e^{-2 zeta r_A}(zeta + 1/r_A))
            + sigma |c><c|                  (level-shift core projector; c is a product function)
    Core exchange is omitted (stated).  Measures ||[h_eff, Lambda(1,1)]||_rel, and for the
    lowest valence eigenstates: <Lambda>, leakage ||(Lambda-<Lambda>)psi||, and the rank-1
    (x (x) y) ratio sigma2/sigma1.  Compare: stage C bare-nucleus MOs on the non-product
    Route C set had leakage 1.8-4.7 and ||[h1,Lambda]||_rel = 2.2e-2.
Output: stdout + debug/data/product_basis_tradeoff_probe.log
"""
from __future__ import annotations

import os
import sys
from typing import Dict, List

import numpy as np
from scipy.linalg import cholesky, eigh, solve_triangular
from scipy.special import eval_genlaguerre, roots_laguerre, roots_legendre

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from prolate_separation_constant import OneE, _assoc_legendre_norm  # noqa: E402

ZETA = 2.6875
LINES: List[str] = []


def log(s: str) -> None:
    print(s, flush=True)
    LINES.append(s)


# ------------------------------------------------------------------ weighted 1D blocks
def _xi_norm(Nx: int, m: int, alpha: float) -> np.ndarray:
    nq = Nx + m + 6
    t, w = roots_laguerre(nq)
    x = t / (2 * alpha)
    u = x * (x + 2.0)
    Lv = np.array([eval_genlaguerre(a, m, t) for a in range(Nx)])
    r2d = (Lv ** 2 * (u ** m * (x + 1) ** 2 * w / (2 * alpha))).sum(axis=1)
    return 1.0 / np.sqrt(r2d)


def xi_weighted(Nx: int, m: int, alpha: float, gam: float, k: int) -> np.ndarray:
    """<f_a| xi^k e^{-gam (xi-1)} |f_b>  (same normalization as xi_blocks)."""
    s = _xi_norm(Nx, m, alpha)
    rate = 2 * alpha + gam
    nq = Nx + m + k + 8
    t, w = roots_laguerre(nq)
    x = t / rate
    u = x * (x + 2.0)
    Lv = np.array([eval_genlaguerre(a, m, 2 * alpha * x) for a in range(Nx)])
    wt = w / rate * u ** m * (1 + x) ** k
    M = (Lv * wt) @ Lv.T
    return M * np.outer(s, s)


def xi_vec(Nx: int, alpha: float, gam: float, k: int) -> np.ndarray:
    """<f_a| xi^k e^{-gam (xi-1)}>  for m = 0."""
    s = _xi_norm(Nx, 0, alpha)
    rate = alpha + gam
    t, w = roots_laguerre(Nx + k + 8)
    x = t / rate
    Lv = np.array([eval_genlaguerre(a, 0, 2 * alpha * x) for a in range(Nx)])
    return (Lv * (w / rate * (1 + x) ** k)).sum(axis=1) * s


_ETA_T, _ETA_W = roots_legendre(400)


def eta_weighted(Ne: int, m: int, beta: float, q: int) -> np.ndarray:
    """<g_b| eta^q e^{-beta (1+eta)} |g_c>."""
    P = _assoc_legendre_norm(m + Ne - 1, m, _ETA_T)
    wt = _ETA_W * _ETA_T ** q * np.exp(-beta * (1 + _ETA_T))
    return (P * wt) @ P.T


def eta_vec(Ne: int, beta: float, q: int) -> np.ndarray:
    P = _assoc_legendre_norm(Ne - 1, 0, _ETA_T)
    return (P * (_ETA_W * _ETA_T ** q * np.exp(-beta * (1 + _ETA_T)))).sum(axis=1)



# ------------------------------------------------------------------ W-frame helpers
def wframe(W: np.ndarray, *mats: np.ndarray):
    L = cholesky(W, lower=True)
    out = []
    for A in mats:
        Ah = solve_triangular(L, solve_triangular(L, A, lower=True).T, lower=True).T
        out.append(0.5 * (Ah + Ah.T))
    return L, out


def rel_comm(A: np.ndarray, B: np.ndarray) -> float:
    C = A @ B - B @ A
    return float(np.linalg.norm(C, 2) / (np.linalg.norm(A, 2) * np.linalg.norm(B, 2)))


# ------------------------------------------------------------------ P1
def p1() -> None:
    log("=" * 78)
    log("P1  core resolution in a single-alpha product basis (m=0)")
    log("    bare h(Z_A=3,Z_B=0): exact E_1s = -4.5 Ha, R-independent;  core-norm")
    log("    outside span for e^{-zeta r_A}, zeta=2.6875")
    log("=" * 78)
    Rs = [2.6, 3.015, 3.45]
    for (Nx, Ne) in [(10, 10), (20, 20), (30, 40)]:
        for amode in ["1.2", "R/2", "ZR/2"]:
            Es, outs = [], []
            for R in Rs:
                al = {"1.2": 1.2, "R/2": 0.5 * R, "ZR/2": 1.5 * R}[amode]
                op = OneE(R, 3.0, 0.0, 0, Nx, Ne, al)
                E = eigh(op.Hm, op.W, eigvals_only=True)[0]
                Es.append(E)
                g = ZETA * R / 2
                b = np.kron(xi_vec(Nx, al, g, 2), eta_vec(Ne, g, 0)) - \
                    np.kron(xi_vec(Nx, al, g, 0), eta_vec(Ne, g, 2))
                # exact core norm in the (xi^2 - eta^2) measure
                t, w = roots_laguerre(60)
                x = t / (2 * g)
                X0 = np.sum(w / (2 * g)); X2 = np.sum(w / (2 * g) * (1 + x) ** 2)
                E0 = np.sum(_ETA_W * np.exp(-2 * g * (1 + _ETA_T)))
                E2 = np.sum(_ETA_W * _ETA_T ** 2 * np.exp(-2 * g * (1 + _ETA_T)))
                Nc = X2 * E0 - X0 * E2
                outs.append(1.0 - b @ np.linalg.solve(op.W, b) / Nc)
            log(f"  Nx={Nx:2d} Ne={Ne:2d} alpha={amode:5s}  E_1s-(-4.5) = "
                + "  ".join(f"{e + 4.5:+.2e}" for e in Es)
                + f"  | R-spread {max(Es) - min(Es):.2e} Ha"
                + "  | core outside span " + " ".join(f"{o:.1e}" for o in outs))


# ------------------------------------------------------------------ P2
def build_eff(R: float, m: int, Nx: int, Ne: int, al: float, with_dV: bool,
              with_P: bool, sigma: float = 20.0) -> Dict[str, np.ndarray]:
    op11 = OneE(R, 1.0, 1.0, m, Nx, Ne, al)          # separable part + Lambda(1,1)
    op31 = OneE(R, 3.0, 1.0, m, Nx, Ne, al)          # bare: Lambda(3,1)
    H = op11.Hm.copy()
    g2 = ZETA * R                                     # 2 zeta r_A = zR (xi + eta)
    if with_dV:
        X0, X1, X2 = (xi_weighted(Nx, m, al, g2, k) for k in (0, 1, 2))
        E0, E1, E2 = (eta_weighted(Ne, m, g2, q) for q in (0, 1, 2))
        D = -2.0 * (ZETA * (np.kron(X2, E0) - np.kron(X0, E2))
                    + (2.0 / R) * (np.kron(X1, E0) - np.kron(X0, E1)))
        H = H + D
    if with_P and m == 0:
        g = ZETA * R / 2
        b = np.kron(xi_vec(Nx, al, g, 2), eta_vec(Ne, g, 0)) - \
            np.kron(xi_vec(Nx, al, g, 0), eta_vec(Ne, g, 2))
        Nc = b @ np.linalg.solve(op11.W, b)          # in-span norm (self-consistent)
        H = H + sigma * np.outer(b, b) / Nc
    return dict(W=op11.W, H=H, N11=op11.N, N31=op31.N, H31=op31.Hm)


def p2() -> None:
    log("=" * 78)
    log("P2  screened Lambda: h_eff = T - 1/r_A - 1/r_B + dV_Hartree(Li 1s^2) + sigma|c><c|")
    log("    (core exchange omitted).  Lambda(1,1) vs bare Lambda(3,1).")
    log("=" * 78)
    Nx, Ne = 30, 40
    for R in [3.015, 2.6, 3.45]:
        al = 1.2
        for m in (0, 1):
            for tag, dV, P in [("sep only", False, False), ("+dV", True, False),
                               ("+dV+P", True, True)]:
                if m == 1 and P:
                    continue
                d = build_eff(R, m, Nx, Ne, al, dV, P)
                L, (Hh, N11h, N31h) = wframe(d["W"], d["H"], d["N11"], d["N31"])
                c11 = rel_comm(Hh, N11h)
                c31 = rel_comm(Hh, N31h)
                E, V = eigh(Hh)
                rows = []
                for i in range(4):
                    v = V[:, i]
                    lam = v @ N11h @ v
                    leak = np.linalg.norm(N11h @ v - lam * v)
                    lam3 = v @ N31h @ v
                    leak3 = np.linalg.norm(N31h @ v - lam3 * v)
                    cf = solve_triangular(L.T, v, lower=False).reshape(Nx, Ne)
                    s = np.linalg.svd(cf, compute_uv=False)
                    rows.append(f"E={E[i]:+.4f} L11={lam:+.3f} leak11={leak:.2e} "
                                f"leak31={leak3:.2e} r1={s[1] / s[0]:.1e}")
                log(f"  R={R} m={m} {tag:8s} ||[h,L11]||rel={c11:.1e}  ||[h,L31]||rel={c31:.1e}")
                for r in rows:
                    log("      " + r)


def p2_perturbation_size() -> None:
    """How big is the non-separable part, measured where the valence lives?"""
    log("=" * 78)
    log("P2b size of the non-separable pieces on the lowest valence state of h_eff (m=0)")
    log("=" * 78)
    Nx, Ne, al = 30, 40, 1.2
    for R in [2.6, 3.015, 3.45]:
        d_full = build_eff(R, 0, Nx, Ne, al, True, True)
        d_sep = build_eff(R, 0, Nx, Ne, al, False, False)
        d_dv = build_eff(R, 0, Nx, Ne, al, True, False)
        L, (Hf, Hs, Hd) = wframe(d_full["W"], d_full["H"], d_sep["H"], d_dv["H"])
        E, V = eigh(Hf)
        v = V[:, 0]
        dV = Hd - Hs
        Pm = Hf - Hd
        Es, Vs = eigh(Hs)
        ov = np.abs(Vs.T @ v) ** 2
        log(f"  R={R}: E_val={E[0]:+.5f}  <dV>={v @ dV @ v:+.4f}  <P>={v @ Pm @ v:+.2e}  "
            f"||dV v||={np.linalg.norm(dV @ v):.3f}  max|<sep-eigvec|v>|^2={ov.max():.3f} "
            f"(sep state #{int(ov.argmax())}, E_sep={Es[int(ov.argmax())]:+.4f})")


if __name__ == "__main__":
    p1()
    p2()
    p2_perturbation_size()
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    with open(os.path.join(HERE, "data", "product_basis_tradeoff_probe.log"), "w") as f:
        f.write("\n".join(LINES))
