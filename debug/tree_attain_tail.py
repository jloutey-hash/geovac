"""tree_attain: the split of the integrand on [1, oo) for a general s-type Coulomb tree, checked numerically.

Hat normalisation: rho = K0, Phi W(x) = (1/x) int_0^x W s^2 ds + int_x^oo W s ds, H = x^2 prod_{w child of root} P_w.
For every non-root vertex u:  P_u = Q_u/x + E_u,  Q_u = int_0^oo W_u s^2 ds,
    E_u(x) = -(1/x) int_x^oo W_u(s) s (s - x) ds   (the exterior tail; theta(theta+1) E_u = -x^2 W_u).
"Full tail": W_u^full = K0 prod_{w child of u} E_w^full, E_u^full = tail of W_u^full; H^top = x^2 prod E_w^full.
Claims checked:
  (a) P_u - Q_u/x equals the tail formula (grid, each vertex);
  (b) e^{(n-1)x} H^top(x) agrees on 25 <= x <= 55 with the local-rule series (tests/test_paper63_local_rule.py Tail
      machinery, 7 terms), whose leading constant is (-1)^(n-1) (pi/2)^((n-1)/2) / prod_u k_u^2, k_u = |subtree(u)|;
  (c) J(eps) = int_1^oo I0((n-1-eps)x) H^top dx has the leading singular term the proof uses:
      C (2 pi (n-1))^(-1/2) Gamma(a+1) eps^(-a-1), a = (4-n)/2 (n odd or n = 4), fitted amplitude vs prediction.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp
from scipy.integrate import cumulative_simpson
from scipy.special import i0e, k0

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "debug"))
import test_paper63_local_rule as LR  # noqa: E402
from tree_attain_algebra import Tree  # noqa: E402

T0 = time.time()


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


U = np.linspace(np.log(1e-10), np.log(160.0), 800001)
R = np.exp(U)


def cum0(f):
    return cumulative_simpson(f * R, x=U, initial=0.0)


def tail(f):
    """int_x^oo f ds, summed from the right (keeps relative accuracy of exponentially small tails)."""
    g = (f * R)[::-1]
    return cumulative_simpson(g, x=-U[::-1], initial=0.0)[::-1]


def run(name: str, n: int, edges: List[Tuple[int, int]]) -> Dict:
    tree = Tree(n, edges)
    res: Dict = {"n": n, "edges": edges}
    sel = (R >= 1.0) & (R <= 30.0)

    # ---- (a) direct potentials and the tail formula
    Wd: Dict[int, np.ndarray] = {}
    Pd: Dict[int, np.ndarray] = {}
    worst = 0.0

    def W_direct(v: int) -> np.ndarray:
        out = k0(R)
        for w in tree.children[v]:
            Wd[w] = W_direct(w)
            Pd[w] = cum0(Wd[w] * R ** 2) / R + tail(Wd[w] * R)
            out = out * Pd[w]
        return out

    for w in tree.children[tree.root]:
        Wd[w] = W_direct(w)
        Pd[w] = cum0(Wd[w] * R ** 2) / R + tail(Wd[w] * R)
    for u in Wd:
        Q = cum0(Wd[u] * R ** 2)[-1]
        E = -(tail(Wd[u] * R ** 2) - R * tail(Wd[u] * R)) / R
        d = np.max(np.abs(Pd[u][sel] - Q / R[sel] - E[sel]) / np.abs(Pd[u][sel]))
        worst = max(worst, float(d))
    res["a_max_rel_defect_P_minus_Q_over_x_minus_tail"] = worst

    # ---- (b) full tails on the grid and from the local rule
    def full(v: int) -> np.ndarray:
        out = k0(R)
        for w in tree.children[v]:
            Wf = full(w)
            out = out * (-(tail(Wf * R ** 2) - R * tail(Wf * R)) / R)
        return out

    Htop = R ** 2
    for w in tree.children[tree.root]:
        Wf = full(w)
        Htop = Htop * (-(tail(Wf * R ** 2) - R * tail(Wf * R)) / R)

    def full_series(v: int) -> LR.Tail:
        out = LR.K_tail(0)
        for w in tree.children[v]:
            out = out * full_series(w).potential(0).times(1 / (4 * sp.pi))
        return out

    Hs = None
    for w in tree.children[tree.root]:
        e = full_series(w).potential(0).times(1 / (4 * sp.pi))
        Hs = e if Hs is None else Hs * e
    Hs = Hs.times(1, shift=2)
    assert Hs.k == n - 1 and Hs.a == sp.Rational(5 - n, 2), (Hs.k, Hs.a)
    Cpred = (-1) ** (n - 1) * (sp.pi / 2) ** sp.Rational(n - 1, 2) / sp.prod(
        [tree.size[u] ** 2 for u in range(n) if u != tree.root])
    assert sp.simplify(Hs.lead() - Cpred) == 0, (Hs.lead(), Cpred)
    res["C_local_rule"] = str(Cpred)
    res["C_local_rule_float"] = float(Cpred)
    res["series_monic_coeffs"] = [str(c) for c in Hs.monic()[:4]]
    pref, coef = float(Hs.pref), [float(c) for c in Hs.c]

    def series_scaled(x: np.ndarray) -> np.ndarray:  # e^{(n-1)x} H^top from the series
        return pref * x ** float(Hs.a) * sum(c * x ** (-j) for j, c in enumerate(coef))

    chk = {}
    for x0 in (25.0, 40.0, 55.0):
        i = int(np.searchsorted(R, x0))
        true = Htop[i] * math.exp((n - 1) * R[i])
        ser = float(series_scaled(np.array([R[i]]))[0])
        lead_only = pref * R[i] ** float(Hs.a) * coef[0]
        chk[str(x0)] = {"true_over_series_minus_1": true / ser - 1, "true_over_leading_minus_1": true / lead_only - 1}
    res["b_Htop_vs_local_rule"] = chk

    # ---- (c) the branch integral
    XC = 50.0
    a = sp.Rational(4 - n, 2)

    def J(eps: float) -> float:
        y = np.linspace(0.0, np.log(1.0 + 60.0 / eps + XC), 400001)
        x = np.exp(y)
        scaled = np.where(x <= XC, np.interp(np.minimum(x, XC), R, Htop * np.exp((n - 1) * R) if n <= 5 else Htop)
                          * (1.0 if n <= 5 else 0.0), 0.0)
        if n > 5:  # avoid overflow: interpolate log
            lg = np.interp(np.minimum(x, XC), R[R >= 0.5], np.log(np.abs(Htop[R >= 0.5])) + (n - 1) * R[R >= 0.5])
            scaled = np.sign(pref) * np.exp(lg)
        scaled = np.where(x <= XC, scaled, series_scaled(x))
        f = i0e((n - 1 - eps) * x) * np.exp(-eps * x) * scaled
        return float(np.trapezoid(f * x, y))

    eps_list = [10 ** (-k / 2) for k in range(2, 9)]
    Jv = np.array([J(e) for e in eps_list])
    e = np.array(eps_list)
    lead_amp = float(Cpred / sp.sqrt(2 * sp.pi * (n - 1)))
    if n % 2 == 1:
        ex = float(-a - 1)
        basis = np.stack([e ** ex, e ** (ex + 1), np.ones_like(e), e ** (ex + 2), e], axis=1)
        pred = lead_amp * math.gamma(float(a) + 1)
        what = f"coefficient of eps^({-a - 1})"
    elif n == 4:
        basis = np.stack([1 / e, np.log(e), np.ones_like(e), e * np.log(e), e], axis=1)
        pred = lead_amp
        what = "residue (coefficient of 1/eps)"
    else:  # n = 6: -log eps
        basis = np.stack([np.log(e), np.ones_like(e), e * np.log(e), e], axis=1)
        pred = -lead_amp
        what = "coefficient of log eps"
    sol = np.linalg.lstsq(basis, Jv, rcond=None)[0]
    res["c_branch"] = {"what": what, "fitted": float(sol[0]), "predicted": pred, "ratio_minus_1": float(sol[0] / pred - 1),
                       "eps": eps_list, "J": Jv.tolist()}
    if n == 4:
        res["c_branch"]["fitted_log_coefficient"] = float(sol[1])
    log(f"{name}: (a) {worst:.1e}; (b) {chk}; (c) {what}: fitted {sol[0]:.6g}, predicted {pred:.6g}")
    return res


TREES = {
    "star4": (4, [(0, 1), (0, 2), (0, 3)]),
    "path3_root_end": (3, [(0, 1), (1, 2)]),
    "path4_root_end": (4, [(0, 1), (1, 2), (2, 3)]),
    "path4_root_inner": (4, [(0, 1), (0, 2), (2, 3)]),
    "star5": (5, [(0, 1), (0, 2), (0, 3), (0, 4)]),
    "path5_root_end": (5, [(0, 1), (1, 2), (2, 3), (3, 4)]),
    "path5_root_centre": (5, [(0, 1), (1, 2), (0, 3), (3, 4)]),
    "fork5_root_end": (5, [(0, 1), (1, 2), (1, 3), (3, 4)]),
    "fork5_root_branch": (5, [(1, 0), (0, 2), (0, 3), (3, 4)]),
    "spider6_root_end": (6, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)]),
}

if __name__ == "__main__":
    out = {}
    for name in (sys.argv[1:] or list(TREES)):
        n, edges = TREES[name]
        out[name] = run(name, n, edges)
        (ROOT / "debug" / "data" / "tree_attain_tail.json").write_text(json.dumps(out, indent=1))
    log("done")
