"""tree_attain_l, step (ii): the top term of G for trees with l > 0 densities (any root).

G^top(x) = x^2 < P_0(-d) prod_{w child of root} E_w(x, omega) >,  E_w the exterior tail of the subtree of w with
every potential replaced by its tail (hat normalisation, rho = P(d) K0, channel-L tail
E_L(x) = (1/(2L+1)) int_x^oo W_L(s) [x^L s^(1-L) - x^(-L-1) s^(L+2)] ds).
Computed by nested Gauss-Laguerre in exponentially scaled variables with channels projected by mu-quadrature
(closed-form densities; no model, no local rule).  Root derivative: channel 0 of P_0(-d) Psi via
  s: Psi_0;  pz: -(1/3) r^-1 (theta+2) Psi_1;  h2: (1/5) r^-2 (theta+1)(theta+3) Psi_2;  lap: r^-2 theta(theta+1) Psi_0,
with theta E_L = (1/(2L+1)) [L x^L A + (L+1) x^(-L-1) B] and theta^2 E_L = -theta E_L + L(L+1) E_L - x^2 W_L.

N4: the constant C in G^top = C x^((5-n)/2) e^(-(n-1)x) (1 + O(1/x)) against the prediction
  C = (-1)^(n-1) (pi/2)^((n-1)/2) / prod_{u != root} k_u^2  x  < P_0((n-1) omega) prod_{v != root} P_v(-omega) >.
N5: the branch at t = -(n-1): J(eps) = int_1^oo I0((n-1-eps) x) G^top dx, leading singular coefficient
  C Gamma(b+1) / sqrt(2 pi (n-1)),  b = (4-n)/2, fitted on eps = 1e-1 .. 1e-3.
Usage: python debug/tree_attain_l_tail.py
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict

import numpy as np
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_legendre, i0e, k0e, k1e

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tree_attain_l_model as M  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
T0 = time.time()
GLX, GLW = laggauss(48)
MU, MW = leggauss(16)
LMAX = 12


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def dens_scaled(name: str, r, mu):
    """e^r P(d) K0(|r|), closed forms (Bessel recurrences)."""
    K0, K1 = k0e(r), k1e(r)
    dz = -mu * K1
    dz2 = mu ** 2 * K0 + (2 * mu ** 2 - 1) * K1 / r
    lap = K0 - K1 / r
    return {"s": K0 + 0 * mu, "pz": dz, "dz2": dz2, "lap": lap, "h2": 1.5 * dz2 - 0.5 * lap,
            "cancel": 2.25 * lap - 3.75 * dz2}[name]


class TailTree:
    def __init__(self, n, edges, root, dens):
        self.t = M.LTree(n, edges, root, {v: M.Pdict(dens[v]) for v in range(n)})
        self.dens, self.n, self.root = dens, n, root

    def k(self, u):
        return self.t.size[u]

    def wsc(self, u, s):
        """e^{k_u s} W_u(s, mu), shape s.shape + (len MU,)."""
        out = dens_scaled(self.dens[u], s[..., None], MU)
        for c in self.t.children[u]:
            Ec, _, _ = self.esc(c, s)
            out = out * Ec
        return out

    def esc(self, u, x):
        """e^{k x} E_u(x, mu), e^{k x} theta E_u(x, mu), e^{k x} theta^2 E_u(x, mu)."""
        k = self.k(u)
        v = GLX / k
        s = x[..., None] + v                                       # x.shape + (nodes,)
        Wv = self.wsc(u, s)                                        # x.shape + (nodes, mu)
        Wx = self.wsc(u, x)                                        # x.shape + (mu,)
        E = np.zeros(x.shape + (len(MU),))
        tE = np.zeros_like(E)
        t2E = np.zeros_like(E)
        for L in range(LMAX + 1):
            PL = eval_legendre(L, MU)
            WL = (2 * L + 1) / 2 * np.sum(Wv * (MW * PL), axis=-1)       # x.shape + (nodes,)
            WLx = (2 * L + 1) / 2 * np.sum(Wx * (MW * PL), axis=-1)
            if np.max(np.abs(WL)) == 0.0:
                continue
            xx = x[..., None]
            # A = int_x^oo s^(1-L) W, B = int_x^oo s^(L+2) W (scaled by e^{kx}); E = (x^L A - x^(-L-1) B)/(2L+1)
            # E written without cancellation: -(1/(2L+1)) int W x^(-L-1) s^(1-L) (s^(2L+1) - x^(2L+1))
            diff = sum(s ** j * xx ** (2 * L - j) for j in range(2 * L + 1)) * v
            EL = -np.sum(GLW / k * WL * xx ** (-L - 1) * s ** (1 - L) * diff, axis=-1) / (2 * L + 1)
            A = np.sum(GLW / k * WL * s ** (1 - L), axis=-1)
            B = np.sum(GLW / k * WL * s ** (L + 2), axis=-1)
            tEL = (L * x ** L * A + (L + 1) * x ** (-L - 1) * B) / (2 * L + 1)
            t2EL = -tEL + L * (L + 1) * EL - x ** 2 * WLx
            E = E + EL[..., None] * PL
            tE = tE + tEL[..., None] * PL
            t2E = t2E + t2EL[..., None] * PL
        return E, tE, t2E

    def Gtop_scaled(self, x):
        """e^{(n-1)x} G^top(x)."""
        ch = self.t.children[self.root]
        parts = [self.esc(w, x) for w in ch]
        Psi = np.ones(x.shape + (len(MU),))
        for E, _, _ in parts:
            Psi = Psi * E
        tPsi = np.zeros_like(Psi)
        t2Psi = np.zeros_like(Psi)
        for a, (E, tE, t2E) in enumerate(parts):
            o = np.ones_like(Psi)
            for b, (Eb, _, _) in enumerate(parts):
                if b != a:
                    o = o * Eb
            tPsi += tE * o
            t2Psi += t2E * o
            for c, (Ec, tEc, _) in enumerate(parts):
                if c > a:
                    o2 = np.ones_like(Psi)
                    for b, (Eb, _, _) in enumerate(parts):
                        if b not in (a, c):
                            o2 = o2 * Eb
                    t2Psi += 2 * tE * tEc * o2

        def proj(F, L):
            return (2 * L + 1) / 2 * np.sum(F * (MW * eval_legendre(L, MU)), axis=-1)
        name = self.dens[self.root]
        if name == "s":
            g = proj(Psi, 0)
        elif name == "pz":
            g = -(1 / 3) / x * (proj(tPsi, 1) + 2 * proj(Psi, 1))
        elif name == "h2":
            g = (1 / 5) / x ** 2 * (proj(t2Psi, 2) + 4 * proj(tPsi, 2) + 3 * proj(Psi, 2))
        elif name == "lap":
            g = 1 / x ** 2 * (proj(t2Psi, 0) + proj(tPsi, 0))
        else:
            raise ValueError(name)
        return x ** 2 * g

    def prediction(self) -> float:
        n = self.n
        poly = [Fraction(1)]
        for v in range(n):
            q = [Fraction(0)] * 3
            s = (n - 1) if v == self.root else -1
            for (a, b), c in M.Pdict(self.dens[v]).items():
                q[a] += c * Fraction(s) ** (a + 2 * b)
            poly = M.pmul(poly, q)
        ang = M.mu_int(poly)
        pk = 1
        for u in range(n):
            if u != self.root:
                pk *= self.k(u) ** 2
        return float((-1) ** (n - 1) * (math.pi / 2) ** ((n - 1) / 2) * float(ang) / pk), float(ang), pk


def top_constant(tt: TailTree, xs=(100.0, 200.0, 400.0, 800.0, 1600.0)):
    a = (5 - tt.n) / 2
    vals = [float(tt.Gtop_scaled(np.array([x]))[0]) * x ** (-a) for x in xs]
    tab, lvl = list(vals), 1
    while len(tab) > 1:
        tab = [(2 ** lvl * tab[i + 1] - tab[i]) / (2 ** lvl - 1) for i in range(len(tab) - 1)]
        lvl += 1
    return tab[0], vals


TREES = {
    "A_pz_s_pz_rootmid": (3, [(0, 1), (0, 2)], 0, {0: "s", 1: "pz", 2: "pz"}),
    "B1_pz_pz_d_rootend": (3, [(0, 1), (1, 2)], 0, {0: "h2", 1: "pz", 2: "pz"}),
    "F_star4_sroot_pz_pz_s": (4, [(0, 1), (0, 2), (0, 3)], 0, {0: "s", 1: "pz", 2: "pz", 3: "s"}),
    "G_star4_pzroot_pz_s_s": (4, [(0, 1), (0, 2), (0, 3)], 0, {0: "pz", 1: "pz", 2: "s", 3: "s"}),
    "J_path4_droot_end": (4, [(0, 1), (1, 2), (2, 3)], 0, {0: "h2", 1: "pz", 2: "pz", 3: "s"}),
    "K_path4_dz2leaf_lap": (4, [(0, 1), (1, 2), (0, 3)], 0, {0: "lap", 1: "s", 2: "dz2", 3: "s"}),
    "L_fork5_droot": (5, [(0, 1), (1, 2), (1, 3), (3, 4)], 0, {0: "h2", 1: "pz", 2: "s", 3: "pz", 4: "s"}),
    "M_star5_pzroot": (5, [(0, 1), (0, 2), (0, 3), (0, 4)], 0, {0: "pz", 1: "pz", 2: "h2", 3: "s", 4: "s"}),
    "N_spider6_dleaf": (6, [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)], 0,
                        {0: "s", 1: "h2", 2: "s", 3: "pz", 4: "pz", 5: "s"}),
}


def branch_fit(tt: TailTree, C: float, eps_list=(1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3)):
    """J(eps) on a log grid in x; least-squares fit of the singular coefficient."""
    n = tt.n
    b = (4 - n) / 2
    xg = np.exp(np.linspace(0.0, np.log(2.0e4), 6001))
    g = np.concatenate([tt.Gtop_scaled(xg[i:i + 500]) for i in range(0, len(xg), 500)])
    Js = []
    for e in eps_list:
        integrand = i0e((n - 1 - e) * xg) * np.exp(-e * xg) * g * xg
        Js.append(float(np.trapezoid(integrand, np.log(xg))))
    eps = np.array(eps_list)
    pred = C * math.gamma(b + 1) / math.sqrt(2 * math.pi * (n - 1))
    if abs(b + 1) < 1e-12:
        basis = [-np.log(eps), np.ones_like(eps), eps * np.log(eps), eps]
    elif abs(b + 1 - round(b + 1)) < 1e-12:
        p = int(round(b + 1))
        basis = [eps ** (-p)] + [eps ** (-p + j) for j in range(1, p)] + [-np.log(eps), np.ones_like(eps), eps]
    else:
        basis = [eps ** (-(b + 1)), eps ** (-(b + 1) + 1), np.ones_like(eps), eps ** (-(b + 1) + 2), eps]
    A = np.vstack(basis).T
    coef = np.linalg.lstsq(A, np.array(Js), rcond=None)[0]
    return {"J": Js, "eps": list(eps_list), "predicted_singular_coef": pred, "fitted": float(coef[0]),
            "rel": abs(coef[0] - pred) / abs(pred)}


if __name__ == "__main__":
    which = [a for a in sys.argv[1:] if a in TREES] or list(TREES)
    fit = [a for a in sys.argv[1:] if a.startswith("fit=")]
    fit = fit[0][4:].split(",") if fit else ["F_star4_sroot_pz_pz_s", "M_star5_pzroot"]
    out = {}
    for name in which:
        n, edges, root, dens = TREES[name]
        tt = TailTree(n, edges, root, dens)
        pred, ang, pk = tt.prediction()
        Cn, vals = top_constant(tt)
        rec = {"n": n, "edges": edges, "densities": dens, "angular_average": ang, "prod_k2": pk,
               "C_predicted": pred, "C_numeric": Cn, "rel": abs(Cn - pred) / abs(pred), "raw": vals}
        log(f"N4 {name}: angular avg {ang}, prod k^2 {pk}, C pred {pred:.12g}, numeric {Cn:.12g}, "
            f"rel {rec['rel']:.1e}")
        if name in fit:
            rec["branch"] = branch_fit(tt, pred)
            log(f"N5 {name}: singular coef predicted {rec['branch']['predicted_singular_coef']:.8g}, fitted "
                f"{rec['branch']['fitted']:.8g}, rel {rec['branch']['rel']:.1e}")
        out[name] = rec
        (ROOT / "debug" / "data" / "tree_attain_l_tail.json").write_text(json.dumps(out, indent=1))
    log("done")
