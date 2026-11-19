"""tree_attain_l, numerical checks (float64 grid; independent of the model's Gaunt/channel code where stated).

 N1  density channels of the model (d_z, d_z^2, h_2, lap, cancel applied to K0) against closed forms obtained by
     hand-differentiating K0(|r|) in Cartesian coordinates: d_z K0 = -mu K1, d_z^2 K0 = mu^2 K0 + (2mu^2-1)K1/r,
     lap K0 = K0 - K1/r  (Bessel recurrences only), at three radii and five mu.
 N2  the root relocation: direct fibre int r^2 <[P_0(d) K0(t r)] Psi> dr (root derivative on K0(tr), channels of
     Psi by numerical mu-quadrature, no model) against int K0(t r) G(r) dr with G the model element evaluated.
 N3  the transfer phi(a) T(G) = T(a G) for an exact (rational) least-order a of shape B1 (d root) and shape A.
Usage: python debug/tree_attain_l_numeric.py
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_simpson
from scipy.special import k0, k1, eval_legendre
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tree_attain_l_model as M  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
T0 = time.time()


def log(*a) -> None:
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


U = np.linspace(np.log(1e-9), np.log(80.0), 400001)
R = np.exp(U)
MU, MW = leggauss(16)


def cum0(f):
    return cumulative_simpson(f * R, x=U, initial=0.0)          # int_0^r f ds  (ds = r dU)


def cumtop(f):
    """int_r^oo f ds, accumulated from the top of the grid (no cancellation)."""
    return cumulative_simpson((f * R)[::-1], x=-U[::-1], initial=0.0)[::-1]


def closed_density(name: str, r, mu):
    """P(d) K0(|r|) at radius r, cos(angle) mu, by hand (independent of the model)."""
    K0, K1 = k0(r), k1(r)
    dz = -mu * K1
    dz2 = mu ** 2 * K0 + (2 * mu ** 2 - 1) * K1 / r
    lap = K0 - K1 / r
    return {"s": K0, "pz": dz, "dz2": dz2, "lap": lap, "h2": 1.5 * dz2 - 0.5 * lap,
            "cancel": 2.25 * lap - 3.75 * dz2}[name]


def model_density(name: str, r, mu):
    ch = M.density_channels(M.Pdict(name))
    out = 0.0
    for L, f in ch.items():
        for (i, q), c in f.items():
            out = out + float(c) * r ** q * (k0(r) if i == 0 else -r * k1(r)) * eval_legendre(L, mu)
    return out


def N1() -> Dict:
    worst = 0.0
    for name in ("s", "pz", "dz2", "h2", "lap", "cancel"):
        for r in (0.3, 1.7, 6.0):
            for mu in (-0.9, -0.2, 0.0, 0.5, 1.0):
                a, b = closed_density(name, r, mu), model_density(name, r, mu)
                worst = max(worst, abs(a - b) / max(1e-300, abs(a) + 1e-12))
    log(f"N1 density channels: worst relative difference {worst:.2e}")
    return {"worst_rel": worst}


# ------------------------------------------------------------------------------------------------ potentials
def potentials(tree: M.LTree, Lmax: int = 10):
    """Phi_{w,L}(r), theta Phi_{w,L}(r) for every non-root w and L <= Lmax, from W(r, mu) on a mu-grid
    (closed-form densities, numerical projection; does not use tree.W)."""
    phi: Dict = {}
    full: Dict = {}            # w -> Phi_w(r, mu) array (len R, len MU)
    for w in tree.post:
        if w == tree.root:
            continue
        name = DENS[w]
        Wrm = closed_density(name, R[:, None], MU[None, :])
        for c in tree.children[w]:
            Wrm = Wrm * full[c]
        tot = np.zeros_like(Wrm)
        for L in range(Lmax + 1):
            PL = eval_legendre(L, MU)
            WL = (2 * L + 1) / 2 * (Wrm * (MW * PL)[None, :]).sum(axis=1)
            if np.max(np.abs(WL)) < 1e-14 * max(1.0, np.max(np.abs(Wrm))):
                continue
            A = cum0(WL * R ** (L + 2))
            B = cumtop(WL * R ** (1 - L))
            Ph = (A * R ** (-L - 1) + B * R ** L) / (2 * L + 1)
            tPh = (-(L + 1) * A * R ** (-L - 1) + L * B * R ** L) / (2 * L + 1)
            phi[(w, L)] = (Ph, tPh)
            tot = tot + Ph[:, None] * PL[None, :]
        full[w] = tot
    return phi, full


def eval_model(vec: Dict, phi: Dict) -> np.ndarray:
    K0 = k0(R)
    D = -R * k1(R)
    out = np.zeros_like(R)
    for (k, i, kept, e), poly in vec.items():
        f = K0 ** (k - i) * D ** i
        for (y, L, st) in kept:
            f = f * phi[(y, L)][st]
        for m, c in enumerate(poly):
            if c:
                out = out + float(c) * R ** M.Model.rpow_key(e, m) * f
    return out


DENS: Dict[int, str] = {}


def N2(label: str, n: int, edges, dens: Dict[int, str], ts=(0.7, 1.5, 2.6)) -> Dict:
    global DENS
    DENS = dens
    tree = M.LTree(n, edges, 0, {v: M.Pdict(dens[v]) for v in range(n)})
    phi, full = potentials(tree)
    Psi = np.ones((len(R), len(MU)))
    for w in tree.children[0]:
        Psi = Psi * full[w]
    G = eval_model(M.Model(tree).G(), phi)
    deg = max(a + 2 * b for (a, b) in M.Pdict(dens[0]))
    rows = []
    for t in ts:
        rho0 = t ** deg * closed_density(dens[0], t * R[:, None], MU[None, :])    # P_0(d) K0(t|r|)
        direct = np.trapezoid(R ** 3 * 0.5 * (rho0 * Psi * MW[None, :]).sum(axis=1), U)
        reloc = np.trapezoid(k0(t * R) * G * R, U)
        rows.append({"t": t, "direct": direct, "relocated": reloc, "rel": abs(direct - reloc) / abs(direct)})
        log(f"N2 {label} t={t}: direct {direct:.12e}  relocated {reloc:.12e}  rel {rows[-1]['rel']:.1e}")
    return {"label": label, "rows": rows}


# ------------------------------------------------------------------------------------------------ N3 transfer
def exact_a(tree: M.LTree):
    mdl = M.Model(tree, None)
    G = mdl.G()
    vs = [G]
    z = sp.symbols("z")

    def tosym(poly):
        return sum(sp.Rational(c.numerator, c.denominator) * z ** m for m, c in enumerate(poly) if c)
    for N in range(1, 30):
        vs.append(mdl.theta(vs[-1]))
        cols = sorted({k for u in vs for k in u if not mdl.pure(k)}, key=str)
        Mx = sp.Matrix([[tosym(u.get(k, [])) for k in cols] for u in vs])
        if Mx.subs(z, sp.Rational(7, 5)).rank() < len(vs):
            break
    dM = DomainMatrix.from_Matrix(Mx.T).to_field()
    ns = dM.nullspace().to_Matrix()
    cvec = [sp.together(e) for e in ns.row(0)]
    den = sp.lcm([sp.fraction(e)[1] for e in cvec])
    cvec = [sp.expand(sp.cancel(e * den)) for e in cvec]
    cont = sp.gcd_list([e for e in cvec if e != 0])
    cvec = [sp.expand(sp.cancel(e / cont)) for e in cvec]
    a = {}
    for j, cj in enumerate(cvec):
        if cj != 0:
            for (k,), cf in sp.Poly(cj, z).terms():
                a[(k, j)] = Fraction(int(sp.fraction(cf)[0]), int(sp.fraction(cf)[1]))
    tot: Dict = {}
    for (k, j), c in a.items():
        for key, poly in vs[j].items():
            mdl.padd(tot, key, poly, k, c)
    tot = mdl.trim(tot)
    assert all(mdl.pure(k) for k in tot) and tot
    return a, tot, N


_KC: Dict = {}


def kernel(k: int, j: int):
    """(y^-2 theta^2)^k (-1 - theta)^j K0(y) = al(y) K0 + be(y) K1  (from debug/tree_attain_algebra.py)."""
    if (k, j) in _KC:
        return _KC[(k, j)]
    yy = sp.symbols("y")
    al, be = sp.Integer(1), sp.Integer(0)

    def th_(al, be):
        return (sp.expand(yy * sp.diff(al, yy) - yy * be), sp.expand(yy * sp.diff(be, yy) - yy * al - be))
    for _ in range(j):
        a1, b1 = th_(al, be)
        al, be = sp.expand(-al - a1), sp.expand(-be - b1)
    for _ in range(k):
        a1, b1 = th_(*th_(al, be))
        al, be = sp.expand(a1 / yy ** 2), sp.expand(b1 / yy ** 2)
    f = (sp.lambdify(yy, al, "numpy"), sp.lambdify(yy, be, "numpy"))
    _KC[(k, j)] = f
    return f


def N3(label: str, n: int, edges, dens: Dict[int, str], ts=(1.5, 2.5)) -> Dict:
    global DENS
    DENS = dens
    tree = M.LTree(n, edges, 0, {v: M.Pdict(dens[v]) for v in range(n)})
    a, aG, N = exact_a(tree)
    w = max(2 * k + j for (k, j) in a)
    tt = sp.symbols("t")
    lead = sp.factor(sum(sp.Rational(c.numerator, c.denominator) * (-tt) ** j
                         for (k, j), c in a.items() if 2 * k + j == w))
    log(f"N3 {label}: exact least order {N}, weight {w}, symbol at (1,-t): {lead}")
    phi, _ = potentials(tree)
    G = eval_model(M.Model(tree).G(), phi)
    g = eval_model(aG, phi)
    rows = []
    for t in ts:
        Y = t * R
        terms = []
        for (k, j), c in a.items():
            fa, fb = kernel(k, j)
            ker = fa(Y) * k0(Y) + fb(Y) * k1(Y)
            terms.append(float(c) * np.trapezoid(ker * R ** (2 * k) * G * R, U))
        lhs = float(np.sum(terms))
        rhs = float(np.trapezoid(k0(Y) * g * R, U))
        rows.append({"t": t, "phi_a_T_G": lhs, "T_aG": rhs, "rel": abs(lhs - rhs) / abs(rhs),
                     "largest_term": float(np.max(np.abs(terms)))})
        log(f"N3 {label} t={t}: phi(a)T(G) {lhs:.12e}  T(aG) {rhs:.12e}  rel {rows[-1]['rel']:.1e}  "
            f"largest term {rows[-1]['largest_term']:.1e}")
    return {"label": label, "order": N, "weight": w, "symbol": str(lead), "rows": rows}


if __name__ == "__main__":
    out = {"N1": N1()}
    out["N2"] = [] if "skipN2" in sys.argv else [
        N2("A_pz_s_pz_rootmid", 3, [(0, 1), (0, 2)], {0: "s", 1: "pz", 2: "pz"}),
        N2("B1_pz_pz_d_rootend", 3, [(0, 1), (1, 2)], {0: "h2", 1: "pz", 2: "pz"}),
        N2("B2_cancel_rootend", 3, [(0, 1), (1, 2)], {0: "cancel", 1: "pz", 2: "pz"}),
        N2("C_pzroot_mid_pz_s", 3, [(0, 1), (0, 2)], {0: "pz", 1: "pz", 2: "s"}),
        N2("J_path4_droot_end", 4, [(0, 1), (1, 2), (2, 3)], {0: "h2", 1: "pz", 2: "pz", 3: "s"}),
        N2("G_star4_pzroot", 4, [(0, 1), (0, 2), (0, 3)], {0: "pz", 1: "pz", 2: "s", 3: "s"}),
        N2("K_path4_dz2leaf_lap", 4, [(0, 1), (1, 2), (0, 3)], {0: "lap", 1: "s", 2: "dz2", 3: "s"}),
    ]
    out["N3"] = [
        N3("B1_pz_pz_d_rootend", 3, [(0, 1), (1, 2)], {0: "h2", 1: "pz", 2: "pz"}),
        N3("A_pz_s_pz_rootmid", 3, [(0, 1), (0, 2)], {0: "s", 1: "pz", 2: "pz"}),
    ]
    jp = ROOT / "debug" / "data" / "tree_attain_l_numeric.json"
    if "skipN2" in sys.argv and jp.exists():
        old = json.loads(jp.read_text())
        out["N2"] = old.get("N2", [])
    jp.write_text(json.dumps(out, indent=1, default=float))
    log("done")
