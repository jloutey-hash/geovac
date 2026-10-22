"""FD1-FD3 step A (pre-registered: docs/predictions/2026-10-21_four_density_star.md, commit 9a33fbef).

Three-leaf star on the slice (t;1,1,1), fingerprint weights:
  S(t) = N t^3 int_0^inf chi(x)^3 K0(t x) x^-1 dx,  chi = int_0^x s K1(s) ds,  N = 16/pi^4
(derivation in step B's docstring).  Mellin: Phi(s) = N A(-3 - s) E(s + 3), A = M[chi^3],
E = M[K0], E(sig + 2) = sig^2 E(sig)  ->  alpha = -3, beta = 3 in density_loop_B_ode.mellin_t_ode.

Controls (registered): the same code with m = 2, alpha = -2 reproduces the chain's order-seven
operator (tests/test_paper59_density_ring.py QC); with K0^3 (alpha = 1, beta = 0) it reproduces
F_4's operator (debug/data/linking_ceiling_B_ode.json, and eq:four_density_op in x = 1/(4 t^2)).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "tests"))
import density_loop_B_ode as B  # noqa: E402

z, s, u = B.z, B.s, B.u
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_star_A.json"


def rebind(e):
    e = sp.sympify(e)
    m = {q: {"t": t, "u": u}[q.name] for q in e.free_symbols if q.name in ("t", "u")}
    return e.xreplace(m)


def same_up_to_scalar(A_, B_):
    B_ = [rebind(b) for b in B_]
    if len(A_) != len(B_):
        return False
    r = sp.cancel(A_[0] / B_[0])
    return bool(r.is_number) and all(sp.expand(a - r * b) == 0 for a, b in zip(A_, B_))


def x_form(Q):
    """t-form sum_i t^(2i) Q_i(theta_t) of an order-4, 3-term operator -> x = 1/(4 t^2) form after
    the gauge x^(1/2), as tests/test_paper59_four_density_rung.py does."""
    th = sp.Symbol("th")
    Lw = [sp.expand(q.subs(u, 2 * th)) for q in Q]
    K = len(Q) - 1
    Lv = {K - i: sp.expand(c.subs(th, -th)) for i, c in enumerate(Lw)}
    Lx = {k: sp.expand(c.subs(th, th + sp.Rational(1, 2)) * 4**k) for k, c in Lv.items()}
    lead = Lx[0].as_poly(th).LC()
    return {k: sp.factor(sp.expand(Lx[k] / lead)) for k in Lx}, th


def main() -> None:
    rec = {}
    t0 = time.time()
    import test_paper59_density_ring as DR
    # control 1: chain, m = 2
    vchi2 = B.sym_power([0, z - 3, 4], 3, 2)
    QC = B.mellin_t_ode(vchi2, -2, 3, lambda sg: sg**2)
    ok1 = same_up_to_scalar(QC, DR.QC)
    ordC = max(sp.Poly(q, u).degree() for q in QC)
    print(f"control chain: order {ordC}, terms {len(QC)}, equals registered QC up to scalar: {ok1}", flush=True)
    rec["control_chain"] = {"order": ordC, "matches_QC": ok1, "Q": [str(q) for q in QC]}
    # control 2: F_4, K0^3
    vK3 = B.sym_power([z, 0], 2, 3)
    Q4 = B.mellin_t_ode(vK3, 1, 0, lambda sg: sg**2)
    LCd = json.loads((HERE / "data" / "linking_ceiling_B_ode.json").read_text())
    Q4ref = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in LCd["4"]["Q"]]
    ok2 = same_up_to_scalar(Q4, Q4ref)
    Lx, th = x_form(Q4)
    target = {0: th**4, 1: -(2 * th + 1) ** 2 * (10 * th**2 + 10 * th + 3), 2: 9 * (2 * th + 1) ** 2 * (2 * th + 3) ** 2}
    ok2b = all(sp.expand(Lx[k] - target[k]) == 0 for k in target)
    print(f"control F_4: Q = {Q4}; equals linking_ceiling Q4: {ok2}; equals eq:four_density_op in x: {ok2b}", flush=True)
    rec["control_F4"] = {"Q": [str(q) for q in Q4], "matches_linking_ceiling": ok2, "matches_eq_four_density_op": ok2b}
    # the star
    t1 = time.time()
    v3 = B.sym_power([0, z - 3, 4], 3, 3)
    print(f"sym cube: order {len(v3) - 1} in theta_x, z-degree {max(sp.Poly(e, z).degree() for e in v3 if e != 0)}"
          f" ({time.time() - t1:.0f}s)", flush=True)
    rec["symcube_v"] = [str(e) for e in v3]
    QS = B.mellin_t_ode(v3, -3, 3, lambda sg: sg**2)
    ordS = max(sp.Poly(q, u).degree() for q in QS)
    print(f"L: order {ordS}, t-degree {2 * (len(QS) - 1)}  ({time.time() - t0:.0f}s total)", flush=True)
    for i, q in enumerate(QS):
        print(f"  t^{2 * i}:", q, flush=True)
    lead = sum(sp.Poly(q, u).coeff_monomial(u**ordS) * t ** (2 * i) for i, q in enumerate(QS))
    print("  leading coefficient:", sp.factor(lead), flush=True)
    rec["L"] = {"order": ordS, "t_degree": 2 * (len(QS) - 1), "Q": [str(q) for q in QS], "lead": str(sp.factor(lead))}
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
