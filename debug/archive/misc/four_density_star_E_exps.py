"""FD1-FD3 step E: local exponents of L at t = 0, 1, 2, 3, infinity, and its log-free Frobenius
solutions at t = 0 with their minimal annihilator orders (mod p, density_loop_D_gauge.minimal_order).
Indicial polynomial at t0 != 0: apply L to (t - t0)^rho, lowest power of (t - t0) (computed here
by a truncated expansion in s = t - t0, exact)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import density_loop_D6_peel as F6  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402

u, t, rho, s = sp.symbols("u t rho s")
OUT = HERE / "data" / "four_density_star_E.json"


def load_L():
    return [sp.sympify(q).subs(sp.Symbol("u"), u) for q in json.loads((HERE / "data" / "four_density_star_A.json").read_text())["L"]["Q"]]


def theta_coeffs(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    L = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            L[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in L]


def indicial(Lc, t0):
    """theta^k (s^rho) = s^(rho - k) * poly_k(s, rho) with t = t0 + s: theta = (t0 + s) d/ds.
    Represent theta^k s^rho = s^rho * sum_m a_{k,m}(rho) s^(m - k).  Keep only what reaches the
    lowest power of the full sum."""
    n = len(Lc) - 1
    # coefficients c_k(t0 + s) as polynomials in s
    cs = [sp.Poly(sp.expand(c.subs(t, t0 + s)), s) for c in Lc]
    vals = [min((m[0] for m in p.monoms()), default=10**9) if not p.is_zero else 10**9 for p in cs]
    # theta^k s^rho: lowest term t0^k rho (rho-1) ... (rho-k+1) s^(rho-k) plus higher.  Need terms
    # up to a shift M beyond the lowest: compute full expansion of theta^k s^rho as dict power->coef
    lowest = min(v - k for k, v in enumerate(vals) if v < 10**9)
    M = max(v - k for k, v in enumerate(vals) if v < 10**9) - lowest + 1
    cur = {0: sp.Integer(1)}       # s^rho * sum cur[e] s^e
    total = {}
    for k in range(n + 1):
        if not cs[k].is_zero:
            for (pw,), cf in cs[k].terms():
                for e, c in cur.items():
                    if pw + e - lowest < M + 1:
                        total[pw + e] = sp.expand(total.get(pw + e, 0) + cf * c)
        new = {}
        for e, c in cur.items():
            # theta (s^(rho+e)) = (t0 + s) (rho + e) s^(rho+e-1)
            new[e - 1] = sp.expand(new.get(e - 1, 0) + t0 * (rho + e) * c)
            new[e] = sp.expand(new.get(e, 0) + (rho + e) * c)
        cur = {e: c for e, c in new.items() if c != 0 and e - k - 1 <= M + 1}
    lo = min(e for e, c in total.items() if sp.expand(c) != 0)
    return sp.factor(total[lo]), lo


def main() -> None:
    Q = load_L()
    rec = {}
    t0 = time.time()
    Lc = theta_coeffs(Q)
    rec["t=0"] = str(sp.factor(Q[0].subs(u, rho)))
    rec["t=inf_top_Q"] = str(sp.factor(Q[-1].subs(u, rho)))
    print("t=0:", rec["t=0"], flush=True)
    print("t=inf (top Q, t^mu):", rec["t=inf_top_Q"], " deg", sp.Poly(Q[-1], u).degree(), "of order", len(Lc) - 1, flush=True)
    for tt in (1, 2, 3, -1):
        ip, lo = indicial(Lc, sp.Integer(tt))
        rec[f"t={tt}"] = str(ip)
        print(f"t={tt}: {ip}   ({time.time()-t0:.0f}s)", flush=True)
    OUT.write_text(json.dumps(rec, indent=1))
    # log-free Frobenius solutions at t = 0 and their minimal orders (mod p)
    roots = sorted(set(sp.roots(sp.Poly(Q[0], u).as_expr(), u)))
    rec["frob"] = {}
    for lam in roots:
        ser, free = F6.frob_free(Q, lam, NT=400)
        if ser is None:
            rec["frob"][str(lam)] = "log forced"
            print(f"lam {lam}: log forced", flush=True)
            continue
        mo = D.minimal_order(None, ser, 10, 24)
        rec["frob"][str(lam)] = {"free": free, "min_order_deg24": mo}
        print(f"lam {lam}: log-free, free slots {free}, minimal order (deg<=24) {mo}  ({time.time()-t0:.0f}s)", flush=True)
        OUT.write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
