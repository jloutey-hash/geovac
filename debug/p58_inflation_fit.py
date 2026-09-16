r"""Exact closed forms + R-limit for the Paper 58 permitted-density inflation.

builder reformulated O(n^5): sum over (l_p,l_q,l_r,l_s) of radial multiplicity
r(l)=n_max-l times a closed-form m-count with the L-range intersection. This
reaches large n_max (the O(M^2) dominance counter could not), so the period-2
quasi-polynomial branches can be fit with held-out verification.

Validated against (i) the two known test anchors 214 (n=2), 7600 (n=3) and
(ii) the O(M^2) dominance-counter values at n=4,5 (116552, 1046658).
"""
from __future__ import annotations
import sympy as sp

P = lambda *a: print(*a, flush=True)
n = sp.Symbol('n')

DOMINANCE_CHECK = {2: 214, 3: 7600, 4: 116552, 5: 1046658, 6: 6537894,
                   7: 31462568, 8: 124447624, 9: 422541514, 10: 1269630606}


def c_count(mu, n_max):
    mu = abs(mu)
    return 0 if mu > n_max - 1 else (n_max - mu) * (n_max - mu + 1) // 2


def sum_d_squared(n_max):
    ms = range(-(n_max - 1), n_max)
    return sum(sum(c_count(m, n_max) * c_count(t - m, n_max) for m in ms) ** 2
               for t in range(-2 * (n_max - 1), 2 * (n_max - 1) + 1))


def P_mcount(la, lb, M):
    """#{(ma,mb): |ma|<=la, |mb|<=lb, mb-ma=M}."""
    lo = max(-la, -lb - M)
    hi = min(la, lb - M)
    return max(0, hi - lo + 1)


def _start(lo, top):
    return lo if (top - lo) % 2 == 0 else lo + 1


def builder_percenter(n_max):
    tot = 0
    L = range(n_max)  # l = 0..n_max-1
    for lp in L:
        rp = n_max - lp
        for lq in L:
            rpq = rp * (n_max - lq)
            top1 = lp + lq
            adl1 = abs(lp - lq)
            par1 = top1 % 2
            for lr in L:
                rr = n_max - lr
                for ls in L:
                    rrs = rrs_top = None
                    top2 = lr + ls
                    par2 = top2 % 2
                    if par1 != par2:
                        continue
                    w = rpq * rr * (n_max - ls)
                    adl2 = abs(lr - ls)
                    tmin = min(top1, top2)
                    # sum over transfer M1 with |M1|=mu; both signs, mu=0 once
                    inner = 0
                    for mu in range(0, tmin + 1):
                        lo1 = max(adl1, mu)
                        lo2 = max(adl2, mu)
                        if max(_start(lo1, top1), _start(lo2, top2)) > tmin:
                            continue
                        for M1 in ((0,) if mu == 0 else (mu, -mu)):
                            inner += P_mcount(lp, lq, M1) * P_mcount(lr, ls, -M1)
                    tot += w * inner
    return tot


def fit_quasi(name, ns, vals, deg=11):
    """Fit degree-`deg` polys separately on even-n and odd-n; verify held-out."""
    out = {}
    for par, tag in ((0, "even"), (1, "odd")):
        sub = [(x, v) for x, v in zip(ns, vals) if x % 2 == par]
        k = deg + 1
        if len(sub) < k + 1:
            P(f"  [{name}/{tag}] need >= {k+1} pts, have {len(sub)} -- SKIP")
            continue
        poly = sp.expand(sp.interpolate(sub[:k], n))
        held = sub[k:]
        ok = all(int(poly.subs(n, x)) == v for x, v in held)
        P(f"  [{name}/{tag}] deg-{deg} on {k} pts, {len(held)} held-out OK={ok}"
          f"  lead={sp.LC(sp.Poly(poly, n))}")
        out[tag] = poly
    return out


def main():
    P("=== validate fast O(n^5) builder vs anchors + dominance counter ===")
    okall = True
    for nm, want in DOMINANCE_CHECK.items():
        got = 2 * builder_percenter(nm)
        ok = (got == want)
        okall &= ok
        P(f"  n={nm}: builder={got} want={want} {'OK' if ok else 'MISMATCH'}")
    if not okall:
        raise SystemExit("builder reformulation disagrees -- stop.")

    NMAX = 30
    ns = list(range(2, NMAX + 1))
    P(f"\n=== exact counts n=2..{NMAX} ===")
    oc = [2 * sum_d_squared(x) for x in ns]
    bld = [2 * builder_percenter(x) for x in ns]
    gen = [8 * (oc[i] // 2) - oc[i] + bld[i] for i in range(len(ns))]  # 7*oc/... check below
    gen = [7 * oc[i] + bld[i] for i in range(len(ns))]
    for i, x in enumerate(ns):
        P(f"  n={x:2d}  R={gen[i]/bld[i]:.6f}")

    P("\n=== one_center closed form (single degree-11 poly) ===")
    oc_poly = sp.expand(sp.interpolate(list(zip(ns[:12], oc[:12])), n))
    ok = all(int(oc_poly.subs(n, x)) == v for x, v in zip(ns, oc))
    P(f"  fits all {len(ns)} pts: {ok}")
    P(f"  one_center = {oc_poly}")
    oc_lead = sp.LC(sp.Poly(oc_poly, n))
    P(f"  lead(one_center) = {oc_lead}")

    P("\n=== builder quasi-polynomial (period 2) ===")
    bpoly = fit_quasi("builder", ns, bld, deg=11)

    P("\n=== R limit ===")
    for tag in ("even", "odd"):
        if tag not in bpoly:
            continue
        b_lead = sp.LC(sp.Poly(bpoly[tag], n))
        R_expr = sp.simplify((7 * oc_poly + bpoly[tag]) / bpoly[tag])
        lim = sp.limit(R_expr, n, sp.oo)
        P(f"  [{tag}] lead(builder)={b_lead}  lim R = {lim} = "
          f"{sp.nsimplify(lim)} (float {float(lim):.6f})")
        P(f"        7*lead(oc)/lead(builder)+1 = "
          f"{sp.nsimplify(1 + 7 * oc_lead / b_lead)} = {float(1 + 7*oc_lead/b_lead):.6f}")


if __name__ == "__main__":
    main()
