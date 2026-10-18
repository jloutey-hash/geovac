"""nosplit step 2: complete (bounded-by-proof) search for a splitting of W = Top o Bot.

Split  <=>  there is S = s0 + s1 delta + s2 delta^2, s_k in Q(x), with  Bot o S = 1  mod Top
(right remainder), i.e. S maps Sol(Top) back into Sol(W) as a section of Bot.  With Top
irreducible, S is unique when it exists.

Pole/degree bound (proved in the memo): for a local fundamental system f_i of Top at p and
y_i = S f_i in Sol(W),  (s_k) = Phi^{-1} y  with Phi the delta-Wronskian matrix; by valuations,
   ord_p s_k >= ceil(mu_W(p) - rhomax_Top(p) + k)   at finite p != 0 where Top is singular,
   ord_p s_k >= ceil(mu_W(p))                        at finite p where Top is ordinary,
   ord_p s_k >= ceil(mu_W(p) - rhomax_Top(p))        at p = 0 and p = oo (delta-basis),
where mu_W = min Re exponent of W at p, rhomax_Top = max Re exponent of Top at p.
Then s_k = N_k / D_k with D_k = prod (x-p)^{-bound}, deg N_k <= deg D_k - bound_oo.
Exhaustive exact linear algebra over Q within that bound."""
from __future__ import annotations
import json, sys, time
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import nosplit_core as N

x, r = N.x, N.r


def annihilator(L, G):
    """Monic annihilator of G(Sol L), L of order 3 (requires the image to have dim 3)."""
    n = len(L) - 1
    _, v = N.rdiv(G, L)
    vecs = [v]
    for _ in range(n):
        _, w = N.rdiv(N.compose([0, 1], vecs[-1]), L)
        vecs.append(w)
    M = sp.Matrix([[vecs[k][i] for k in range(n + 1)] for i in range(n)])
    ns = M.nullspace()
    assert len(ns) == 1, len(ns)
    return [sp.factor(sp.cancel(c / ns[0][n])) for c in ns[0]]


def points(*ops):
    fs = set()
    for L in ops:
        for f in N.singular_points(L):
            if sp.Poly(f, x).degree() > 0:
                fs.add(sp.Poly(f, x).monic().as_expr())
    rat, alg = set(), []
    for f in fs:
        pf = sp.Poly(f, x)
        if pf.degree() == 1:
            rat.add(-pf.all_coeffs()[-1])
        else:
            alg.append(("alg", sp.expand(f)))
    return sorted(rat | {sp.Integer(0)}) + alg


def to_d_form(L):
    """delta-form -> d/dx-form coefficient list."""
    out = [sp.Integer(0)] * len(L)
    dk = [sp.Integer(1)]
    for c in L:
        for j, a in enumerate(dk):
            out[j] += c * a
        nxt = [sp.Integer(0)] * (len(dk) + 1)
        for j, a in enumerate(dk):
            nxt[j] += x * sp.diff(a, x)
            nxt[j + 1] += x * a
        dk = [sp.expand(a) for a in nxt]
    return [sp.cancel(c) for c in out]


def exponents_alg(L, q):
    """Exact indicial roots of L at a root alpha of the irreducible polynomial q, computed in
    Q[x]/(q).  Asserts the indicial polynomial has rational coefficients (then all conjugate
    roots share the exponents)."""
    a = to_d_form(L)
    n = len(a) - 1
    a = [sp.cancel(c / a[-1]) for c in a]
    qp = sp.Poly(q, x)
    dq = qp.diff(x).as_expr()

    def at_alpha(R):
        num, den = sp.fraction(sp.cancel(R))
        inv = sp.invert(sp.expand(den), q, x)
        return sp.rem(sp.expand(num * inv), q, x)

    regular, I = True, sp.prod([r - j for j in range(n)])
    for k in range(n):
        num, den = sp.fraction(sp.cancel(a[k]))
        m, dd = 0, sp.Poly(den, x)
        while dd.degree() > 0 and sp.rem(dd, qp).is_zero:
            dd = sp.quo(dd, qp); m += 1
        if m > n - k:
            regular = False
        if m == n - k:
            # (x - alpha)^m a_k = a_k q^m / qtilde^m and qtilde(alpha) = q'(alpha)
            v = sp.rem(sp.expand(at_alpha(sp.cancel(a[k] * q**m)) * sp.invert(sp.expand(dq**m), q, x)), q, x)
            vp = sp.Poly(v, x)
            assert vp.degree() <= 0, ("indicial coefficient not rational", v)
            I += vp.as_expr() * sp.prod([r - j for j in range(k)])
    I = sp.expand(I)
    roots = sp.roots(sp.Poly(I, r))
    assert sum(roots.values()) == n, I
    return I, regular, roots


def ordinary_alg(L, q):
    Lm = [sp.cancel(c / L[-1]) for c in L]
    return all(sp.gcd(sp.Poly(sp.denom(sp.together(c)), x), sp.Poly(q, x)).degree() == 0 for c in Lm)


def ordinary(L, p):
    Lm = [sp.cancel(c / L[-1]) for c in L]
    return all(sp.denom(sp.together(c)).subs(x, p) != 0 for c in Lm)


def bounds(Top, Bot):
    W = N.compose(Top, Bot)
    n = len(Top) - 1
    out = {}
    for p in points(W, Top) + ["oo"]:
        if isinstance(p, tuple):
            _, regW, eW = exponents_alg(W, p[1])
            _, regT, eT = exponents_alg(Top, p[1])
            tord = ordinary_alg(Top, p[1])
        else:
            _, regW, eW = N.exponents(W, p)
            _, regT, eT = N.exponents(Top, p)
            tord = p not in (0, "oo") and ordinary(Top, p)
        assert regW and regT, ("irregular", p)
        mu = min(sp.re(e) for e in eW)
        rho = max(sp.re(e) for e in eT)
        if p == 0 or p == "oo":
            b = [sp.ceiling(mu - rho)] * n
        elif tord:
            b = [sp.ceiling(mu)] * n
        else:
            b = [sp.ceiling(mu - rho + k) for k in range(n)]
        out[str(p)] = {"p": p, "mu_W": mu, "rhomax_Top": rho,
                       "exps_W": sorted(eW, key=lambda z: sp.re(z)),
                       "exps_Top": sorted(eT, key=lambda z: sp.re(z)), "Top_ordinary": tord,
                       "ord_bound": b}
    return W, out


def search(Top, Bot, B, pad=0, rhs=None):
    """Exact search within the bound B (optionally padded by `pad` in every pole order and
    degree). Returns (consistent, solution S or None, nunknowns, ranks, dens, degs)."""
    n = len(Top) - 1
    rhs = rhs or [1] + [0] * (n - 1)
    dens, degs = [], []
    for k in range(n):
        D = sp.Integer(1)
        for key, v in B.items():
            if key == "oo":
                continue
            m = max(0, -int(v["ord_bound"][k])) + pad
            D *= (v["p"][1] if isinstance(v["p"], tuple) else (x - v["p"])) ** m
        D = sp.expand(D)
        dens.append(D)
        degs.append(sp.Poly(D, x).degree() - int(B["oo"]["ord_bound"][k]) + pad)
    cols, labels = [], []
    for k in range(n):
        for m in range(degs[k] + 1):
            S = [sp.Integer(0)] * n
            S[k] = sp.cancel(x**m / dens[k])
            _, rem = N.rdiv(N.compose(Bot, S), Top)
            cols.append(rem); labels.append((k, m))
    target = [sp.Integer(c) for c in rhs]
    den = sp.Integer(1)
    for c in [c for col in cols for c in col] + target:
        den = sp.lcm(den, sp.denom(sp.together(c)))
    rows = []
    for i in range(n):
        polys = [sp.expand(sp.cancel(col[i] * den)) for col in cols]
        tp = sp.expand(sp.cancel(target[i] * den))
        top = max(sp.Poly(p, x).degree() for p in polys + [tp] if p != 0)
        for e in range(top + 1):
            rows.append([sp.Poly(p, x).coeff_monomial(x**e) if p != 0 else 0 for p in polys]
                        + [sp.Poly(tp, x).coeff_monomial(x**e) if tp != 0 else 0])
    Mx = DomainMatrix.from_list_sympy(len(rows), len(cols) + 1, rows).convert_to(sp.QQ)
    A = Mx[:, :len(cols)]
    rA, rAb = A.rank(), Mx.rank()
    sol = None
    if rA == rAb:
        Am = A.to_Matrix(); bm = Mx[:, len(cols):].to_Matrix()
        vsol, params = Am.gauss_jordan_solve(bm)
        vsol = vsol.subs({p: 0 for p in params})
        S = [sp.Integer(0)] * n
        for c, (k, m) in zip(vsol, labels):
            S[k] += c * x**m / dens[k]
        sol = [sp.factor(sp.cancel(s)) for s in S]
    return rA == rAb, sol, len(cols), (rA, rAb), [str(sp.factor(d)) for d in dens], degs


def verify(Top, Bot, S):
    n = len(Top) - 1
    _, rem = N.rdiv(N.compose(Bot, S), Top)
    return all(sp.cancel(a - b) == 0 for a, b in zip(rem, [1] + [0] * (n - 1)))


def ser(B):
    return {k: {kk: (str(vv) if not isinstance(vv, list) else [str(z) for z in vv]) for kk, vv in v.items()}
            for k, v in B.items()}


def run(name, Top, Bot, pads=(0, 2), rec=None):
    t0 = time.time()
    W, B = bounds(Top, Bot)
    print(f"== {name}: ord Top {len(Top)-1}, ord Bot {len(Bot)-1}", flush=True)
    for k, v in B.items():
        print(f"   p={k:6s} exps W {v['exps_W']}  exps Top {v['exps_Top']}  Top ordinary {v['Top_ordinary']}"
              f"  bound {v['ord_bound']}", flush=True)
    res = {"bounds": ser(B)}
    for pad in pads:
        ok, S, nunk, ranks, dens, degs = search(Top, Bot, B, pad=pad)
        good = verify(Top, Bot, S) if ok else None
        print(f"   pad {pad}: unknowns {nunk}, denominators {dens}, numerator degrees {degs}, rank A/[A|b] {ranks},"
              f" split found: {ok}" + (f", verified exactly: {good}, S = {S}" if ok else ""), flush=True)
        res[f"pad{pad}"] = {"unknowns": nunk, "dens": dens, "degs": degs, "ranks": list(ranks), "split": ok,
                            "S": [str(s) for s in S] if S else None, "verified": good}
    res["seconds"] = round(time.time() - t0, 1)
    if rec is not None:
        rec[name] = res
    return res


def main():
    ops = N.chain_ops()
    R3, R1, T2, L3 = ops["R3"], ops["R1"], ops["T2"], ops["L3"]
    E = N.compose(R1, T2)
    rec = {}
    only = sys.argv[1] if len(sys.argv) > 1 else "all"
    if only in ("all", "controls"):
        # positive controls: planted splits Top_c = ann(Bot(Sol L3)), so W_c = lclm(Bot, L3)
        for cname, Bot in [("control: ann(delta Sol L3) o delta", [sp.Integer(0), sp.Integer(1)]),
                           ("control: ann(delta^2 Sol L3) o delta^2", [sp.Integer(0), sp.Integer(0), sp.Integer(1)]),
                           ("control: ann(R1 Sol L3) o R1", R1),
                           ("control: ann(E Sol L3) o E", E)]:
            Top = annihilator(L3, Bot)
            try:
                run(cname, Top, Bot, rec=rec)
            except ValueError as ex:
                print(cname, "skipped:", ex, flush=True)
                rec[cname] = {"skipped": str(ex)}
    if only in ("all", "chain"):
        run("chain partial: R3 o R1", R3, R1, rec=rec)
        run("chain full: R3 o R1 o (theta-2) = W", R3, E, rec=rec)
    (N.HERE / "data" / f"nosplit_search_{only}.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
