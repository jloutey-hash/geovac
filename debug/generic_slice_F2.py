"""GS4/GS5 local exponents, fast route (replaces generic_slice_F.py's symbolic-rho route, stopped for
time on W7): convert sum_k c_k theta^k to sum_k p_k(t) D^k (D = d/dt, polynomial p_k after clearing
denominators); at t0 the indicial polynomial is sum over k with ord_t0(p_k) - k minimal of
p_k^(m_k)(t0)/m_k! * rho(rho-1)...(rho-k+1).  At infinity, the same in s = 1/t.
Checked on R_V against step F's printed exponents."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import generic_slice_C as GC  # noqa: E402
from generic_slice_B import compose, coeff_list, t, u  # noqa: E402

HERE = Path(__file__).parent
rho, s = sp.symbols("rho s")


def to_D(Lc):
    """theta^k = sum_j S(k, j) t^j D^j (Stirling numbers of the second kind)."""
    n = len(Lc) - 1
    p = [sp.Integer(0)] * (n + 1)
    for k, c in enumerate(Lc):
        for j in range(k + 1):
            p[j] += c * sp.functions.combinatorial.numbers.stirling(k, j) * t**j
    den = sp.lcm([sp.fraction(sp.cancel(x))[1] for x in p])
    return [sp.expand(sp.cancel(x * den)) for x in p]


def falling(k):
    return sp.prod([rho - i for i in range(k)])


def indicial_D(p, var, x0):
    best, terms = None, []
    for k, pk in enumerate(p):
        if pk == 0:
            continue
        P = sp.Poly(pk, var)
        m = 0
        val = P.eval(x0)
        while val == 0:
            P = P.diff(var)
            m += 1
            val = P.eval(x0) / sp.factorial(m)
        terms.append((m - k, k, val))
    lo = min(e for e, _, _ in terms)
    return sp.factor(sum(v * falling(k) for e, k, v in terms if e == lo))


def exps(Lc):
    p = to_D(Lc)
    out = {}
    for t0 in [0, 1, 2, 3]:
        ip = indicial_D(p, t, t0)
        out[str(t0)] = {str(k): v for k, v in sp.roots(sp.Poly(ip, rho)).items()} if ip.has(rho) else {}
    # infinity: t = 1/s; operator in s: theta_t = -theta_s
    Ls = [sp.cancel(c.subs(t, 1 / s) * (-1) ** k) for k, c in enumerate(Lc)]
    ps = to_D([c.subs(s, t) for c in Ls])   # reuse t as the variable s
    ip = indicial_D(ps, t, 0)
    out["oo"] = {str(k): v for k, v in sp.roots(sp.Poly(ip, rho)).items()}
    return out


def main():
    Dd = json.loads((HERE / "data" / "generic_slice_D2.json").read_text())
    Cd = json.loads((HERE / "data" / "generic_slice_C.json").read_text())
    R = [sp.sympify(c, locals={"t": t}) for c in Dd["R"]]
    E2 = [sp.sympify(c, locals={"t": t}) for c in Cd["E2"]]
    E3 = compose(E2, [sp.Integer(2), sp.Integer(1)])
    W7 = compose(R, E3)
    rec = {}
    for name, Lc in [("R_V", coeff_list(GC.QV)), ("R", R), ("E3", E3), ("W7", W7)]:
        e = exps(Lc)
        for k, v in e.items():
            rec[f"{name}|{k}"] = v
        print(name, e, flush=True)
    (HERE / "data" / "generic_slice_F.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
