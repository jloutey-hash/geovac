"""MG1-MG2 Legendre map search (docs/predictions/2026-10-23_gamma1_mechanism.md; copy of cf_g2_map.py, outputs renamed): an exact map, rational in t, from the
gauged Legendre operator L_Leg(gamma) (degree -1/2 in z = (t^2 + 1 - gamma^2)/(2t), twisted by t^(1/2)) to
an order-two target operator T.  Generalises debug/c2b_feyn_B_map.py (gamma = 3/2, target C2b) to any
gamma and any target given as a theta-form operator in x = t^2 (correlated_triangle_g32 conventions).

Over Q(t): T and L_Leg have coefficients in Q(x) and are invariant under t -> -t, so any Q(t)-map splits
into an even part (a Q(x)-map from L_Leg) and an odd part (t times a Q(x)-map, i.e. a Q(x)-map from L_Leg
gauged by t = x^(1/2)).  Both twists, x^0 and x^(1/2), are searched; that is exhaustive for the twist.
Maps of order <= 1 suffice (the source has order two).

A = (N0 + N1 d/dx) / E,  E = x^a ((x - x1)(x - x2))^k (x - xa)^m  (x1, x2 = source singular points
(1 - gamma)^2, (1 + gamma)^2; xa = 1 - gamma^2, the source's apparent point).
Identity T(A y) = 0 for y in Sol(source), exactly, as a linear system over Q (EC stage E linear algebra).

usage: python cf_g2_map.py control            (gamma = 3/2: recover the CF2 map into C2b)
       python cf_g2_map.py run <pieces json> <key> <gamma> <kmax> <amax> <mmax> <dNmax>"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import flint
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G  # noqa: E402
from c2b_ident_E_map import fq, int_nullspace, apply_c2b_to_image, gauge  # noqa: E402
from c2b_feyn_A import leg_op, t  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"


def leg_x(g):
    c0, c1, c2 = leg_op(sp.Rational(g))
    e2 = sp.expand(sp.cancel(c2 / t).subs(t, sp.sqrt(x)))
    e1 = sp.expand(c1.subs(t, sp.sqrt(x)))
    e0 = sp.expand(sp.cancel(c0 / t).subs(t, sp.sqrt(x)))
    for e in (e2, e1, e0):
        assert not e.has(sp.sqrt(x)) and sp.Poly(e, x) is not None
    b = [e0, 2 * e2 + 2 * e1, 4 * x * e2]
    gg = sp.gcd_list([sp.expand(v) for v in b])
    return [sp.expand(sp.cancel(v / gg)) for v in b]


def source(g, twist):
    src0 = leg_x(g)
    if not twist:
        return src0
    tw = gauge(src0, sp.Rational(1, 2) / x)
    gg = sp.gcd_list(tw)
    return [sp.expand(sp.cancel(v / gg)) for v in tw]


def target_dform(L):
    """theta-form operator (list of fmpz_poly in x) -> [a0, a1, a2] fmpq_poly, d/dx form"""
    a = G.dform(L)
    return [flint.fmpq_poly([flint.fmpq(int(c)) for c in aj.coeffs()]) for aj in a]


def search(a, src, Eexpr, dN):
    b = [fq(v) for v in src]
    E = fq(Eexpr)
    X = flint.fmpq_poly([0, 1])
    zero = flint.fmpq_poly([])
    cols, labels = [], []
    for j in range(dN + 1):
        cols.append(apply_c2b_to_image(a, b, X**j, zero, E)); labels.append(("N0", j))
    for j in range(dN + 2):
        cols.append(apply_c2b_to_image(a, b, zero, X**j, E)); labels.append(("N1", j))
    deg = max(max(c[0].degree(), c[1].degree()) for c in cols) + 1
    rows = []
    for r in range(deg):
        rows.append([c[0].coeffs()[r] if r < len(c[0].coeffs()) else 0 for c in cols])
        rows.append([c[1].coeffs()[r] if r < len(c[1].coeffs()) else 0 for c in cols])
    NS, nul = int_nullspace(rows, len(cols))
    out = {"E": str(sp.factor(Eexpr)), "dN": dN, "nullity": nul}
    if nul >= 1:
        v = [int(NS[i, 0]) for i in range(len(cols))]
        N0 = sum(c * x**j for (nm, j), c in zip(labels, v) if nm == "N0")
        N1 = sum(c * x**j for (nm, j), c in zip(labels, v) if nm == "N1")
        gg = sp.gcd(sp.Poly(N0, x).content() if N0 != 0 else 0, sp.Poly(N1, x).content() if N1 != 0 else 0)
        if gg not in (0, 1):
            N0, N1 = sp.expand(N0 / gg), sp.expand(N1 / gg)
        n0, n1 = fq(N0), fq(N1)
        U, V = apply_c2b_to_image(a, b, n0, n1, E)
        out.update({"N0": str(N0), "N1": str(N1), "exact": U.is_zero() and V.is_zero(), "order_one": N1 != 0})
        cs = list(n0.coeffs()) if not n0.is_zero() else None
        if cs:
            lead = max(range(len(cs)), key=lambda i: abs(cs[i]))
            cs[lead] = cs[lead] * flint.fmpq(1001, 1000)
            Up, Vp = apply_c2b_to_image(a, b, flint.fmpq_poly(cs), n1, E)
            out["perturbed_1e-3_breaks"] = not (Up.is_zero() and Vp.is_zero())
    return out


def control():
    D = json.load(open(DATA / "correlated_triangle_g32_series_final.json"))
    C2b = G.parse_op(D["factor_C2b"])
    a = target_dform(C2b)
    src = source(sp.Rational(3, 2), True)
    ref = json.load(open(DATA / "c2b_feyn_B_map.json"))
    same_src = [sp.expand(sp.sympify(s)) for s in ref["source"]]
    ratio = sp.cancel(src[2] / same_src[2])
    print("control source equals CF2's source up to a constant:",
          all(sp.expand(p - ratio * q) == 0 for p, q in zip(src, same_src)), flush=True)
    q = (4 * x - 1) * (4 * x - 25)
    runs = []
    for k, m, dN in [(3, 1, 20), (5, 1, 22)]:
        r = search(a, src, q**k * (4 * x + 5) ** m, dN)
        print(f"control k = {k}, m = {m}, dN = {dN}: nullity {r['nullity']}",
              {kk: r[kk] for kk in ("exact", "order_one", "perturbed_1e-3_breaks") if kk in r}, flush=True)
        runs.append(r)
    r = runs[-1]
    if r["nullity"] == 1:
        refrun = [rr for rr in ref["runs"] if rr["nullity"] == 1][0]
        N0r, N1r = sp.sympify(refrun["N0"]), sp.sympify(refrun["N1"])
        N0, N1 = sp.sympify(r["N0"]), sp.sympify(r["N1"])
        c = sp.cancel(N0 / N0r)
        print("control map equals the CF2 map up to scale:", c.is_constant() and sp.expand(N1 - c * N1r) == 0, flush=True)
        r["equals_CF2_map_up_to_scale"] = bool(c.is_constant() and sp.expand(N1 - c * N1r) == 0)
    json.dump({"source": [str(v) for v in src], "runs": runs}, open(DATA / "mg_map_control_g1.5.json", "w"), indent=1)


def run(piece_file, key, g, kmax, amax, mmax, dNmax):
    D = json.load(open(piece_file))
    T = G.parse_op(D[key])
    assert len(T) - 1 == 2, "target must have order two"
    a = target_dform(T)
    g = sp.Rational(g)
    x1, x2, xa = (1 - g) ** 2, (1 + g) ** 2, 1 - g**2
    res = {"target": key, "gamma": str(g), "runs": []}
    found = None
    for twist in (True, False):
        src = source(g, twist)
        print("source (twist x^(1/2)):" if twist else "source (no twist):", [str(sp.factor(v)) for v in src], flush=True)
        for k in range(kmax + 1):
            for m in range(mmax + 1):
                for aa in range(amax + 1):
                    E = x**aa * ((x - x1) * (x - x2)) ** k * (x - xa) ** m
                    r = search(a, src, sp.expand(E), dNmax)
                    r.update({"twist": "x^(1/2)" if twist else "1", "k": k, "m": m, "a": aa})
                    print(f"  twist {r['twist']}, x^{aa} ((x-{x1})(x-{x2}))^{k} (x-({xa}))^{m}, dN {dNmax}: nullity {r['nullity']}",
                          {kk: r[kk] for kk in ("exact", "order_one", "perturbed_1e-3_breaks") if kk in r}, flush=True)
                    res["runs"].append(r)
                    if r["nullity"] >= 1:
                        found = r
                        break
                if found:
                    break
            if found:
                break
        res.setdefault("twists_searched", []).append("x^(1/2)" if twist else "1")
        if found:
            # minimise the numerator degree at the found denominator
            src = source(g, twist)
            E = x**found["a"] * ((x - x1) * (x - x2)) ** found["k"] * (x - xa) ** found["m"]
            dN = found["dN"]
            while dN > 0:
                r2 = search(a, src, sp.expand(E), dN - 1)
                if r2["nullity"] >= 1:
                    dN -= 1
                    r2.update({"twist": found["twist"], "k": found["k"], "m": found["m"], "a": found["a"]})
                    found = r2
                else:
                    break
            print("minimal numerator degree:", found["dN"], "nullity", found["nullity"],
                  {kk: found[kk] for kk in ("exact", "order_one", "perturbed_1e-3_breaks") if kk in found}, flush=True)
            res["map"] = found
            break
    json.dump(res, open(DATA / f"mg_map_{key}_g{str(g).replace('/', '_')}.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    if sys.argv[1] == "control":
        control()
    else:
        run(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]))
