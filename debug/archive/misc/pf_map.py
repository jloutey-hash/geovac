"""PF1 Legendre-map search (docs/predictions/2026-10-23_price_formula_slots.md): an exact map, rational in t, from the
gauged Legendre operator of degree -1/2 in z = (t^2 + B^2 - gamma^2)/(2 B t) (the one-function L(t, B; i gamma) =
pi^2 P_(-1/2)(z)/(4 sqrt(tB)), so H = t^(-1/2) w(z)), optionally twisted by a further x^(1/2), into an order-two target
factor.  Generalises debug/c2b_feyn_A.leg_op (B = 1) and drives debug/cf_g2_map.py's search (imported, not modified).

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).

Sources (record): crossed z_1: (B, gamma) = (2, 3/2), z_1 = (4t^2+7)/(16t); crossed z_2: (1, 3), z_2 = (t^2-8)/(2t);
direct d_1: (1, 3/2), z = (t^2-5/4)/(2t); direct d_2: (2, 3), z = (t^2-5)/(4t).
A = (N0 + N1 d/dx)/E, E = x^a ((x - x1)(x - x2))^k (x - xa)^m, x1 = (B-gamma)^2, x2 = (B+gamma)^2 (the source's
singular points), xa = B^2 - gamma^2 (the source's apparent point, dz/dt = 0).  Both twists.  Identity
T(A y) = 0 for y in Sol(source), exactly over Q (cf_g2_map.search); 1/1000 perturbation control inside.

Control: the same `run` on the stored gamma = 3 slice factor (tests/data/paper59_mg_2026_10_23.json, gammas/3/C2b_theta)
with source (B, gamma) = (1, 3) must recover the stored map (E = ((x-4)(x-16))^4 (x+8), N0, N1 up to scale).
usage: python pf_map.py control
       python pf_map.py run <pieces json> <key> <kmax> <amax> <mmax> <dNmax> <out tag> [sources: z1,z2,d1,d2]"""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G  # noqa: E402
import cf_g2_map as CM  # noqa: E402
from c2b_ident_E_map import gauge  # noqa: E402

x, t = sp.symbols("x t")
DATA = Path(__file__).parent / "data"
SOURCES = {"z1": (sp.Integer(2), sp.Rational(3, 2)), "z2": (sp.Integer(1), sp.Integer(3)),
           "d1": (sp.Integer(1), sp.Rational(3, 2)), "d2": (sp.Integer(2), sp.Integer(3))}


def leg_op_general(B, g, twist=True):
    """[c0, c1, c2] polynomials in t with c2 H'' + c1 H' + c0 H = 0 for H = t^(-1/2) w(z), z = (t^2 + B^2 - g^2)/(2 B t),
    w a Legendre function of degree -1/2: (1 - z^2) w'' - 2 z w' - (1/4) w = 0."""
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    z = (t**2 + B**2 - g**2) / (2 * B * t)
    rt = sp.sqrt(t)
    if twist:
        W, W1, W2 = rt * h0, h0 / (2 * rt) + rt * h1, -h0 / (4 * t * rt) + h1 / rt + rt * h2
    else:
        W, W1, W2 = h0, h1, h2
    zp, zpp = sp.diff(z, t), sp.diff(z, t, 2)
    w1 = W1 / zp
    w2 = (W2 - zpp * w1) / zp**2
    e = sp.together(sp.expand(((1 - z**2) * w2 - 2 * z * w1 - sp.Rational(1, 4) * W) / (rt if twist else 1)))
    num = sp.expand(sp.numer(e))
    cs = [num.coeff(h) for h in (h0, h1, h2)]
    g0 = sp.gcd_list([c for c in cs if c != 0])
    return [sp.expand(sp.cancel(c / g0)) for c in cs]


def leg_x_general(B, g):
    c0, c1, c2 = leg_op_general(B, g)
    e2 = sp.expand(sp.cancel(c2 / t).subs(t, sp.sqrt(x)))
    e1 = sp.expand(c1.subs(t, sp.sqrt(x)))
    e0 = sp.expand(sp.cancel(c0 / t).subs(t, sp.sqrt(x)))
    for e in (e2, e1, e0):
        assert not e.has(sp.sqrt(x)) and sp.Poly(e, x) is not None
    b = [e0, 2 * e2 + 2 * e1, 4 * x * e2]
    gg = sp.gcd_list([sp.expand(v) for v in b])
    return [sp.expand(sp.cancel(v / gg)) for v in b]


def source_general(B, g, twist):
    src0 = leg_x_general(B, g)
    if not twist:
        return src0
    tw = gauge(src0, sp.Rational(1, 2) / x)
    gg = sp.gcd_list(tw)
    return [sp.expand(sp.cancel(v / gg)) for v in tw]


def run(piece_file, key, kmax, amax, mmax, dNmax, tag, sources=("z1", "z2", "d1", "d2"), T=None):
    if T is None:
        D = json.load(open(piece_file))
        T = G.parse_op(D[key])
    assert len(T) - 1 == 2, "target must have order two"
    a = CM.target_dform(T)
    res = {"target": key, "pieces": str(piece_file), "runs": [], "hits": {}}
    for sname in sources:
        B, g = SOURCES[sname]
        x1, x2, xa = (B - g) ** 2, (B + g) ** 2, B**2 - g**2
        found = None
        for twist in (True, False):
            src = source_general(B, g, twist)
            for k in range(kmax + 1):
                for m in range(mmax + 1):
                    for aa in range(amax + 1):
                        E = x**aa * ((x - x1) * (x - x2)) ** k * (x - xa) ** m
                        r = CM.search(a, src, sp.expand(E), dNmax)
                        r.update({"source": sname, "B": str(B), "gamma": str(g), "twist": "x^(1/2)" if twist else "1", "k": k, "m": m, "a": aa})
                        res["runs"].append({kk: r[kk] for kk in ("source", "twist", "k", "m", "a", "dN", "nullity")})
                        if r["nullity"] >= 1:
                            found = r
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                dN = found["dN"]
                while dN > 0:
                    r2 = CM.search(a, src, sp.expand(E), dN - 1)
                    if r2["nullity"] >= 1:
                        dN -= 1
                        r2.update({"source": sname, "B": str(B), "gamma": str(g), "twist": found["twist"], "k": found["k"], "m": found["m"], "a": found["a"]})
                        found = r2
                    else:
                        break
                found["source_coefficients_d_dx_form"] = [str(v) for v in src]
                res["hits"][sname] = found
                print(f"  source {sname} (B={B}, gamma={g}): HIT twist {found['twist']}, E = {found['E']}, dN {found['dN']}, nullity {found['nullity']}, "
                      f"exact {found.get('exact')}, order_one {found.get('order_one')}, 1/1000 breaks {found.get('perturbed_1e-3_breaks')}", flush=True)
                break
        if not found:
            print(f"  source {sname} (B={B}, gamma={g}): no map in the box (k<={kmax}, m<={mmax}, a<={amax}, dN<={dNmax}, both twists)", flush=True)
    res["box"] = {"kmax": kmax, "amax": amax, "mmax": mmax, "dNmax": dNmax, "twists": ["x^(1/2)", "1"], "sources": list(sources)}
    json.dump(res, open(DATA / f"pf_map_{tag}.json", "w"), indent=1)
    return res


def control():
    """the stored gamma = 3 slice map must be recovered from source (B, gamma) = (1, 3) by the same code path"""
    ref = json.load(open(Path(__file__).parent.parent / "tests" / "data" / "paper59_mg_2026_10_23.json"))["gammas"]["3"]
    T = G.parse_op(ref["C2b_theta"])
    res = run(None, "mg gamma=3 C2b_theta", 5, 2, 2, 30, "ctrl_mg_g3", sources=("z2",), T=T)
    m = res["hits"].get("z2")
    ok = False
    if m:
        rm = ref["map"]
        N0r, N1r = sp.sympify(rm["N0"]), sp.sympify(rm["N1"])
        N0, N1 = sp.sympify(m["N0"]), sp.sympify(m["N1"])
        c = sp.cancel(N0 / N0r)
        sameE = sp.cancel(sp.sympify(m["E"]) / sp.sympify(rm["E"])).is_constant()
        ok = bool(c.is_constant() and sp.expand(N1 - c * N1r) == 0 and sameE and m["twist"] == rm.get("twist", "x^(1/2)"))
    print("control: the stored gamma = 3 map is recovered from the general (B, gamma) = (1, 3) source (E, twist, N0, N1 up to scale):", ok, flush=True)
    res["equals_stored_mg_g3_map_up_to_scale"] = ok
    # the general source at B = 1 equals c2b_feyn_A's leg_op up to a constant
    from c2b_feyn_A import leg_op
    old = leg_op(sp.Integer(3))
    new = leg_op_general(sp.Integer(1), sp.Integer(3))
    ratio = sp.cancel(new[2] / old[2])
    same = all(sp.expand(n - ratio * o) == 0 for n, o in zip(new, old))
    print("control: leg_op_general(1, 3) equals c2b_feyn_A.leg_op(3) up to a constant:", same, flush=True)
    res["general_source_equals_c2b_feyn_A_at_B1"] = bool(same)
    json.dump(res, open(DATA / "pf_map_ctrl_mg_g3.json", "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "control":
        control()
    else:
        srcs = tuple(sys.argv[9].split(",")) if len(sys.argv) > 9 else ("z1", "z2", "d1", "d2")
        run(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8], srcs)
