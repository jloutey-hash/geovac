"""GA3 Legendre-map search (docs/predictions/2026-10-23_geminal_shapes_g2.md): an exact map, rational in t, from
the gauged Legendre operator L_Leg(gamma = 2) (degree -1/2 in z = (t^2 - 3)/(2t), twisted by t^(1/2)) to an
order-two target factor.  Copy of debug/ga1_map.py at gamma = 2 (driver over debug/cf_g2_map.py's functions,
imported, not modified).  Both twists x^0 and x^(1/2) are searched; A = (N0 + N1 d/dx)/E,
E = x^a ((x-1)(x-9))^k (x+3)^m.

Control (the record's positive control): the SAME `run` on the exponential shape's order-two factor C2b at
gamma = 2 (debug/data/cf_g2_pieces.json, key factor_C2b) must find the CF5 map, E = ((x-1)(x-9))^4 (x+3),
numerator degree 20, equal to the stored map_CF5_Leg_to_C2b up to scale.
usage: python ga3_map.py control
       python ga3_map.py run <pieces json> <key> <kmax> <amax> <mmax> <dNmax> <out tag>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G  # noqa: E402
import cf_g2_map as CM  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"
GAMMA = sp.Integer(2)


def run(piece_file, key, kmax, amax, mmax, dNmax, tag):
    D = json.load(open(piece_file))
    T = G.parse_op(D[key])
    assert len(T) - 1 == 2, "target must have order two"
    a = CM.target_dform(T)
    g = GAMMA
    x1, x2, xa = (1 - g) ** 2, (1 + g) ** 2, 1 - g**2
    res = {"target": key, "gamma": str(g), "runs": []}
    found = None
    for twist in (True, False):
        src = CM.source(g, twist)
        print("source (twist x^(1/2)):" if twist else "source (no twist):", [str(sp.factor(v)) for v in src], flush=True)
        for k in range(kmax + 1):
            for m in range(mmax + 1):
                for aa in range(amax + 1):
                    E = x**aa * ((x - x1) * (x - x2)) ** k * (x - xa) ** m
                    r = CM.search(a, src, sp.expand(E), dNmax)
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
            src = CM.source(g, twist)
            E = x**found["a"] * ((x - x1) * (x - x2)) ** found["k"] * (x - xa) ** found["m"]
            dN = found["dN"]
            while dN > 0:
                r2 = CM.search(a, src, sp.expand(E), dN - 1)
                if r2["nullity"] >= 1:
                    dN -= 1
                    r2.update({"twist": found["twist"], "k": found["k"], "m": found["m"], "a": found["a"]})
                    found = r2
                else:
                    break
            print("minimal numerator degree:", found["dN"], "nullity", found["nullity"],
                  {kk: found[kk] for kk in ("exact", "order_one", "perturbed_1e-3_breaks") if kk in found}, flush=True)
            res["map"] = found
            res["source"] = [str(v) for v in src]
            break
    json.dump(res, open(DATA / f"ga3_map_{tag}.json", "w"), indent=1)
    return res


def control():
    """the same `run` on the exponential shape's C2b at gamma = 2 must recover the CF5 map"""
    res = run(DATA / "cf_g2_pieces.json", "factor_C2b", 5, 2, 2, 60, "ctrl_C2b_g2")
    ref = json.load(open(DATA / "cf_g2_pieces.json"))["map_CF5_Leg_to_C2b"]
    m = res.get("map")
    ok = False
    if m:
        N0r, N1r = sp.sympify(ref["N0"]), sp.sympify(ref["N1"])
        N0, N1 = sp.sympify(m["N0"]), sp.sympify(m["N1"])
        c = sp.cancel(N0 / N0r)
        sameE = sp.cancel(sp.sympify(m["E"]) / sp.sympify(ref["E"])).is_constant()
        ok = bool(c.is_constant() and sp.expand(N1 - c * N1r) == 0 and sameE and m["twist"] == ref["twist"])
    print("control: the gamma = 2 run on C2b recovers the CF5 map (E, twist, N0, N1 up to scale):", ok, flush=True)
    res["equals_CF5_map_up_to_scale"] = ok
    json.dump(res, open(DATA / "ga3_map_ctrl_C2b_g2.json", "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "control":
        control()
    else:
        run(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8])
