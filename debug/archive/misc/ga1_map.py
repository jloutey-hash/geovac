"""GA1 Legendre-map search: an exact map, rational in t, from the gauged Legendre operator L_Leg(gamma)
(degree -1/2 in z = (t^2 + 1 - gamma^2)/(2t), twisted by t^(1/2)) to an order-two target factor, at
gamma = 3/2.  Driver over debug/cf_g2_map.py's functions (imported, not modified; its `run` writes cf_g2_*
files, so the loop is repeated here with ga1_* outputs).  Both twists x^0 and x^(1/2) are searched
(exhaustive for the twist over Q(t)); A = (N0 + N1 d/dx)/E, E = x^a ((x-1/4)(x-25/4))^k (x+5/4)^m.
usage: python ga1_map.py control                                   (gamma = 3/2: recover the CF2 map into C2b)
       python ga1_map.py run <pieces json> <key> <kmax> <amax> <mmax> <dNmax> <out tag>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G  # noqa: E402
import cf_g2_map as CM  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"


def run(piece_file, key, kmax, amax, mmax, dNmax, tag):
    D = json.load(open(piece_file))
    T = G.parse_op(D[key])
    assert len(T) - 1 == 2, "target must have order two"
    a = CM.target_dform(T)
    g = sp.Rational(3, 2)
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
    json.dump(res, open(DATA / f"ga1_map_{tag}.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    if sys.argv[1] == "control":
        CM.control()
    else:
        run(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8])
