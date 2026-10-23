"""CF5/CF6 checks at gamma = 2: (a) 1/1000 perturbation of the F_3 map breaks F o A = 0 mod L3;
(b) apparent-point checks (mod p) for the irreducible leading factors of C2b and F;
(c) symbolic j = 64 (z^2+3)^3/(1-z^2)^2 = 16 (x+3)^6/(x (x-1)^2 (x-9)^2) at gamma = 2 (the record's item ii);
(d) the maps are stored into data/cf_g2_pieces.json as strings."""
import json, sys
from fractions import Fraction
from pathlib import Path
import flint
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib2 as H
import cf_g2_lib as C

DATA = C.DATA
P = json.loads((DATA / "cf_g2_pieces.json").read_text())
F = G.parse_op(P["factor_F"])
C2b = G.parse_op(P["factor_C2b"])
out = {}
k3 = json.loads((DATA / "cf_g2_10_k3map.json").read_text())
A = G.parse_op(k3["A_numerators_theta_0_1_2"])
s, b = k3["t_power"], -k3["x_minus_4_power"]
nu, a = (s, 0) if s > 0 else (s % 2, (s % 2 - s) // 2)
den = [int(c) for c in (flint.fmpz_poly([-4, 1]) ** b).coeffs()]
ok = H.verify_map(F, nu, a, den, A)
cs = [int(v) for v in A[0].coeffs()]
j = max(range(len(cs)), key=lambda i: abs(cs[i]))
cs2 = list(cs); cs2[j] = cs2[j] * 1001
A2 = [flint.fmpz_poly(cs2)] + [c * 1000 for c in A[1:]]
bad = H.verify_map(F, nu, a, den, A2)
out["k3map"] = {"identity_exact": ok, "perturbed_1e-3_identity": bad}
print("F_3 map: exact", ok, "; with one coefficient times 1.001:", bad, flush=True)
out["apparent"] = {}
for name, L in (("C2b", C2b), ("F", F)):
    for f, e in C.lead_factors(L):
        if f.degree() >= 2:
            r = C.apparent_check(L, f, nprimes=2)
            out["apparent"][f"{name}: degree {f.degree()} (mult {e})"] = r
            print(name, "factor degree", f.degree(), "mult", e, r, flush=True)
x, t = sp.symbols("x t")
g = 2
z2 = (x + 1 - g**2) ** 2 / (4 * x)
jj = sp.cancel(64 * (z2 + 3) ** 3 / (1 - z2) ** 2)
target = 16 * (x + 3) ** 6 / (x * (x - 1) ** 2 * (x - 9) ** 2)
out["j_gamma2"] = {"j": str(sp.factor(jj)), "equals_record": sp.cancel(jj - target) == 0}
print("j at gamma 2:", sp.factor(jj), "equals record:", out["j_gamma2"]["equals_record"], flush=True)
lm = json.loads((DATA / "cf_g2_map_factor_C2b_g2.json").read_text())["map"]
P["map_CF5_Leg_to_C2b"] = {"form": "A = (N0 + N1 d/dx)/E from L_Leg(gamma=2) in x gauged by x^(1/2) (i.e. t^(1/2) twist on top of L_Leg's t^(-1/2)) into factor_C2b; C2b(A y) = 0 for y in Sol(source), exact over Q",
                           "source": json.loads((DATA / "cf_g2_map_factor_C2b_g2.json").read_text()).get("source", None),
                           "E": lm["E"], "N0": lm["N0"], "N1": lm["N1"], "twist": lm["twist"]}
P["map_CF6_F3_to_F"] = {"form": k3["A"], "A_numerators_theta_0_1_2": k3["A_numerators_theta_0_1_2"]}
(DATA / "cf_g2_pieces.json").write_text(json.dumps(P, indent=1))
(DATA / "cf_g2_11_checks.json").write_text(json.dumps(out, indent=1, default=str))
