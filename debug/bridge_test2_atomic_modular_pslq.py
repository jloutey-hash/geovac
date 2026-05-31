"""Bridge Test 2: is the Track-5 atomic correlation entropy modular / log-integer?

Re-PSLQ S_full(GS) against a FROZEN basis combining (a) the modular / q-series
ring at tau=i (q=e^{-2pi}, lemniscatic constants) and (b) log-of-small-integers
(motivated by Test 1, where the KMS entropy = log(g1) + q-corrections).

PASS = lands with small coeffs in He AND Li+  -> bridge reaches correlation entropy.
FAIL = null -> atomic correlation entropy is its own class (von Neumann of
       algebraic occupations), the bridge BREAKS at geometric confinement.

Audit discipline: basis frozen here BEFORE inspecting results; test He, Li+,
He_n4, and a null-control of comparable magnitude; report coeff sizes.
"""
import json
from mpmath import mp, mpf, pslq, log, pi, gamma, exp, sqrt, catalan, zeta

mp.dps = 60

d = json.load(open("debug/data/sprint_td_track5.json"))
hp = d["s_full_high_precision"]
He  = mpf(hp["He_n3"]["S_full_GS_str"])
Li  = mpf(hp["Li+_n3"]["S_full_GS_str"])
He4 = mpf(hp["He_n4"]["S_full_GS_str"]) if "S_full_GS_str" in hp.get("He_n4", {}) else None

G14   = gamma(mpf(1)/4)
varpi = G14**2 / (2*sqrt(2*pi))
BASIS = {
    "1": mpf(1),
    "log2": log(2), "log3": log(3), "log5": log(5), "log7": log(7),
    "logpi": log(pi), "logGamma14": log(G14), "log_varpi": log(varpi),
    "pi": pi, "pi2": pi**2, "zeta3": zeta(3), "zeta3_over_pi2": zeta(3)/pi**2,
    "varpi": varpi, "G14": G14, "catalan": catalan,
    "e^-pi": exp(-pi), "e^-2pi": exp(-2*pi),
}
names = list(BASIS); vals = [BASIS[n] for n in names]

def test(label, x):
    if x is None:
        return {"label": label, "skipped": True, "verdict": "SKIP"}
    rel = pslq([x] + vals, maxcoeff=10**8, maxsteps=10**6)
    if rel is None or rel[0] == 0:
        return {"label": label, "verdict": "NULL"}
    c0 = rel[0]
    recon = -sum(rel[i]*vals[i-1] for i in range(1, len(rel))) / c0
    return {"label": label, "verdict": ("HIT" if abs(x-recon) < mpf(10)**(-45) else "WEAK"),
            "max_abs_coeff": max(abs(int(c)) for c in rel),
            "residual": mp.nstr(abs(x-recon), 5),
            "terms": {names[i-1]: int(rel[i]) for i in range(1, len(rel)) if rel[i] != 0}}

null_ctrl = mpf("0.0407")*log(3)/log(2) + exp(-2*pi)  # arbitrary, frozen
out = {"dps": mp.dps, "maxcoeff": 10**8, "basis": names,
       "He_n3": test("He_n3", He), "Li+_n3": test("Li+_n3", Li),
       "He_n4": test("He_n4", He4), "null_control": test("null_control", null_ctrl)}
json.dump(out, open("debug/data/bridge_test2.json", "w"), indent=1, default=str)
for k in ("He_n3", "Li+_n3", "He_n4", "null_control"):
    r = out[k]; extra = ""
    if r.get("max_abs_coeff"):
        extra = f"  coeff<= {r['max_abs_coeff']}  resid={r['residual']}"
    print(f"{k:>14}: {r['verdict']}{extra}")
