"""Bridge Test 2 (v2, CALIBRATED): atomic correlation entropy vs modular/log-integer basis.

Fixes v1 flaws:
  (1) dps 60 -> 220 (v1 too low: ~n_basis*log10(maxcoeff) ~ 150 digits needed).
  (2) Adds POSITIVE CONTROLS (numbers KNOWN to be in the span). If PSLQ fails to
      recover these, the machinery is broken and any NULL is meaningless. v1 had
      only a null-control, which proves nothing.

Gate: positive controls MUST be HIT. Then:
  targets HIT  -> correlation entropy IS in this ring (headline lives).
  targets NULL -> not in this basis at this ceiling (Track 5 hardened; bridge
                  stops at correlation entropy).
"""
import json
from mpmath import mp, mpf, pslq, log, pi, gamma, exp, sqrt, catalan, zeta

mp.dps = 220

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
    "pi": pi, "pi2": pi**2, "zeta3": zeta(3),
    "varpi": varpi, "G14": G14, "catalan": catalan,
    "e^-pi": exp(-pi), "e^-2pi": exp(-2*pi),
}
names = list(BASIS); vals = [BASIS[n] for n in names]

pc1 = 3*BASIS["1"] - 2*BASIS["log2"] + 5*BASIS["zeta3"]
pc2 = BASIS["log2"] + BASIS["log3"]                       # log 6
pc3 = 7*BASIS["catalan"] - 4*BASIS["zeta3"] + 2*BASIS["log5"]
nullc = mpf("0.0407")*log(3)/log(2) + exp(-2*pi)

def test(label, x, maxcoeff):
    if x is None:
        return {"label": label, "verdict": "SKIP"}
    rel = pslq([x] + vals, maxcoeff=maxcoeff, maxsteps=10**6)
    if rel is None or rel[0] == 0:
        return {"label": label, "verdict": "NULL"}
    c0 = rel[0]
    recon = -sum(rel[i]*vals[i-1] for i in range(1, len(rel))) / c0
    return {"label": label,
            "verdict": ("HIT" if abs(x-recon) < mpf(10)**(-90) else "WEAK"),
            "max_abs_coeff": max(abs(int(c)) for c in rel),
            "residual": mp.nstr(abs(x-recon), 5),
            "terms": {("x" if i == 0 else names[i-1]): int(rel[i])
                      for i in range(len(rel)) if rel[i] != 0}}

out = {"dps": mp.dps, "basis": names, "runs": {}}
for mc in (10**4, 10**6):
    run = {
        "POS_ctrl_1": test("POS_ctrl_1", pc1, mc),
        "POS_ctrl_2_log6": test("POS_ctrl_2", pc2, mc),
        "POS_ctrl_3": test("POS_ctrl_3", pc3, mc),
        "He_n3": test("He_n3", He, mc),
        "Li+_n3": test("Li+_n3", Li, mc),
        "He_n4": test("He_n4", He4, mc),
        "NULL_ctrl": test("NULL_ctrl", nullc, mc),
    }
    out["runs"][f"maxcoeff_{mc}"] = run
    print(f"\n===== maxcoeff = {mc}, dps = {mp.dps} =====")
    for k, r in run.items():
        extra = ""
        if r.get("max_abs_coeff") is not None:
            extra = f"  coeff<= {r['max_abs_coeff']}  resid={r['residual']}"
        print(f"  {k:>16}: {r['verdict']}{extra}")

json.dump(out, open("debug/data/bridge_test2_v2.json", "w"), indent=1, default=str)
valid = all(out["runs"]["maxcoeff_10000"][p]["verdict"] == "HIT"
            for p in ("POS_ctrl_1", "POS_ctrl_2_log6", "POS_ctrl_3"))
print(f"\nMACHINERY VALID (all positive controls HIT at maxcoeff 1e4): {valid}")
print("WROTE debug/data/bridge_test2_v2.json")
