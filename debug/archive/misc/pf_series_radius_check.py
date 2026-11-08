"""PF1-PF4: radius and slice controls for the Cauchy series (pf_4_series.py).
 (i) two radii (4.5 and 6) of the off-slice kite at Mp = 256, dps = 120 must agree coefficient by coefficient where
     both are accurate (aliasing (rho/9)^Mp and round-off 10^-dps (9/rho)^n; the bound used is 100 times their sum); the ratio c_n/c_{n+1} tends to -9;
 (ii) the slice members (1,1;2,2) and (1,1;3,3) from the same code must reproduce the stored 600-digit series
     debug/data/cf_g2_series_g2_M1600_d600_rho2.json and mg_series_g3_M1600_d600_rho2.json coefficient by coefficient.
usage: python pf_series_radius_check.py"""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp

DATA = Path(__file__).parent / "data"
mp.mp.dps = 130
J = lambda n: json.loads((DATA / n).read_text())
out = {}
A = J("pf_series_rad45_M256_d120_rho4.5.json")
B = J("pf_series_rad6_M256_d120_rho6.json")
rows = []
worst = mp.mpf(0)
for n in range(0, 128):
    a, b = mp.mpf(A["phi_n"][n]), mp.mpf(B["phi_n"][n])
    rel = abs(a / b - 1)
    # expected accuracy: rho = 4.5: 1e-120 * 2^n ; rho = 6: max(1e-120 * 1.5^n, (2/3)^256)
    # round-off 10^-dps (R/rho)^n plus the Cauchy aliasing (rho/R)^Mp, R = 9
    exp45 = mp.mpf(10) ** -120 * mp.mpf(2) ** n + (mp.mpf(1) / 2) ** 256
    exp6 = mp.mpf(10) ** -120 * mp.mpf(1.5) ** n + (mp.mpf(2) / 3) ** 256
    bound = 100 * (exp45 + exp6)
    if n in (0, 1, 2, 5, 10, 20, 40, 60, 80, 100, 120):
        rows.append({"n": n, "rel": mp.nstr(rel, 3), "expected_bound": mp.nstr(bound, 3),
                     "ratio_rho45": mp.nstr(a / mp.mpf(A["phi_n"][n + 1]), 12) if n < 127 else None})
    worst = max(worst, rel / bound)
out["two_radii"] = {"rows": rows, "worst_rel_over_expected_bound": mp.nstr(worst, 3), "agree": bool(worst < 1)}
print("two radii: worst rel / expected bound =", mp.nstr(worst, 3), "; agree:", worst < 1)
for r in rows:
    print("   ", r)
for tag, ref, nmax in (("slice22", "cf_g2_series_g2_M1600_d600_rho2.json", 120), ("slice33", "mg_series_g3_M1600_d600_rho2.json", 120)):
    S = J(f"pf_series_{tag}_M256_d120_rho2.json")
    R = J(ref)
    w = mp.mpf(0)
    for n in range(nmax):
        a, b = mp.mpf(S["phi_n"][n]), mp.mpf(R["phi_n"][n])
        rel = abs(a / b - 1)
        bound = 100 * (mp.mpf(10) ** -120 * mp.mpf(2) ** n + (mp.mpf(1) / 2) ** 256)   # round-off + aliasing of the 120-digit, Mp = 256 run at rho = 2, R = 4
        w = max(w, rel / bound)
    out[tag] = {"reference": ref, "coefficients_compared": nmax, "worst_rel_over_expected_bound": mp.nstr(w, 3), "agree": bool(w < 1),
                "rel_n0": mp.nstr(abs(mp.mpf(S["phi_n"][0]) / mp.mpf(R["phi_n"][0]) - 1), 3),
                "rel_n100": mp.nstr(abs(mp.mpf(S["phi_n"][100]) / mp.mpf(R["phi_n"][100]) - 1), 3)}
    print(tag, "vs", ref, ":", out[tag])
(DATA / "pf_series_radius_check.json").write_text(json.dumps(out, indent=1))
