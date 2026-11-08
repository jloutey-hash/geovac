"""Scored quantities for DR1, DR2, DR4, DR5 (DR3 in deuteron_radius_dr3.py).

Requires data/deuteron_radius_controls.json (all controls PASS) and data/deuteron_radius_ode.json.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction as Fr

import mpmath as mp

import deuteron_radius_laguerre as L
import deuteron_radius_ho as HO
from deuteron_radius_controls import lam_opt

HB = mp.mpf("41.47")
HW_SCAN = [2, 4, 6, 8, 10, 15, 20, 30, 40, 60]


def s(x, n=16):
    return mp.nstr(x, n)


def lag_point(N, lam):
    r = L.solve(N, lam, L.MT3, L.HB2M, 50, want_vec=True)
    B = -r["E"]
    return {"lam": str(lam), "B": B, "r2": r["r2"], "r_m": mp.sqrt(r["r2"]) / 2,
            "kappa": mp.sqrt(B / HB), "c": r["c"]}


def AS_profile(c, lam, kap, rs):
    mp.mp.dps = 50
    vals = []
    sign = None
    for r in rs:
        u, _ = L.u_of_r(c, lam, r)
        if sign is None:
            sign = 1 if u > 0 else -1
        vals.append(sign * u * mp.e ** (kap * r))
    mean = mp.fsum(vals) / len(vals)
    spread = (max(vals) - min(vals)) / mean
    return vals, mean, spread


if __name__ == "__main__":
    t0 = time.time()
    ctrl = json.load(open("data/deuteron_radius_controls.json"))
    assert ctrl["hydrogen"]["PASS"] and ctrl["B_PASS_basis"] and ctrl["B_PASS_ode"]
    assert ctrl["a_t_PASS_basis"] and ctrl["a_t_PASS_ode"] and ctrl["singlet"]["all_positive"]
    ode = json.load(open("data/deuteron_radius_ode.json"))
    rm_ode = mp.mpf(ode["r_m"])
    AS_ode = mp.mpf(ode["AS_corr_20_30"]["richardson"])
    out = {"ode": {"r_m": ode["r_m"], "B": ode["B"]["richardson"], "a_t": ode["a_t"]["richardson"],
                   "A_S_20_30_wallcorr": ode["AS_corr_20_30"]["richardson"],
                   "A_S_20_30_raw": ode["AS_raw_20_30"]["richardson"],
                   "A_S_8_16_spread_ode": ode["AS_816_spread_ode"]}}
    # 1. lambda_N table
    lamtab = {}
    pts = {}
    for N in (5, 10, 15, 20, 25, 30, 35, 40):
        if N == 40:
            lam = Fr(ctrl["lambda_40"])
        else:
            lam, _, _ = lam_opt(N)
        p = lag_point(N, lam)
        pts[N] = p
        lamtab[N] = {"lambda_N": str(lam), "B": s(p["B"], 20), "r_m": s(p["r_m"], 20),
                     "rel_err_r_m_vs_ode": s((p["r_m"] - rm_ode) / rm_ode, 6),
                     "rel_err_B_vs_ode": s((p["B"] - mp.mpf(ode["B"]["richardson"])) / mp.mpf(ode["B"]["richardson"]), 6)}
        print("lambda_N", N, lamtab[N], round(time.time() - t0))
    out["lambda_N_table"] = lamtab
    # 2. kappa-multiple scan at N = 30, 40
    scan = {}
    for N in (30, 40):
        kap = pts[N]["kappa"]
        row = {}
        for m in (0.5, 1, 2, 3, 4, 6, 8):
            lam = Fr(round(float(m * kap) * 10000), 10000)
            p = lag_point(N, lam)
            row[str(m)] = {"lambda": str(lam), "B": s(p["B"], 15), "r_m": s(p["r_m"], 15),
                           "rel_err_r_m": s((p["r_m"] - rm_ode) / rm_ode, 5)}
        scan[N] = row
        print("scan", N, row, round(time.time() - t0))
    out["kappa_scan"] = scan
    # 3. DR2: lambda = 2 kappa, 3 kappa (kappa from N=40 at lambda_40), N_r = 40
    kap40 = pts[40]["kappa"]
    rs = [mp.mpf(8) + mp.mpf(k) / 10 for k in range(81)]
    dr2 = {}
    for m in (2, 3, "opt"):
        lam = pts[40]["lam"] if m == "opt" else Fr(round(float(m * kap40) * 10000), 10000)
        lam = Fr(lam)
        p = pts[40] if m == "opt" else lag_point(40, lam)
        vals, mean, spread = AS_profile(p["c"], lam, p["kappa"], rs)
        dr2[str(m)] = {"lambda": str(lam), "B_N": s(p["B"], 15), "kappa_N": s(p["kappa"], 15),
                       "r_m": s(p["r_m"], 15), "AS_mean_8_16": s(mean, 12), "AS_spread_rel": s(spread, 5),
                       "AS_mean_vs_ode_rel": s((mean - AS_ode) / AS_ode, 5),
                       "AS_at_8_12_16": [s(vals[0], 12), s(vals[40], 12), s(vals[80], 12)]}
        print("DR2", m, dr2[str(m)], round(time.time() - t0))
    out["DR2_laguerre"] = dr2
    # 4. HO basis
    ho = {}
    for N in (10, 20, 30, 40):
        row = {}
        for hw in HW_SCAN:
            r = HO.solve(N, hw)
            B = -r["E"]
            rm = mp.sqrt(r["r2"]) / 2
            Nq = 2 * (N - 1)
            L2 = mp.sqrt(2 * (Nq + mp.mpf(3) / 2 + 2)) * r["b"]
            row[hw] = {"B": s(B, 12), "r_m": s(rm, 12), "rel_err_r_m": s((rm - rm_ode) / rm_ode, 5),
                       "b": s(r["b"], 10), "L2": s(L2, 10), "S_err": s(r["S_err"], 3),
                       "K_diag_err": s(r["K_diag_err"], 3)}
            if N == 20:
                row[hw]["_c"] = r["c"]
                row[hw]["_b"] = r["b"]
                row[hw]["_E"] = r["E"]
        ho[N] = row
        print("HO", N, {hw: (row[hw]["B"], row[hw]["rel_err_r_m"]) for hw in HW_SCAN}, round(time.time() - t0))
    # HO kinetic off-diagonal check
    chk = HO.build(10, 10)
    off = max(abs(abs(chk["K"][n, n + 1]) - mp.sqrt((n + 1) * (n + mp.mpf(3) / 2))) for n in range(9))
    out["HO_K_offdiag_err"] = s(off, 3)
    # DR2 HO control: N_r = 20, energy-optimal hw in scan
    hw_opt = min(HW_SCAN, key=lambda hw: ho[20][hw]["_E"])
    c, b = ho[20][hw_opt]["_c"], ho[20][hw_opt]["_b"]
    kapHO = mp.sqrt(-ho[20][hw_opt]["_E"] / HB)
    mp.mp.dps = 30
    rates = []
    for r in rs:
        u, du = HO.u_of_r(c, b, r)
        rates.append(-du / u)
    devs = [abs(q - kapHO) / kapHO for q in rates]
    devs_exact = [abs(q - mp.mpf(ode["kappa"])) / mp.mpf(ode["kappa"]) for q in rates]
    out["DR2_HO_control"] = {"N_r": 20, "hw_opt": hw_opt, "B": s(-ho[20][hw_opt]["_E"], 12),
                             "kappa_N_HO": s(kapHO, 10), "max_rel_dev": s(max(devs), 5),
                             "max_rel_dev_vs_exact_kappa": s(max(devs_exact), 5),
                             "r_at_max": s(rs[devs.index(max(devs))], 4),
                             "rate_at_8_12_16": [s(rates[0], 8), s(rates[40], 8), s(rates[80], 8)]}
    print("HO control", out["DR2_HO_control"])
    for N in ho:
        for hw in ho[N]:
            for k in ("_c", "_b", "_E"):
                ho[N][hw].pop(k, None)
    out["HO"] = ho
    # DR1 summary
    lag30 = (pts[30]["r_m"] - rm_ode) / rm_ode
    ho30 = {hw: mp.mpf(ho[30][hw]["rel_err_r_m"]) for hw in HW_SCAN}
    out["DR1"] = {"laguerre_N30_rel_err": s(lag30, 6), "laguerre_PASS": abs(lag30) <= mp.mpf("1e-6"),
                  "HO_N30_min_abs_rel_err": s(min(abs(v) for v in ho30.values()), 5),
                  "HO_PASS": all(abs(v) >= mp.mpf("1e-4") for v in ho30.values())}
    # DR4 and DR5 raw numbers (r_ch in the memo with re-read r_p, r_n^2)
    rm40 = pts[40]["r_m"]
    err40 = (rm40 - rm_ode) / rm_ode
    resid = abs(rm40 - mp.mpf("1.969")) / mp.mpf("1.969")
    out["DR4"] = {"r_m_N40": s(rm40, 15), "basis_err_N40": s(err40, 5), "residual_vs_1.969": s(resid, 6),
                  "ratio_residual_over_basis": s(resid / abs(err40), 5),
                  "r_m_in_window_1.99_2.05": bool(mp.mpf("1.99") <= rm40 <= mp.mpf("2.05"))}
    kap = pts[40]["kappa"]
    for route, at in (("ode", mp.mpf(ode["a_t"]["richardson"])),
                      ("kohn", mp.mpf(ctrl["kohn_a_t"]["40"]["a_t"]))):
        rho = 2 * (kap - 1 / at) / kap ** 2
        rere = mp.sqrt(1 / (8 * kap ** 2 * (1 - kap * rho)))
        AS_ere = mp.sqrt(2 * kap / (1 - kap * rho))
        out[f"DR5_{route}"] = {"kappa_N40": s(kap, 15), "a_t": s(at, 15), "rho": s(rho, 10),
                               "r_ERE": s(rere, 12), "A_S_ERE": s(AS_ere, 10),
                               "ratio": s((rere - rm40) / rm40, 6),
                               "PASS": bool(0 < (rere - rm40) / rm40 < mp.mpf("0.02"))}
    out["runtime_s"] = time.time() - t0
    json.dump(out, open("data/deuteron_radius_score.json", "w"), indent=1, default=str)
    print(json.dumps({k: out[k] for k in ("DR1", "DR4", "DR5_ode", "DR5_kohn", "DR2_laguerre", "DR2_HO_control")},
                     indent=1, default=str))
