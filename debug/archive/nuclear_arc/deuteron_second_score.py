"""Scored quantities DS1-DS5 (docs/predictions/2026-10-20_deuteron_second.md) plus unscored numbers.

Requires data/deuteron_second_controls.json (all controls PASS) and data/deuteron_second_ode.json.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction as Fr

import mpmath as mp

import deuteron_second_laguerre as L
import deuteron_second_ho as HO

HW_SCAN = [2, 4, 6, 8, 10, 15, 20, 30, 40, 60]
RP, RN2, RDF2 = mp.mpf("0.84075"), mp.mpf("-0.1155"), mp.mpf("0.0331")


def s(x, n=16):
    return mp.nstr(x, n)


def lag_point(N, lam, name):
    pot, hb = L.POTS[name]
    r = L.solve(N, lam, pot, hb, 50, want_vec=True)
    B = -r["E"]
    return {"lam": lam, "B": B, "r2": r["r2"], "r_m": mp.sqrt(r["r2"]) / 2,
            "kappa": mp.sqrt(B / L.mpq(hb)), "c": r["c"]}


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
    return vals, mean, (max(vals) - min(vals)) / mean


def ere(kap, at):
    rho = 2 * (kap - 1 / at) / kap ** 2
    return rho, mp.sqrt(1 / (8 * kap ** 2 * (1 - kap * rho))), mp.sqrt(2 * kap / (1 - kap * rho))


def rel(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    return (a - b) / b


if __name__ == "__main__":
    t0 = time.time()
    ctrl = json.load(open("data/deuteron_second_controls.json"))
    assert ctrl["hydrogen"]["PASS"] and ctrl["regression_MT3"]["PASS"] and ctrl["P2_B"]["PASS"]
    assert all(g["PASS_1e-30"] for g in ctrl["gauss_quad"])
    odeall = json.load(open("data/deuteron_second_ode.json"))
    out = {}
    rs = [mp.mpf(8) + mp.mpf(k) / 10 for k in range(81)]
    pts = {"P1": {}, "P2": {}}
    for name, Ns in (("P1", (20, 30, 40)), ("P2", (10, 15, 20, 25, 30, 35, 40))):
        ode = odeall[name]
        rm_ode, B_ode = mp.mpf(ode["r_m"]), mp.mpf(ode["B"]["richardson"])
        tab = {}
        for N in Ns:
            if name == "P2" and N == 40:
                lam, edge = Fr(ctrl["P2_B"]["lambda_40"]), ctrl["P2_B"]["edge"]
            else:
                pot, hb = L.POTS[name]
                lam, _, _, edge = L.lam_opt(N, pot, hb)
            p = lag_point(N, lam, name)
            pts[name][N] = p
            tab[N] = {"lambda_N": str(lam), "edge_of_coarse_grid": edge, "B": s(p["B"], 20),
                      "r_m": s(p["r_m"], 20), "rel_err_r_m_vs_ode": s(rel(p["r_m"], rm_ode), 6),
                      "rel_err_B_vs_ode": s(rel(p["B"], B_ode), 6)}
            print(name, "N", N, tab[N], round(time.time() - t0), flush=True)
        out[name + "_lambda_N_table"] = tab
        p = pts[name][40]
        pot, hb = L.POTS[name]
        at_k, alpha = L.kohn_scattering_length(40, p["lam"], pot, hb)
        vals, mean, spread = AS_profile(p["c"], p["lam"], p["kappa"], rs)
        AS_ode = mp.mpf(ode["AS_corr_20_30"]["richardson"])
        at_o = mp.mpf(ode["a_t"]["richardson"])
        summ = {"B_basis": s(p["B"], 15), "B_ode": ode["B"]["richardson"],
                "r_m_basis": s(p["r_m"], 15), "r_m_ode": ode["r_m"],
                "a_t_kohn": s(at_k, 15), "a_t_kohn_alpha": s(alpha, 15), "a_t_ode": ode["a_t"]["richardson"],
                "kappa_N40": s(p["kappa"], 15), "kappa_ode": ode["kappa"],
                "A_S_mean_8_16_basis": s(mean, 12), "A_S_spread_8_16_basis": s(spread, 5),
                "A_S_ode_20_30_wallcorr": ode["AS_corr_20_30"]["richardson"],
                "A_S_ode_8_16_spread": ode["AS_816_spread_ode"],
                "A_S_mean_vs_ode": s(rel(mean, AS_ode), 5),
                "A_S_at_8_12_16": [s(vals[0], 12), s(vals[40], 12), s(vals[80], 12)]}
        for route, at in (("ode", at_o), ("kohn", at_k)):
            rho, rere, ASere = ere(p["kappa"], at)
            summ["ERE_" + route] = {"rho": s(rho, 10), "r_ERE": s(rere, 12), "A_S_ERE": s(ASere, 10),
                                    "ratio_(rERE-rm)/rm": s((rere - p["r_m"]) / p["r_m"], 8)}
        rch = mp.sqrt(p["r_m"] ** 2 + RP ** 2 + RN2 + RDF2)
        summ["r_ch"] = s(rch, 10)
        summ["r_ch_vs_CODATA2018_2.12799"] = s(rel(rch, "2.12799"), 5)
        summ["r_ch_vs_CODATA2022_2.12778"] = s(rel(rch, "2.12778"), 5)
        out[name] = summ
        print(name, summ, round(time.time() - t0), flush=True)
    # HO basis for P2, N_r = 40 (DS3) and N_r = 20, 30 (unscored)
    ho = {}
    B2o = mp.mpf(odeall["P2"]["B"]["richardson"])
    rm2o = mp.mpf(odeall["P2"]["r_m"])
    for N in (20, 30, 40):
        row = {}
        for hw in HW_SCAN:
            r = HO.solve(N, hw, "P2")
            B = -r["E"]
            rm = mp.sqrt(r["r2"]) / 2
            Nq = 2 * (N - 1)
            L2 = mp.sqrt(2 * (Nq + mp.mpf(3) / 2 + 2)) * r["b"]
            row[hw] = {"B": s(B, 14), "rel_err_B_vs_ode": s(rel(B, B2o), 5), "r_m": s(rm, 12),
                       "rel_err_r_m": s(rel(rm, rm2o), 5), "b": s(r["b"], 10), "L2": s(L2, 8),
                       "S_err": s(r["S_err"], 3), "above_exact": bool(B <= B2o * (1 + mp.mpf("1e-12")))}
        ho[N] = row
        print("HO", N, {hw: (row[hw]["B"], row[hw]["rel_err_B_vs_ode"]) for hw in HW_SCAN}, round(time.time() - t0), flush=True)
    out["P2_HO"] = ho
    # ---- scoring -------------------------------------------------------------------------
    P1 = out["P1"]
    ds1 = {}
    for q, tgt, gate, vb, vo in (("B", "2.27203", "2e-3", P1["B_basis"], P1["B_ode"]),
                                 ("a_t", "5.4739", "2e-3", P1["a_t_kohn"], P1["a_t_ode"]),
                                 ("r_m", "1.97625", "1e-3", P1["r_m_basis"], P1["r_m_ode"])):
        rb, ro = rel(vb, tgt), rel(vo, tgt)
        ds1[q] = {"target": tgt, "gate": gate, "basis": vb, "ode": vo, "rel_basis": s(rb, 5), "rel_ode": s(ro, 5),
                  "PASS": bool(abs(rb) <= mp.mpf(gate) and abs(ro) <= mp.mpf(gate))}
    ds1["HELD"] = all(ds1[q]["PASS"] for q in ("B", "a_t", "r_m"))
    gap = rel("1.9897134950799908", "1.97625")  # first-run MT-III r_m vs thesis
    ds1["gap_first_run_vs_thesis"] = s(gap, 6)
    ds1["fraction_of_gap_explained"] = s(1 - abs(rel(P1["r_m_ode"], "1.97625")) / abs(gap), 6)
    out["DS1"] = ds1
    p30 = pts["P2"][30]
    e30 = rel(p30["r_m"], rm2o)
    out["DS2"] = {"N_r": 30, "lambda_30": str(p30["lam"]), "r_m_30": s(p30["r_m"], 15), "r_m_ode": odeall["P2"]["r_m"],
                  "rel_err": s(e30, 6), "gate": "1e-6", "HELD": bool(abs(e30) <= mp.mpf("1e-6"))}
    hw_opt = max(HW_SCAN, key=lambda hw: mp.mpf(ho[40][hw]["B"]))
    e3 = mp.mpf(ho[40][hw_opt]["rel_err_B_vs_ode"])
    out["DS3"] = {"N_r": 40, "hw_opt": hw_opt, "B_HO": ho[40][hw_opt]["B"], "B_ode": odeall["P2"]["B"]["richardson"],
                  "rel_err": s(e3, 6), "gate": "1e-4", "HELD": bool(abs(e3) <= mp.mpf("1e-4")),
                  "any_hw_within_1e-4": [hw for hw in HW_SCAN if abs(mp.mpf(ho[40][hw]["rel_err_B_vs_ode"])) <= mp.mpf("1e-4")]}
    P2 = out["P2"]
    sp, mv = mp.mpf(P2["A_S_spread_8_16_basis"]), mp.mpf(P2["A_S_mean_vs_ode"])
    out["DS4"] = {"N_r": 40, "lambda_40": ctrl["P2_B"]["lambda_40"], "spread": s(sp, 5), "mean_vs_ode": s(mv, 5),
                  "A_S_mean": P2["A_S_mean_8_16_basis"], "A_S_ode": P2["A_S_ode_20_30_wallcorr"], "gate": "1e-4",
                  "HELD": bool(abs(sp) <= mp.mpf("1e-4") and abs(mv) <= mp.mpf("1e-4"))}
    rat = mp.mpf(P2["ERE_ode"]["ratio_(rERE-rm)/rm"])
    ratk = mp.mpf(P2["ERE_kohn"]["ratio_(rERE-rm)/rm"])
    out["DS5"] = {"ratio_ode_a_t": s(rat, 8), "ratio_kohn_a_t": s(ratk, 8), "gate": "(0, 0.02005)",
                  "HELD": bool(0 < rat < mp.mpf("0.02005") and 0 < ratk < mp.mpf("0.02005"))}
    out["runtime_s"] = time.time() - t0
    json.dump(out, open("data/deuteron_second_score.json", "w"), indent=1, default=str)
    print(json.dumps({k: out[k] for k in ("DS1", "DS2", "DS3", "DS4", "DS5")}, indent=1, default=str))
