"""Controls for DS1-DS5 (run before any scored quantity).

1. Regression: MT-III B, a_t, r_m (Laguerre N_r = 40 at its energy-optimal lambda, re-searched by
   the same lam_opt used for P1/P2; Kohn a_t; ODE route) vs the first run's stored values, 1e-8.
2. Hydrogen exact limit (first run's control, rerun in the new code).
3. Gaussian elements (P2): closed form (erfc + recurrence) vs quadrature at dps 50 of
   <phi_m| V_P2 |phi_n>, N_r <= 5, lambda in {1/2, 1} (two lambdas fixed before computing); gate 1e-30.
   Plus: erfc route vs parabolic-cylinder route, M_q and the full N_r = 40 V matrix; and the
   mp-transform precision check (S on the mp path vs exact).
4. P2 B vs the published 2.202 MeV within 1e-3 (both routes).
5. HO: Gaussian closed form vs quadrature; MT-III HO regression (N_r = 40, 60 MeV) vs the first run.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction as Fr

import mpmath as mp

import deuteron_second_laguerre as L
import deuteron_second_ho as HO


def s(x, n=20):
    return mp.nstr(x, n)


def hydrogen():
    pot = [("y", Fr(-1), Fr(0))]
    hb = Fr(1, 2)
    out = {}
    M = L.build(1, Fr(1), pot, hb)
    assert L.check_orthogonality(M["S"], Fr(1)) and M["Vg"] is None
    E1 = M["Hx"][0][0] / M["S"][0][0]
    out["N1_lam1_E_exact"] = str(E1)
    M2 = L.build(2, Fr(1, 2), pot, hb)
    assert L.check_orthogonality(M2["S"], Fr(1, 2))
    E = Fr(-1, 8)
    H, S = M2["Hx"], M2["S"]
    det = (H[0][0] - E * S[0][0]) * (H[1][1] - E * S[1][1]) - (H[0][1] - E * S[0][1]) * (H[1][0] - E * S[1][0])
    out["N2_lam_half_det_at_minus_1_8"] = str(det)
    r = L.solve(2, Fr(1, 2), pot, hb, 50, want_vec=False)
    out["N2_lam_half_evals_mp50"] = [s(e, 50) for e in r["evals"]]
    out["PASS"] = (E1 == Fr(-1, 2)) and det == 0
    return out


def gauss_quad(lam: Fr, Nmax: int = 5):
    pot, hb = L.POTS["P2"]
    M = L.build(Nmax, lam, pot, hb)
    mp.mp.dps = 50
    Vg = L.to_orthonormal(M["Vg"], M["S"])
    lamm = L.mpq(lam)
    worst = mp.mpf(0)
    rows = []
    for m in range(Nmax):
        for n in range(m, Nmax):
            f = lambda r: L.phi_n(m, lam, r) * L.phi_n(n, lam, r) * (200 * mp.exp(-mp.mpf("1.487") * r * r)
                                                                   - 178 * mp.exp(-mp.mpf("0.639") * r * r))
            pts = [0] + [mp.mpf(k) / (2 * lamm) for k in range(1, 41)] + [mp.inf]
            q = mp.quad(f, pts)
            d = abs(q - Vg[m, n])
            worst = max(worst, d)
            rows.append((m, n, s(Vg[m, n], 35), s(q, 35), s(d, 4)))
    return {"lambda": str(lam), "max_abs_diff": s(worst, 5), "PASS_1e-30": bool(worst <= mp.mpf("1e-30")),
            "elements": rows}


def route_compare(lam: Fr, N: int = 40):
    pot, hb = L.POTS["P2"]
    mp.mp.dps = 120
    worstM = mp.mpf(0)
    for _, g, a in pot:
        b = L.mpq(a) / (4 * L.mpq(lam) ** 2)
        Me = L.M_erfc(2 * N + 2, b)
        for q in range(2 * N + 3):
            Mp = L.M_pcfd(q, b)
            worstM = max(worstM, abs(Me[q] - Mp) / abs(Mp))
    Ma = L.build(N, lam, pot, hb, route="erfc")
    Mb = L.build(N, lam, pot, hb, route="pcfd")
    mp.mp.dps = 60
    Va = L.to_orthonormal(Ma["Vg"], Ma["S"])
    Vb = L.to_orthonormal(Mb["Vg"], Mb["S"])
    worstV = max(abs(Va[i, j] - Vb[i, j]) for i in range(N) for j in range(N))
    # precision of the mp transform: S through the mp path vs exact
    mp.mp.dps = L.WORK_DPS
    from math import factorial
    Sm = [[mp.mpf(factorial(k + l + 2)) / (2 * L.mpq(lam)) for l in range(N)] for k in range(N)]
    St = L._transform_mp(Ma["C"], Sm)
    worstS = max(abs(St[i][j] - L.mpq(Ma["S"][i][j])) / L.mpq(Ma["S"][i][i]) for i in range(N) for j in range(N))
    # transform cancellation: rebuild at WORK_DPS + 150 and compare the orthonormal Gaussian V
    keep = L.WORK_DPS
    L.WORK_DPS = keep + 150
    Mc = L.build(N, lam, pot, hb, route="erfc")
    L.WORK_DPS = keep
    mp.mp.dps = 60
    Vc = L.to_orthonormal(Mc["Vg"], Mc["S"])
    worstP = max(abs(Va[i, j] - Vc[i, j]) for i in range(N) for j in range(N))
    return {"lambda": str(lam), "N": N, "M_q_rel_diff_erfc_vs_pcfd_max": s(worstM, 5),
            "V_orthonormal_dps200_vs_dps350_max_abs": s(worstP, 5),
            "V_orthonormal_erfc_vs_pcfd_max_abs": s(worstV, 5),
            "S_mp_path_vs_exact_max_rel": s(worstS, 5)}


if __name__ == "__main__":
    t0 = time.time()
    out = {}
    out["hydrogen"] = hydrogen()
    print("hydrogen", out["hydrogen"], flush=True)
    # 1. regression
    prev_s = json.load(open("data/deuteron_radius_score.json"))
    prev_c = json.load(open("data/deuteron_radius_controls.json"))
    prev_o = json.load(open("data/deuteron_radius_ode.json"))
    pot, hb = L.POTS["MT3"]
    lam40, E40, coarse, edge = L.lam_opt(40, pot, hb)
    r = L.solve(40, lam40, pot, hb)
    B = -r["E"]
    rm = mp.sqrt(r["r2"]) / 2
    at, _ = L.kohn_scattering_length(40, lam40, pot, hb)
    ode = json.load(open("data/deuteron_second_ode.json"))["MT3"]
    reg = {"lambda_40": str(lam40), "lambda_40_prev": prev_c["lambda_40"], "edge": edge,
           "B_40": s(B), "B_40_prev": prev_c["B_40"], "B_rel": s((B - mp.mpf(prev_c["B_40"])) / B, 5),
           "r_m_40": s(rm), "r_m_40_prev": prev_s["lambda_N_table"]["40"]["r_m"],
           "r_m_rel": s((rm - mp.mpf(prev_s["lambda_N_table"]["40"]["r_m"])) / rm, 5),
           "a_t_kohn": s(at), "a_t_kohn_prev": prev_c["kohn_a_t"]["40"]["a_t"],
           "a_t_rel": s((at - mp.mpf(prev_c["kohn_a_t"]["40"]["a_t"])) / at, 5),
           "ode_B_rel": (ode["B"]["richardson"] - prev_o["B"]["richardson"]) / prev_o["B"]["richardson"],
           "ode_r_m_rel": (ode["r_m"] - prev_s["ode"]["r_m"]) / prev_s["ode"]["r_m"],
           "ode_a_t_rel": (ode["a_t"]["richardson"] - prev_o["a_t"]["richardson"]) / prev_o["a_t"]["richardson"]}
    tol = 1e-8
    reg["PASS"] = (str(lam40) == prev_c["lambda_40"] and abs(float(reg["B_rel"])) <= tol
                   and abs(float(reg["r_m_rel"])) <= tol and abs(float(reg["a_t_rel"])) <= tol
                   and abs(reg["ode_B_rel"]) <= tol and abs(reg["ode_r_m_rel"]) <= tol
                   and abs(reg["ode_a_t_rel"]) <= tol)
    out["regression_MT3"] = reg
    print("regression", reg, round(time.time() - t0), flush=True)
    # 3. Gaussian elements
    gq = [gauss_quad(Fr(1, 2)), gauss_quad(Fr(1))]
    out["gauss_quad"] = gq
    print("gauss quad", [(g["lambda"], g["max_abs_diff"], g["PASS_1e-30"]) for g in gq], round(time.time() - t0), flush=True)
    rc = [route_compare(Fr(1, 2)), route_compare(Fr(1)), route_compare(Fr(936, 625))]
    out["route_compare"] = rc
    print("route compare", rc, round(time.time() - t0), flush=True)
    # 4. P2 B vs 2.202
    pot2, hb2 = L.POTS["P2"]
    lam40p, E40p, coarse2, edge2 = L.lam_opt(40, pot2, hb2)
    B2 = -E40p
    ode2 = json.load(open("data/deuteron_second_ode.json"))["P2"]
    out["P2_B"] = {"lambda_40": str(lam40p), "edge": edge2, "coarse": coarse2, "B_40": s(B2),
                   "B_40_rel_vs_2.202": s((B2 - mp.mpf("2.202")) / mp.mpf("2.202"), 6),
                   "B_ode": ode2["B"]["richardson"],
                   "B_ode_rel_vs_2.202": (ode2["B"]["richardson"] - 2.202) / 2.202}
    out["P2_B"]["PASS"] = (abs(float(out["P2_B"]["B_40_rel_vs_2.202"])) <= 1e-3
                           and abs(out["P2_B"]["B_ode_rel_vs_2.202"]) <= 1e-3)
    print("P2 B", out["P2_B"], round(time.time() - t0), flush=True)
    # 5. HO checks
    out["HO_gauss_quad"] = HO.gauss_quad_check(6, 20, "P2", [(0, 0), (0, 3), (2, 5), (5, 5)])
    print("HO gauss quad", out["HO_gauss_quad"], flush=True)
    rho = HO.solve(40, 60, "MT3")
    out["HO_MT3_regression_N40_hw60"] = {"B": s(-rho["E"], 12), "B_prev": prev_s["HO"]["40"]["60"]["B"],
                                         "S_err": s(rho["S_err"], 3)}
    print("HO MT3", out["HO_MT3_regression_N40_hw60"], flush=True)
    out["runtime_s"] = time.time() - t0
    json.dump(out, open("data/deuteron_second_controls.json", "w"), indent=1, default=str)
