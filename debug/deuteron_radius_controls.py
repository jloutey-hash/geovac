"""Controls for DR1-DR5 (run before any scored quantity).

1. Hydrogen exact limit in the same Laguerre code (V = -1/r, a.u.): E = -1/2 exactly at N_r = 1,
   lambda = 1; -1/8 an exact eigenvalue at N_r = 2, lambda = 1/2.
2. lambda_40 (energy-optimal on the grid k/10000) for MT-III; B_40 vs 2.23069 MeV (1e-5 relative).
3. a_t by the basis route (zero-energy Kohn/Schwinger functional in the same Laguerre basis plus
   s(r) = 1 - e^{-lambda r}; DEVIATION: the registration names no basis method for a_t, so this
   method was fixed before computing it) at N_r = 40, lambda_40, vs 5.5132 fm.
4. MT-I singlet: lowest Laguerre eigenvalue > 0 at every N_r <= 40 (at lambda_40 of the triplet,
   and at the kappa-multiple scan lambdas), decreasing toward 0.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction as Fr

import mpmath as mp

import deuteron_radius_laguerre as L


def mpf2s(x, n=20):
    return mp.nstr(x, n)


def hydrogen():
    pot = [(Fr(-1), Fr(0))]
    hb = Fr(1, 2)
    out = {}
    M = L.exact_mats(1, Fr(1), pot, hb)
    assert L.check_orthogonality(M["S"], Fr(1))
    E1 = M["H"][0][0] / M["S"][0][0]
    out["N1_lam1_E_exact"] = str(E1)
    M2 = L.exact_mats(2, Fr(1, 2), pot, hb)
    assert L.check_orthogonality(M2["S"], Fr(1, 2))
    E = Fr(-1, 8)
    H, S = M2["H"], M2["S"]
    det = (H[0][0] - E * S[0][0]) * (H[1][1] - E * S[1][1]) - (H[0][1] - E * S[0][1]) * (H[1][0] - E * S[1][0])
    out["N2_lam_half_det_at_minus_1_8"] = str(det)
    out["N2_lam_half_Hu"] = [[str(x) for x in row] for row in H]
    out["N2_lam_half_Su"] = [[str(x) for x in row] for row in S]
    r = L.solve(2, Fr(1, 2), pot, hb, 50, want_vec=False)
    out["N2_lam_half_evals_mp50"] = [mpf2s(e, 50) for e in r["evals"]]
    out["PASS"] = (E1 == Fr(-1, 2)) and det == 0
    return out


def lam_opt(N: int, pot=L.MT3):
    coarse = [Fr(k, 10) for k in range(2, 31, 2)]
    Es = [L.lowest_E(N, lam, pot) for lam in coarse]
    i = min(range(len(Es)), key=lambda j: Es[j])
    lo = int(coarse[max(i - 1, 0)] * 10000)
    hi = int(coarse[min(i + 1, len(coarse) - 1)] * 10000)
    lam, E = L.lambda_search(N, lo, hi, 10000, pot)
    return lam, E, [(str(c), mpf2s(e, 15)) for c, e in zip(coarse, Es)]


if __name__ == "__main__":
    t0 = time.time()
    out = {}
    out["hydrogen"] = hydrogen()
    print("hydrogen", out["hydrogen"]["PASS"], out["hydrogen"]["N1_lam1_E_exact"],
          out["hydrogen"]["N2_lam_half_det_at_minus_1_8"], out["hydrogen"]["N2_lam_half_evals_mp50"])
    lam40, E40, coarse = lam_opt(40)
    B40 = -E40
    out["lambda_40"] = str(lam40)
    out["coarse_scan_N40"] = coarse
    out["B_40"] = mpf2s(B40, 25)
    out["B_40_rel_vs_2.23069"] = mpf2s((B40 - mp.mpf("2.23069")) / mp.mpf("2.23069"), 6)
    out["B_PASS_basis"] = abs((B40 - mp.mpf("2.23069")) / mp.mpf("2.23069")) <= mp.mpf("1e-5")
    print("lambda_40", lam40, "B_40", out["B_40"], out["B_40_rel_vs_2.23069"], time.time() - t0)
    # a_t by the basis (Kohn) route
    at = {}
    for N in (10, 20, 30, 40):
        a, alpha = L.kohn_scattering_length(N, lam40)
        at[N] = {"a_t": mpf2s(a, 20), "alpha": mpf2s(alpha, 20)}
        print("Kohn", N, at[N])
    out["kohn_a_t"] = at
    a40 = mp.mpf(at[40]["a_t"])
    out["a_t_rel_vs_5.5132"] = mpf2s((a40 - mp.mpf("5.5132")) / mp.mpf("5.5132"), 6)
    out["a_t_PASS_basis"] = abs((a40 - mp.mpf("5.5132")) / mp.mpf("5.5132")) <= mp.mpf("1e-5")
    # ODE route (already run by deuteron_radius_ode.py)
    ode = json.load(open("data/deuteron_radius_ode.json"))
    Bo = ode["B"]["richardson"]
    ao = ode["a_t"]["richardson"]
    out["ode_B"] = Bo
    out["ode_a_t"] = ao
    out["ode_B_rel"] = (Bo - 2.23069) / 2.23069
    out["ode_a_t_rel"] = (ao - 5.5132) / 5.5132
    out["B_PASS_ode"] = abs(out["ode_B_rel"]) <= 1e-5
    out["a_t_PASS_ode"] = abs(out["ode_a_t_rel"]) <= 1e-5
    # MT-I singlet
    sing = {}
    lows = []
    for N in range(1, 41):
        e = L.lowest_E(N, lam40, L.MT1S)
        lows.append(e)
    sing["lambda_40_triplet"] = [mpf2s(e, 12) for e in lows]
    sing["all_positive"] = all(e > 0 for e in lows)
    sing["monotone_nonincreasing"] = all(lows[i + 1] <= lows[i] for i in range(len(lows) - 1))
    kap = mp.sqrt(B40 / mp.mpf("41.47"))
    scan = {}
    for m in (0.5, 1, 2, 3, 4, 6, 8):
        lam = Fr(round(float(m * kap) * 10000), 10000)
        scan[str(lam)] = [mpf2s(L.lowest_E(N, lam, L.MT1S), 12) for N in (10, 20, 30, 40)]
    sing["kappa_multiple_scan_N10_20_30_40"] = scan
    sing["scan_all_positive"] = all(mp.mpf(v) > 0 for vs in scan.values() for v in vs)
    sing["ode_a_s"] = ode["singlet_MT1"]["a_s"]
    sing["ode_zero_energy_nodes"] = ode["singlet_MT1"]["nodes_zero_energy"]
    out["singlet"] = sing
    print("singlet", sing["all_positive"], sing["monotone_nonincreasing"], sing["scan_all_positive"])
    out["runtime_s"] = time.time() - t0
    json.dump(out, open("data/deuteron_radius_controls.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k not in ("singlet", "coarse_scan_N40")}, indent=1, default=str))
