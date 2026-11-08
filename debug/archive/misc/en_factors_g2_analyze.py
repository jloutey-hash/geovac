r"""Gate G2 analysis: analytic e-n elements (debug/data/en_factors_<t>.npz [+ _phase2.npz] from
en_factors_run.py --save [--phase2]) vs the real-space MC (debug/data/en_factors_g2_<t>_s*.npz).
Columns A_k psi, A in [1, F, W1, W2, Fcc, Fc1].  Compared: transition <chi_I|A_k psi>, <chi_I|H A_k psi>
and reference <A_k psi|A_l psi>, <A_k psi|H A_l psi> for every pair touching a NEW column (W1, W2,
Fcc, Fc1).  Control variates (variance reduction only; exactly-known means): <chi_I|psi> - d_I and
<psi|psi> - 1.  Errors from the scatter of per-seed means (independent chains).
Pass: every element within 3 sigma; |rel dev| <= 1e-3 on leading terms (|value| >= 5% of the block max)
-- reported, together with the MC relative error that bounds what 1e-3 can mean.
--both (option b, PI 2026-10-09): analytic debug/data/en_factors_both_<t>.npz (en_factors_run.py
--switch both) vs en_factors_g2b_<t>_s*.npz (en_factors_g2_mc.py --both); columns [1, F, W1, W2, A1, A2],
A_k = the amplitude-switch columns.
Usage: python debug/en_factors_g2_analyze.py 1e-2 [--cv-overlaps] [--both]
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
COLS = ["1", "F", "W1", "W2", "Fcc", "Fc1"]


def cv_est(x, C, seed):
    """control-variate mean; error from the scatter of the per-SEED means (independent Markov chains --
    block groups within one chain are autocorrelated for the core-weighted quantities)."""
    Cc = C - C.mean(0)
    beta = np.linalg.lstsq(Cc, x - x.mean(), rcond=None)[0]
    y = x - C @ beta                          # C has KNOWN mean zero (estimate - exact)
    us = np.unique(seed)
    g = np.array([y[seed == k].mean() for k in us])
    return float(y.mean()), float(g.std(ddof=1) / np.sqrt(len(us)))


def analytic(t, both=False):
    z1 = np.load(os.path.join(DATA, f"en_factors_{'both_' if both else ''}{t:.0e}.npz"), allow_pickle=True)
    p2 = os.path.join(DATA, f"en_factors_{t:.0e}_phase2.npz")
    d = z1['d']
    na = z1['S_I'].shape[1]                                # 2 (dens) or 4 (both: W1, W2, A1, A2)
    if both:
        nan0 = np.zeros((len(d), 0))
        z2 = dict(blk=json.dumps({}), S_I=nan0, H_I=nan0, SPG=[], HPG=[], parts_I=json.dumps({}))
    elif os.path.exists(p2):
        z2 = np.load(p2, allow_pickle=True)
    else:                                                  # phase 1 only: G columns -> NaN (skipped)
        nan2 = np.full((len(d), 2), np.nan)
        z2 = dict(blk=json.dumps({}), S_I=nan2, H_I=nan2, SPG=[np.nan] * 2, HPG=[np.nan] * 2,
                  parts_I=json.dumps({"0": {"Vee": [np.nan] * len(d)}, "1": {"Vee": [np.nan] * len(d)}}))
    N = len(d); S1 = z1['S']; H1 = z1['H']
    H_CI = H1[:N, :N]; t_ = H1[:N, N]; s_ = S1[:N, N]; sig2 = S1[N, N]; gg = H1[N, N]
    gold = json.load(open(os.path.join(os.path.dirname(HERE), "tests", "data", "lih_golden", "golden.json"),
                          encoding="utf-8"))
    items = gold.get("items", gold)
    Fbar = float(items["pi_marriage_exp1.5_1e-3" if N > 20 else "pi_marriage_exp1.5_1e-2"]["values"]
                 ["pieces"]["value"]["Fbar"]["repr"])
    pF = json.loads(str(z1['parts_F'])); pab = json.loads(str(z1['parts_ab']))
    blk = {tuple(k.split("|")): v for k, v in json.loads(str(z2['blk'])).items()}
    E0 = float(d @ H_CI @ d); dt = float(d @ t_)
    HCd = H_CI @ d
    tS = np.column_stack([d, s_ + Fbar * d, z1['S_I'], z2['S_I']])
    tH = np.column_stack([HCd, t_ + Fbar * HCd, z1['H_I'], z2['H_I']])
    K = 6; S = np.full((K, K), np.nan); H = np.full((K, K), np.nan)
    S[0, 0] = 1.0; H[0, 0] = E0
    S[0, 1] = Fbar; H[0, 1] = dt + Fbar * E0
    S[1, 1] = sig2 + Fbar ** 2; H[1, 1] = gg + 2 * Fbar * dt + Fbar ** 2 * E0
    S_P = z1['S_I'].T @ d; HPW = z1['H_I'].T @ d
    for a in range(na):
        S[0, 2 + a] = S_P[a]; H[0, 2 + a] = HPW[a]
        p = pF[str(a)]
        S[1, 2 + a] = p['FW']; H[1, 2 + a] = p['kin'] + p['FWv'] + p['FWee']
    for (a, b) in [(a_, b_) for a_ in range(na) for b_ in range(a_, na)]:
        p = pab[str((a, b))]
        S[2 + a, 2 + b] = p['S']; H[2 + a, 2 + b] = 0.5 * p['grad'] - 0.25 * p['lap'] + p['Vne'] + p['Vee']
    for b in range(len(z2['SPG'])):
        S[0, 4 + b] = z2['SPG'][b]; H[0, 4 + b] = z2['HPG'][b]
    rowof = {'F': 1, 'W1': 2, 'W2': 3, 'Fcc': 4, 'Fc1': 5}
    for (nm, b), v in blk.items():
        k = rowof[nm]; l = 4 + int(b)
        S[k, l] = v['S']; H[k, l] = v['H']
    # V_ee-only parts
    pI = json.loads(str(z2['parts_I']))
    VG = (np.column_stack([np.array(pI[str(b)]['Vee'], float) for b in range(2)]) if not both
          else np.zeros((N, 0)))
    tV = np.column_stack([np.full(N, np.nan), np.full(N, np.nan), z1['Vee_I'], VG])
    Vm = np.full((K, K), np.nan)
    for a in range(na):
        Vm[0, 2 + a] = d @ z1['Vee_I'][:, a]; Vm[1, 2 + a] = pF[str(a)]['FWee']
    for a in range(VG.shape[1]):
        Vm[0, 4 + a] = d @ VG[:, a]
    for (a, b) in [(a_, b_) for a_ in range(na) for b_ in range(a_, na)]:
        Vm[2 + a, 2 + b] = pab[str((a, b))]['Vee']
    for (nm, b), v in blk.items():
        Vm[rowof[nm], 4 + int(b)] = v['Vee']
    # kinetic and V_ne parts (reference side)
    Tm = np.full((K, K), np.nan); Nm = np.full((K, K), np.nan)
    pI2 = {b: pI[str(b)] for b in range(VG.shape[1])}
    for a in range(na):
        Tm[0, 2 + a] = d @ z1['T_I'][:, a]; Nm[0, 2 + a] = d @ z1['Vne_I'][:, a]
        Tm[1, 2 + a] = pF[str(a)]['kin']; Nm[1, 2 + a] = pF[str(a)]['FWv']
    for a in pI2:
        if 'T' in pI2[a]:
            Tm[0, 4 + a] = d @ np.array(pI2[a]['T'], float); Nm[0, 4 + a] = d @ np.array(pI2[a]['Vne'], float)
    for (a, b) in [(a_, b_) for a_ in range(na) for b_ in range(a_, na)]:
        p_ = pab[str((a, b))]
        Tm[2 + a, 2 + b] = 0.5 * p_['grad'] - 0.25 * p_['lap']; Nm[2 + a, 2 + b] = p_['Vne']
    for (nm, b), v in blk.items():
        Tm[rowof[nm], 4 + int(b)] = 0.5 * v['grad'] - 0.25 * (v['lapAB'] + v['lapBA'])
        Nm[rowof[nm], 4 + int(b)] = v['Vne']
    for k in range(K):                                    # symmetric fill
        for l in range(k):
            S[k, l] = S[l, k]; H[k, l] = H[l, k]; Vm[k, l] = Vm[l, k]; Tm[k, l] = Tm[l, k]; Nm[k, l] = Nm[l, k]
    return tS, tH, S, H, d, tV, Vm, Tm, Nm


def main() -> int:
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-2
    both = "--both" in sys.argv
    global COLS
    if both:
        COLS = ["1", "F", "W1", "W2", "A1", "A2"]
    files = sorted(glob.glob(os.path.join(DATA, f"en_factors_g2{'b' if both else ''}_{t:.0e}_s*.npz")))
    zs = [np.load(f) for f in files]
    keys = ["tS", "tH", "S", "H", "tV", "V", "c0", "PP"] + (["Tg", "Vn"] if all('Tg' in z.files for z in zs) else [])
    B = {k: np.concatenate([z[k] for z in zs]) for k in keys}
    seed = np.concatenate([np.full(len(z['c0']), k) for k, z in enumerate(zs)])
    tS, tH, S, H, d, tV, Vm, Tm, Nm = analytic(t, both)
    have_T = all('Tg' in z.files for z in zs)
    N = len(d)
    C = np.column_stack([B['c0'] - d[None, :], B['PP'] - 1.0])
    if "--cv-overlaps" in sys.argv:
        # the transition W-overlaps are certified INDEPENDENTLY (debug/en_factors_onebody_check.py: fine
        # Li-centred quadrature + direct 2x2-minor contraction, agreement 5e-8 rel) -> usable as control
        # variates with known means; they are then excluded from the z test below.
        C = np.column_stack([C, B['tS'][:, :, 2] - tS[None, :, 2], B['tS'][:, :, 3] - tS[None, :, 3]])
    rows = []
    for k in range(2, 6):
        for I in range(N):
            if not ("--cv-overlaps" in sys.argv and k in (2, 3)):
                rows.append(("T", f"<chi_{I}|{COLS[k]} psi>", tS[I, k], B['tS'][:, I, k]))
            rows.append(("T", f"<chi_{I}|H {COLS[k]} psi>", tH[I, k], B['tH'][:, I, k]))
            rows.append(("TV", f"<chi_{I}|Vee {COLS[k]} psi>", tV[I, k], B['tV'][:, I, k]))
    for k in range(6):
        for l in range(6):
            if max(k, l) < 2:
                continue
            rows.append(("P", f"<{COLS[k]} psi|{COLS[l]} psi>", S[k, l], B['S'][:, k, l]))
            rows.append(("P", f"<{COLS[k]} psi|H {COLS[l]} psi>", H[k, l], B['H'][:, k, l]))
            rows.append(("PV", f"<{COLS[k]} psi|Vee {COLS[l]} psi>", Vm[k, l], B['V'][:, k, l]))
            if have_T:
                rows.append(("PT", f"<{COLS[k]} psi|T {COLS[l]} psi>g", Tm[k, l], B['Tg'][:, k, l]))
                rows.append(("PN", f"<{COLS[k]} psi|Vne {COLS[l]} psi>", Nm[k, l], B['Vn'][:, k, l]))
    print(f"files={len(files)} blocks={len(B['c0'])}  N={N}")
    print(f"{'element':28s} {'analytic':>14s} {'MC':>14s} {'err':>9s} {'z':>6s} {'rel':>9s} {'MCrel':>8s}")
    worst_z = 0.0; worst_rel = 0.0; ok_z = True; nlead = 0; nlead_ok = 0
    for blk in ("T", "P", "TV", "PV", "PT", "PN"):
        sub = [r for r in rows if r[0] == blk]
        for kind in (("S", "H") if blk in ("T", "P") else ("V",)):
            sel = [r for r in sub if (kind == "V" or ("|H " in r[1]) == (kind == "H")) and np.isfinite(r[2])]
            if not sel:
                continue
            vmax = max(abs(r[2]) for r in sel)
            for (_b, name, av, x) in sel:
                m, e = cv_est(x, C, seed)
                z = (av - m) / e
                rel = (av - m) / av
                lead = abs(av) >= 0.05 * vmax
                worst_z = max(worst_z, abs(z))
                if abs(z) >= 3:
                    ok_z = False
                if lead:
                    nlead += 1
                    worst_rel = max(worst_rel, abs(rel)); nlead_ok += abs(rel) <= 1e-3
                flag = "  <-- z" if abs(z) >= 3 else ""
                print(f"{name:28s} {av:+14.8f} {m:+14.8f} {e:9.2e} {z:+6.2f} {rel:+9.2e} {e/abs(av):8.1e}"
                      f"{' L' if lead else ''}{flag}")
    print(f"\nworst |z| = {worst_z:.2f} over {len(rows)} elements;  leading terms: {nlead_ok}/{nlead} within "
          f"1e-3 (worst |rel| {worst_rel:.2e})")
    print(f"G2 (3 sigma): {'PASS' if ok_z else 'FAIL'}")
    return 0 if ok_z else 1


if __name__ == "__main__":
    sys.exit(main())
