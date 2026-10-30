"""Scorer for docs/predictions/2026-10-29_ccm_yukawa_coefficient_window.md (YC, commit 0295f53a).

Order: YC-0 (a), (b), (c) first; stop (all items void) if any fails; then the
eight scored points tau in {1, 4} x n_max in {2, 3, 4, 5}, t = tau/(n_max+1/2)^2.
Numerators: continuum_ref window functions (gauge_eps2, kinetic, T = whole window)
on continuum_ref.Sphere with window |lambda| <= n_max + 1/2.  Denominators: (R1)
of continuum_ref (recomputed here at full precision and checked against the
logged 14-decimal values).
Outputs: debug/graded_ac/yc_score.log, debug/graded_ac/yc_score.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from continuum_ref import (build, commutator_D, gauge_eps2, kinetic,  # noqa: E402
                           tail_bound_gauge, tail_bound_kinetic)

LOG: List[str] = []


def log(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


# logged reference values (continuum_ref.log), used by YC-0(a) and as the
# denominator cross-check
LOGGED_FIXED_T = {  # n_max: (kA/kA_inf, kPhi/kPhi_inf, c)
    10: (0.985460145757, 0.999851197287, 3.942427226896),
    12: (0.999671955860, 0.999997716120, 3.998696955984),
    14: (0.999996866872, 0.999999984059, 3.999987531251),
}
LOGGED_DEN = {  # (tau, n_max): (kA_inf, kPhi_inf)
    (1.0, 2): (0.43583439020830, 0.18528852326763),
    (1.0, 3): (0.32401620835791, 0.14922412172569),
    (1.0, 4): (0.25615525714397, 0.12188142097549),
    (1.0, 5): (0.21131126973803, 0.10221061698675),
    (4.0, 2): (0.67297925634500, 0.17085394539082),
    (4.0, 3): (0.57078847739415, 0.20395388070461),
    (4.0, 4): (0.47497707474282, 0.19426247267132),
    (4.0, 5): (0.40189469658341, 0.17579458134149),
}
LOGGED_REF_T01 = (0.35529107367446, 0.16060377028857)


def main() -> None:
    log("yc_score.py -- scorer for YC (docs/predictions/2026-10-29_ccm_yukawa_coefficient_window.md, 0295f53a)")
    J_amb2 = 54
    S, a, b, A, Apg, phi, C = build(J_amb2)
    N_ref = J_amb2 - 3
    T = S.level <= N_ref
    NORM, DELTA = 6.0, 3

    def reference(t: float) -> Dict[str, float]:
        f, s2 = gauge_eps2(S.lam, A, t, T)
        fp, sp2 = gauge_eps2(S.lam, Apg, t, T)
        return dict(kA=f + s2, kP=kinetic(S.lam, C, t, T), kpg=fp + sp2,
                    bA=tail_bound_gauge(N_ref, t, NORM, DELTA), bP=tail_bound_kinetic(N_ref, t, NORM, DELTA))

    def window(n: int, t: float) -> Dict[str, float]:
        idx = np.nonzero(np.abs(S.lam) <= n + 0.5)[0]
        aw = a[idx][:, idx].tocsr()
        bw = b[idx][:, idx].tocsr()
        lw = S.lam[idx]
        Db = commutator_D(lw, bw)
        om = (aw @ Db).tocsr()
        Aw = (om + om.getH()).tocsr()
        Cw = commutator_D(lw, (aw @ bw).tocsr())
        Pw = (Db + Db.getH()).tocsr()
        Tw = np.ones(len(idx), bool)
        f, s2 = gauge_eps2(lw, Aw, t, Tw)
        fp, sp2 = gauge_eps2(lw, Pw, t, Tw)
        return dict(dim=len(idx), kA=f + s2, kP=kinetic(lw, Cw, t, Tw), kpg=fp + sp2)

    out: Dict[str, object] = {"record": "docs/predictions/2026-10-29_ccm_yukawa_coefficient_window.md",
                              "commit": "0295f53a"}

    # ---------------- YC-0 ----------------
    log("=" * 78)
    log("YC-0 (a): fixed t = 0.1 window values for n_max = 10, 12, 14, reproduced to 1e-12")
    ref01 = reference(0.1)
    log(f"   reference t=0.1: kA_inf = {ref01['kA']:.14f} (logged {LOGGED_REF_T01[0]:.14f}), "
        f"kPhi_inf = {ref01['kP']:.14f} (logged {LOGGED_REF_T01[1]:.14f})")
    a_ok = True
    a_rows = []
    for n, (rA, rP, cl) in LOGGED_FIXED_T.items():
        w = window(n, 0.1)
        qA, qP = w["kA"] / ref01["kA"], w["kP"] / ref01["kP"]
        c = 4 * qA / qP
        dev = max(abs(qA - rA), abs(qP - rP), abs(c - cl))
        ok = dev <= 1e-12
        a_ok &= ok
        a_rows.append(dict(n_max=n, kA_ratio=qA, kPhi_ratio=qP, c=c, max_dev_vs_log=dev, pass_=ok))
        log(f"   n_max={n}: kA/kA_inf = {qA:.12f}, kPhi/kPhi_inf = {qP:.12f}, c = {c:.12f}; "
            f"max |dev| vs logged = {dev:.1e}  {'PASS' if ok else 'FAIL'}")
    log(f"   YC-0(a): {'PASS' if a_ok else 'FAIL'}")

    pts = [(tau, n) for tau in (1.0, 4.0) for n in (2, 3, 4, 5)]
    refs = {}
    b_ok, c_ok = True, True
    log("YC-0 (b): tail bounds < 1e-12 relative;  YC-0 (c): |k_pg^inf| < 1e-12")
    for tau, n in pts:
        t = tau / (n + 0.5) ** 2
        r = reference(t)
        refs[(tau, n)] = r
        rb = max(r["bA"] / abs(r["kA"]), r["bP"] / abs(r["kP"]))
        okb = rb < 1e-12
        okc = abs(r["kpg"]) < 1e-12
        b_ok &= okb
        c_ok &= okc
        lA, lP = LOGGED_DEN[(tau, n)]
        log(f"   tau={tau} n_max={n} t={t:.8f}: kA_inf = {r['kA']:.15f} (log dev {abs(r['kA'] - lA):.1e}), "
            f"kPhi_inf = {r['kP']:.15f} (log dev {abs(r['kP'] - lP):.1e}); tail rel {rb:.1e} "
            f"{'PASS' if okb else 'FAIL'}; k_pg_inf = {r['kpg']:+.1e} {'PASS' if okc else 'FAIL'}")
    log(f"   YC-0(b): {'PASS' if b_ok else 'FAIL'};  YC-0(c): {'PASS' if c_ok else 'FAIL'}")
    out["YC0"] = {"a": a_ok, "a_rows": a_rows, "b": b_ok, "c": c_ok}
    if not (a_ok and b_ok and c_ok):
        log("YC-0 FAILED: every item is void. Stopping before any scored point.")
        out["void"] = True
        (HERE / "yc_score.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
        (HERE / "yc_score.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
        return

    # ---------------- scored points ----------------
    log("=" * 78)
    log("Scored points")
    log("   tau  n_max  t           dim   kA/kA_inf        kPhi/kPhi_inf    c                 |c/4-1|      rho_pg       0.1*min(1,miss)")
    rows = []
    for tau, n in pts:
        t = tau / (n + 0.5) ** 2
        r = refs[(tau, n)]
        w = window(n, t)
        qA, qP = w["kA"] / r["kA"], w["kP"] / r["kP"]
        c = 4 * qA / qP
        miss = abs(c / 4 - 1)
        rho = abs(w["kpg"] / w["kA"])
        thr = 0.1 * min(1.0, miss)
        rows.append(dict(tau=tau, n_max=n, t=t, dim=w["dim"], kA=w["kA"], kPhi=w["kP"], kpg=w["kpg"],
                         kA_inf=r["kA"], kPhi_inf=r["kP"], kA_ratio=qA, kPhi_ratio=qP, c=c, miss=miss,
                         rho_pg=rho, yc2_threshold=thr))
        log(f"   {tau:<4} {n:<5}  {t:.8f}  {w['dim']:<4}  {qA:+.10f}  {qP:+.10f}  {c:+.10e}  {miss:.4e}   {rho:.4e}   {thr:.4e}")

    yc1_fail = [r for r in rows if r["miss"] <= 1e-3]
    yc1 = len(yc1_fail) == 0
    yc2_fail = [r for r in rows if r["miss"] > 1e-3 and r["rho_pg"] < r["yc2_threshold"]]
    yc2 = len(yc2_fail) == 0
    yc3 = {}
    for tau in (1.0, 4.0):
        ms = [r["miss"] for r in rows if r["tau"] == tau]
        yc3[tau] = all(ms[i + 1] > ms[i] for i in range(len(ms) - 1))
    log("=" * 78)
    log(f"YC1 (|c/4-1| > 1e-3 at every point): {'HELD' if yc1 else 'MISSED'}"
        + ("" if yc1 else f"  failing points: {[(r['tau'], r['n_max']) for r in yc1_fail]}"))
    log(f"YC2 (rho_pg >= 0.1 min(1, |c/4-1|) where miss > 1e-3): {'HELD' if yc2 else 'MISSED'}; "
        f"failures: {len(yc2_fail)} {[(r['tau'], r['n_max']) for r in yc2_fail]}")
    log(f"YC3 (|c/4-1| strictly increasing in n_max at each tau): tau=1 {'HELD' if yc3[1.0] else 'MISSED'}, "
        f"tau=4 {'HELD' if yc3[4.0] else 'MISSED'}  (reported, not gated)")
    # gate, verbatim
    if yc1 and yc2:
        gate = "GO"
    elif yc1 and 1 <= len(yc2_fail) <= 2:
        gate = "BORDERLINE"
    else:
        gate = "STOP"
    reason = ("YC1 failed (the named interesting outcome)" if not yc1 else
              ("YC2 failed at three or more points" if len(yc2_fail) >= 3 else ""))
    log(f"GATE: {gate} {('-- ' + reason) if reason else ''}")
    out.update(void=False, rows=rows, YC1=yc1, YC2=yc2, YC2_failures=len(yc2_fail),
               YC3={"tau1": yc3[1.0], "tau4": yc3[4.0]}, gate=gate, gate_reason=reason)
    (HERE / "yc_score.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    (HERE / "yc_score.json").write_text(json.dumps(out, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
