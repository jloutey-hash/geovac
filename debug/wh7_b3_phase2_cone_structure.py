# -*- coding: utf-8 -*-
"""B3 Phase 2: exact cone ratios + causal-form signature + reverse-triangle probe
(2026-06-10, v3.114.0).

Part A (exact, discrete-skeleton route -- sympy Clebsch-Gordan, no PSLQ):
  Frobenius matrix elements of band compressions on the j <= 1 window are CG
  products:  <psi^{j1}|D^b_{m'm}|psi^{j2}> = sqrt((2b+1)(2j2+1)/(2j1+1)) CG' CG,
  so  N(b,m',m) := ||C^b_{m'm}||_F^2 = (2b+1) sum_{j1 j2} ((2j2+1)/(2j1+1))
                    S(j1,j2,b,m') S(j1,j2,b,m),
      S(j1,j2,b,w) := sum_mu CG(j2 mu; b w | j1 mu+w)^2     (exact rationals).
  Weight-grading kills the J_x/J_y cross terms in Frobenius norm, so
      q_F(b,m',m) = [m'^2 N - (c+^2 N+ + c-^2 N-)/2] / [m'^2 N + (c+^2 N+ + c-^2 N-)/2]
  is an EXACT RATIONAL for every system element. Verified against float quadrature.

Part B (numeric): inertia (n+, n0, n-) of the causal quadratic form
  Q(F) = ||[Jz,F]||_F^2 - ||[Jx,F]||_F^2 - ||[Jy,F]||_F^2  on the 55-dim system and
  per (b,|m'|) class span. Lorentzian signature (1, k) on a timelike sector implies
  the reverse (Minkowski) triangle inequality there -- the operator-level twin
  paradox -- as a THEOREM (reverse Cauchy-Schwarz), not a sampling claim.
  Plus: op-norm tau super-additivity sampling as an independent check, and op-norm
  top-weight ratio identifications (b=1 null member; ||C_top||^2/||C_top-1||^2).

Companion: debug/wh7_b3_boost_seminorm_probe.py (Phase 1). Frozen falsifier:
tests/test_wh7_b3_phase2.py.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from sympy import Rational, S
from sympy.physics.quantum.cg import CG

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import wh7_b1_joint_product_gh as b1  # noqa: E402
import wh7_b3_boost_seminorm_probe as b3  # noqa: E402

WINDOW = [S(0), S(1) / 2, S(1)]
BANDS = [S(0), S(1) / 2, S(1), S(3) / 2, S(2)]
RNG = np.random.default_rng(20260612)


# ---------------- Part A: exact rational cone table --------------------------
def S_sum(j1, j2, b, w):
    tot = S(0)
    mu = -j2
    while mu <= j2:
        if abs(mu + w) <= j1:
            c = CG(j2, mu, b, w, j1, mu + w).doit()
            tot += c ** 2
        mu += 1
    return sp.nsimplify(tot)


def N_exact(b, mp, m):
    if abs(mp) > b or abs(m) > b:
        return S(0)
    tot = S(0)
    for j1 in WINDOW:
        for j2 in WINDOW:
            if abs(j1 - j2) <= b <= j1 + j2:
                tot += (2 * j2 + 1) / (2 * j1 + 1) * S_sum(j1, j2, b, mp) * S_sum(j1, j2, b, m)
    return sp.simplify((2 * b + 1) * tot)


def q_exact(b, mp, m):
    N0 = N_exact(b, mp, m)
    if N0 == 0:
        return None
    cp2 = b * (b + 1) - mp * (mp + 1)
    cm2 = b * (b + 1) - mp * (mp - 1)
    legs = (cp2 * N_exact(b, mp + 1, m) + cm2 * N_exact(b, mp - 1, m)) / 2
    tz = mp ** 2 * N0
    return sp.simplify((tz - legs) / (tz + legs))


def part_A():
    table = {}
    max_dev = 0.0
    Jx, Jy, Jz = b1.A_GENS
    for b in BANDS:
        dim = int(2 * b + 1)
        for a in range(dim):
            mp = b - a
            for c in range(dim):
                m = b - c
                qe = q_exact(b, mp, m)
                key = f"b={b},mp={mp},m={m}"
                table[key] = str(qe)
                # float cross-check
                D = b3.dcache(float(b))
                f = np.sqrt(float(2 * b + 1)) * D[:, a, c]
                C = b1.compress_spatial(f)
                lz = np.linalg.norm(Jz @ C - C @ Jz)
                lx = np.linalg.norm(Jx @ C - C @ Jx)
                ly = np.linalg.norm(Jy @ C - C @ Jy)
                den = lz ** 2 + lx ** 2 + ly ** 2
                qn = (lz ** 2 - lx ** 2 - ly ** 2) / den if den > 1e-18 else None
                if qe is not None and qn is not None:
                    max_dev = max(max_dev, abs(float(qe) - qn))
    return table, max_dev


# ---------------- Part B: signature + reverse triangle -----------------------
def part_B():
    Jx, Jy, Jz = b1.A_GENS
    ops, labels = b3.build_system()
    U, rank = b3.orthobasis(ops)
    n = U.shape[1]
    basis_ops = [U[:, k].reshape(b1.NW, b1.NW) for k in range(n)]

    def gram(G):
        comms = [ (G @ B - B @ G).flatten() for B in basis_ops ]
        Cm = np.array(comms)
        return Cm.conj() @ Cm.T

    Q = gram(Jz) - gram(Jx) - gram(Jy)
    Q = (Q + Q.conj().T) / 2
    ev = np.linalg.eigvalsh(Q)
    tol = 1e-9 * max(1.0, np.max(np.abs(ev)))
    inertia_global = (int(np.sum(ev > tol)), int(np.sum(np.abs(ev) <= tol)),
                      int(np.sum(ev < -tol)))

    # per-(b,|m'|) class spans
    per_class = {}
    for key in sorted({(float(b), abs(float(mp))) for (b, mp, m) in labels}):
        idx = [i for i, (b, mp, m) in enumerate(labels)
               if (float(b), abs(float(mp))) == key]
        X = np.array([ops[i].flatten() for i in idx]).T
        Uc, sc, _ = np.linalg.svd(X, full_matrices=False)
        rc = int(np.sum(sc > 1e-10 * sc[0]))
        Uc = Uc[:, :rc]
        sub_ops = [Uc[:, k].reshape(b1.NW, b1.NW) for k in range(rc)]

        def gram_c(G):
            Cm = np.array([(G @ B - B @ G).flatten() for B in sub_ops])
            return Cm.conj() @ Cm.T

        Qc = gram_c(Jz) - gram_c(Jx) - gram_c(Jy)
        Qc = (Qc + Qc.conj().T) / 2
        evc = np.linalg.eigvalsh(Qc)
        tc = 1e-9 * max(1.0, np.max(np.abs(evc)))
        per_class[str(key)] = {
            "dim": rc,
            "inertia": (int(np.sum(evc > tc)), int(np.sum(np.abs(evc) <= tc)),
                        int(np.sum(evc < -tc)))}

    # op-norm tau super-additivity sampling on timelike classes
    def legs_op(F):
        lz = np.linalg.norm(Jz @ F - F @ Jz, 2)
        lx = np.linalg.norm(Jx @ F - F @ Jx, 2)
        ly = np.linalg.norm(Jy @ F - F @ Jy, 2)
        return lz ** 2 - lx ** 2 - ly ** 2

    def tau(F):
        v = legs_op(F)
        return np.sqrt(v) if v > 0 else None

    violations, trials = 0, 0
    worst = 0.0
    for key in [(1.5, 1.5), (2.0, 2.0)]:
        idx = [i for i, (b, mp, m) in enumerate(labels)
               if (float(b), float(mp)) == key]          # co-oriented: same m' sign
        for _ in range(60):
            w1 = RNG.random(len(idx))
            w2 = RNG.random(len(idx))
            F = sum(a * ops[i] for a, i in zip(w1, idx))
            G = sum(a * ops[i] for a, i in zip(w2, idx))
            tF, tG, tFG = tau(F), tau(G), tau(F + G)
            if tF is None or tG is None or tFG is None:
                continue
            trials += 1
            gap = tFG - (tF + tG)
            worst = min(worst, gap)
            if gap < -1e-9 * max(1.0, tFG):
                violations += 1

    # op-norm top-weight ratio identifications
    ratios = {}
    for b in (0.5, 1.0, 1.5, 2.0):
        D = b3.dcache(b)
        dim = int(round(2 * b + 1))
        r_top, r_next = [], []
        for c in range(dim):
            f_t = np.sqrt(2 * b + 1) * D[:, 0, c]         # m' = b
            f_n = np.sqrt(2 * b + 1) * D[:, 1, c]         # m' = b-1
            Ct = b1.compress_spatial(f_t)
            Cn = b1.compress_spatial(f_n)
            r_top.append(np.linalg.norm(Ct, 2) ** 2)
            r_next.append(np.linalg.norm(Cn, 2) ** 2)
        ratios[b] = [t / nx for t, nx in zip(r_top, r_next)]

    # b=1 top-weight null member: ||C_{1,m}||_op vs ||C_{0,m}||_op
    D1 = b3.dcache(1.0)
    null_resid = {}
    for c, m in zip(range(3), (1, 0, -1)):
        Ct = b1.compress_spatial(np.sqrt(3) * D1[:, 0, c])
        C0 = b1.compress_spatial(np.sqrt(3) * D1[:, 1, c])
        null_resid[m] = float(legs_op(Ct))
    return {"inertia_global": inertia_global, "per_class": per_class,
            "tau_superadd": {"trials": trials, "violations": violations,
                             "worst_gap": worst},
            "top_ratio_op": {str(k): v for k, v in ratios.items()},
            "b1_top_null_residual": null_resid}


if __name__ == "__main__":
    tab, dev = part_A()
    resB = part_B()
    out = {"A_exact_q_table": tab, "A_max_float_dev": dev, **resB}
    (ROOT / "data" / "wh7_b3_phase2_cone.json").write_text(
        json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(f"A: exact rational q for all elements; max |exact - float| = {dev:.2e}")
    shown = set()
    for k, v in tab.items():
        b = k.split(",")[0]
        mp = k.split(",")[1]
        if v != "None" and (b, mp) not in shown:
            shown.add((b, mp))
            print(f"   {k}: q_F = {v}")
    print(f"B: global inertia (n+, n0, n-) = {resB['inertia_global']}")
    for k, v in resB["per_class"].items():
        print(f"   class {k}: dim {v['dim']}, inertia {v['inertia']}")
    t = resB["tau_superadd"]
    print(f"B: tau super-additivity: {t['violations']}/{t['trials']} violations, "
          f"worst gap {t['worst_gap']:.2e}")
    print(f"B: op-norm top ratios ||C_b,b||^2/||C_b,b-1||^2: "
          f"{ {k: [round(x,6) for x in v] for k,v in resB['top_ratio_op'].items()} }")
    print(f"B: b=1 top null residual per m: "
          f"{ {k: f'{v:.2e}' for k,v in resB['b1_top_null_residual'].items()} }")
