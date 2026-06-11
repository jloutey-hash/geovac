# -*- coding: utf-8 -*-
"""B3 Phase 1: boost/modular-flow seminorm probe (2026-06-10, v3.113.0).

Question (Paper 45 Q1, updated): can the wedge boost -- the modular flow of the
four-witness theorem (Paper 42, K = J_polar doubled to integer spectrum) -- serve as
the TEMPORAL generator of a quantum-metric structure, so that Lorentzian signature
enters through the generator of the seminorm?

Substrate: the B1 Peter-Weyl window j <= 1 (dim 14) with the full reachable
multiplier system, bands b in {0, 1/2, 1, 3/2, 2} (predicted rank 55 = sum (2b+1)^2,
the n_max = 3 system dimension of Paper 38). Boost generator K = 2 J_z (Paper 42
two_m_j convention: doubling makes the conjugation spectrum integer).

Checks:
  T0  system rank = 55; ad-invariance of the system under all generators
  T1  four-witness compatibility: sigma_{2pi}(F) = F bit-exact for ALL F;
      band-parity grading at HALF period: sigma_pi(F) = (-1)^{2b} F bit-exact
      (finite-cutoff spin-statistics shadow of the modular flow)
  T2  boost-alone kernel: dim ker(ad_K | system) = 9 = sum_{b in {0,1,2}} (2b+1)
      -- the boost-invariant multipliers; structured partial kernel, NOT the P45
      annihilation (dim 55) and NOT the metric kernel condition (dim 1)
  T3  frame completion: joint kernel of {ad_Jx, ad_Jy, ad_Jz} = C 1 (dim 1)
  T4  causal classifier: Q(F) = ||[J_z, F]||^2 - ||[J_x, F]||^2 - ||[J_y, F]||^2
      per band-weight element (Frobenius and operator norms); idealized symbol
      prediction sign(2 m'^2 - b(b+1)): b = 1 top weight on the cone, b >= 3/2
      top weights timelike, m' = 0 spacelike.

Companion: debug/wh7_b1_joint_product_gh.py (machinery), geovac/modular_hamiltonian.py
(Paper 42 conventions). Frozen falsifier: tests/test_wh7_b3_boost.py.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import wh7_b1_joint_product_gh as b1  # noqa: E402

BANDS = [0.0, 0.5, 1.0, 1.5, 2.0]


def dcache(b):
    if b in b1.DCACHE:
        return b1.DCACHE[b]
    return np.array([b1.wigner_D(b, *p) for p in b1.PTS])


def build_system():
    """All band-multiplier compressions C^b_{m'm}; returns ops, labels."""
    ops, labels = [], []
    for b in BANDS:
        D = dcache(b)
        dim = int(round(2 * b + 1))
        for a in range(dim):           # row index: m' = b - a (left weight)
            for c in range(dim):
                f = np.sqrt(2 * b + 1) * D[:, a, c]
                ops.append(b1.compress_spatial(f))
                labels.append((b, b - a, b - c))
    return ops, labels


def orthobasis(ops):
    X = np.array([op.flatten() for op in ops]).T          # (196, n_ops)
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    rank = int(np.sum(s > 1e-10 * s[0]))
    return U[:, :rank], rank


def ad_matrix(G, U):
    """Matrix of ad_G = [G, .] in the orthonormal system coordinates U;
    also returns max projection residual (system-invariance check)."""
    n = U.shape[1]
    M = np.zeros((n, n), dtype=complex)
    resid = 0.0
    for jcol in range(n):
        B = U[:, jcol].reshape(b1.NW, b1.NW)
        C = (G @ B - B @ G).flatten()
        coeff = U.conj().T @ C
        M[:, jcol] = coeff
        resid = max(resid, float(np.linalg.norm(C - U @ coeff)))
    return M, resid


def nullity(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False)
    smax = s[0] if s[0] > 0 else 1.0
    return int(np.sum(s < tol * smax))


def run():
    out = {}
    Jx, Jy, Jz = b1.A_GENS
    K_boost = 2.0 * Jz                                    # Paper 42 two_m_j convention

    ops, labels = build_system()
    U, rank = orthobasis(ops)
    out["T0_rank"] = {"rank": rank, "predicted": 55}

    # T1 — flow closure and half-period grading (diagonal phases, exact)
    phases_2pi = np.exp(1j * 2 * np.pi * np.diag(K_boost))
    phases_pi = np.exp(1j * np.pi * np.diag(K_boost))
    max_2pi, max_grade = 0.0, 0.0
    for op, (b, mp, m) in zip(ops, labels):
        s2 = (phases_2pi[:, None] * op) * np.conj(phases_2pi)[None, :]
        max_2pi = max(max_2pi, float(np.max(np.abs(s2 - op))))
        sp = (phases_pi[:, None] * op) * np.conj(phases_pi)[None, :]
        sign = (-1) ** int(round(2 * b))
        max_grade = max(max_grade, float(np.max(np.abs(sp - sign * op))))
    out["T1_flow"] = {"sigma_2pi_residual": max_2pi,
                      "band_parity_grading_residual": max_grade}

    # T2 — boost-alone kernel
    MK, rK = ad_matrix(K_boost, U)
    out["T2_boost_kernel"] = {"dim_ker": nullity(MK), "predicted": 9,
                              "system_invariance_residual": rK}

    # T3 — frame completion
    Ms = [ad_matrix(G, U)[0] for G in (Jx, Jy, Jz)]
    stacked = np.vstack(Ms)
    out["T3_frame_kernel"] = {"dim_ker": nullity(stacked), "predicted": 1}

    # T4 — causal classifier per (b, m') class
    classes = {}
    for op, (b, mp, m) in zip(ops, labels):
        key = (b, abs(mp))
        for norm, tag in ((lambda A: np.linalg.norm(A, 2), "op"),
                          (lambda A: np.linalg.norm(A), "fro")):
            lz = norm(Jz @ op - op @ Jz)
            lx = norm(Jx @ op - op @ Jx)
            ly = norm(Jy @ op - op @ Jy)
            denom = lz ** 2 + lx ** 2 + ly ** 2
            q = (lz ** 2 - lx ** 2 - ly ** 2) / denom if denom > 1e-20 else 0.0
            classes.setdefault((key, tag), []).append(q)
    table = []
    for ((b, amp), tag), qs in sorted(classes.items()):
        ideal = 2 * amp ** 2 - b * (b + 1)
        table.append({"b": b, "abs_mp": amp, "norm": tag,
                      "q_mean": float(np.mean(qs)), "q_min": float(np.min(qs)),
                      "q_max": float(np.max(qs)),
                      "ideal_sign": int(np.sign(ideal)) if abs(ideal) > 1e-12 else 0})
    out["T4_causal_classes"] = table
    return out


if __name__ == "__main__":
    res = run()
    (ROOT / "data" / "wh7_b3_boost_probe.json").write_text(
        json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(f"T0 system rank            : {res['T0_rank']['rank']} (predicted 55)")
    t1 = res["T1_flow"]
    print(f"T1 sigma_2pi closure      : {t1['sigma_2pi_residual']:.2e}; "
          f"half-period band grading : {t1['band_parity_grading_residual']:.2e}")
    t2 = res["T2_boost_kernel"]
    print(f"T2 boost-alone kernel     : dim {t2['dim_ker']} (predicted 9); "
          f"invariance resid {t2['system_invariance_residual']:.2e}")
    print(f"T3 frame kernel           : dim {res['T3_frame_kernel']['dim_ker']} (predicted 1)")
    print("T4 causal classes (op norm):")
    for r in res["T4_causal_classes"]:
        if r["norm"] == "op":
            verdict = ("timelike" if r["q_min"] > 1e-10 else
                       "spacelike" if r["q_max"] < -1e-10 else "mixed/null")
            print(f"   b={r['b']:.1f} |m'|={r['abs_mp']:.1f}: "
                  f"q in [{r['q_min']:+.3f}, {r['q_max']:+.3f}]  ideal sign "
                  f"{r['ideal_sign']:+d}  -> {verdict}")
