# -*- coding: utf-8 -*-
"""WH7 Step-1 probe (2026-06-10): is time metrically visible to the translation
seminorm when the temporal algebra is built from genuine time-dependent
multipliers (Toeplitz-compressed e^{i omega t}) instead of the P45
momentum-diagonal functions g(D_t)?

Setup: circle of circumference T, Fourier window k in {-K..K} (N_t = 2K+1).
  S_q  = P_K M_{e_q} P_K, e_q(t) = exp(2 pi i q t / T)  -> shift matrix.
  U_s  = time translation, diagonal exp(-i omega_k s), omega_k = 2 pi k / T.
  L(F) = sup_s ||U_s F U_s* - F|| / d_T(s, 0)   (translation seminorm;
         the s->0 limit equals ||[D_t, F]|| with D_t = diag(omega_k)).

Five checks:
  A  single-mode exactness:      L(S_q) = 2 pi q / T = Lip(e_q)  (machine precision)
  B  P45 architecture control:   F = g(D_t) diagonal  ->  L(F) = 0 despite g non-constant
  C  kernel condition:           L(F) = 0  <=>  f constant  (band-limited, q_max <= K)
  D  Lipschitz domination:       L(F) <= Lip(f); ratio -> 1 as window K grows
  E  de-compactification:        fixed physical frequency omega, T -> large: L = omega constant

Verdict semantics for WH7 (CLAUDE.md S1.7): A+C+D = compact-carrier temporal
metric REBUILT on the honest Toeplitz algebra (P45's invisibility was the
algebra, not non-compactness); E = visibility does not degrade as the carrier
de-compactifies at fixed bandwidth, i.e. the primary falsifier currently leans
"weakens-to-convention" on the visibility leg, pending the pointed-proper
properness bookkeeping (Step 2, named follow-on).
"""
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260610)


def shift_matrix(K: int, q: int) -> np.ndarray:
    """S_q = P_K M_{e_q} P_K in the Fourier basis: maps mode k -> k+q inside the window."""
    N = 2 * K + 1
    S = np.zeros((N, N), dtype=complex)
    for k in range(-K, K + 1):
        if -K <= k + q <= K:
            S[k + q + K, k + K] = 1.0
    return S


def omega(K: int, T: float) -> np.ndarray:
    return 2 * np.pi * np.arange(-K, K + 1) / T


def compress(coeffs: dict, K: int) -> np.ndarray:
    """F = sum_q f_q S_q for a band-limited f given by Fourier coefficients {q: f_q}."""
    N = 2 * K + 1
    F = np.zeros((N, N), dtype=complex)
    for q, c in coeffs.items():
        F += c * shift_matrix(K, q)
    return F


def translation_seminorm(F: np.ndarray, K: int, T: float, n_grid: int = 400) -> float:
    """L(F) = sup_s ||U_s F U_s* - F||_op / d_T(s,0), d_T = circle geodesic distance.
    Includes the analytic s->0 limit ||[D_t, F]||_op."""
    w = omega(K, T)
    Dt = np.diag(w)
    small_s = np.linalg.norm(Dt @ F - F @ Dt, 2)
    best = small_s
    for s in np.linspace(T / (2 * n_grid), T / 2, n_grid):
        U = np.exp(-1j * w * s)
        Fs = (U[:, None] * F) * np.conj(U)[None, :]
        best = max(best, np.linalg.norm(Fs - F, 2) / s)
    return float(best)


def lipschitz_constant(coeffs: dict, T: float, n_grid: int = 4000) -> float:
    """Lip(f) = max_t |f'(t)| for f = sum_q f_q e_q."""
    t = np.linspace(0, T, n_grid, endpoint=False)
    fp = np.zeros_like(t, dtype=complex)
    for q, c in coeffs.items():
        wq = 2 * np.pi * q / T
        fp += c * 1j * wq * np.exp(1j * wq * t)
    return float(np.max(np.abs(fp)))


def run() -> dict:
    out = {}

    # A — single-mode exactness
    panel_A = []
    for T in (1.0, 2 * np.pi, 10.0):
        for K in (4, 8, 16):
            for q in (1, 2, 3, 4):
                L = translation_seminorm(shift_matrix(K, q), K, T)
                target = 2 * np.pi * q / T
                panel_A.append({"T": T, "K": K, "q": q, "L": L, "target": target,
                                "abs_err": abs(L - target)})
    out["A_single_mode"] = {"max_abs_err": max(r["abs_err"] for r in panel_A),
                            "n_cells": len(panel_A)}

    # B — P45 momentum-diagonal control: g(D_t) commutes with U_s -> L = 0
    panel_B = []
    for K in (4, 8, 16):
        g = np.diag(RNG.normal(size=2 * K + 1))      # non-constant function of D_t
        panel_B.append(translation_seminorm(g, K, T=2 * np.pi))
    out["B_p45_control"] = {"max_L": max(panel_B)}

    # C — kernel condition: L(F)=0 iff f constant (band q_max <= K)
    K, T = 8, 2 * np.pi
    const_L = translation_seminorm(compress({0: 3.7}, K), K, T)
    nonconst_L = []
    for _ in range(20):
        qm = int(RNG.integers(1, K + 1))
        coeffs = {0: RNG.normal()}
        for q in range(1, qm + 1):
            c = RNG.normal() + 1j * RNG.normal()
            coeffs[q] = c
            coeffs[-q] = np.conj(c)                  # real f
        nonconst_L.append(translation_seminorm(compress(coeffs, K), K, T))
    out["C_kernel"] = {"const_L": const_L, "min_nonconst_L": min(nonconst_L)}

    # D — Lipschitz domination + window convergence
    fcoeffs = {1: 0.5, -1: 0.5, 2: 0.25, -2: 0.25, 3: 0.125, -3: 0.125}
    lip = lipschitz_constant(fcoeffs, T)
    ratios = {}
    for K in (3, 4, 8, 16, 32):
        L = translation_seminorm(compress(fcoeffs, K), K, T)
        ratios[K] = L / lip
        assert L <= lip * (1 + 1e-9), "Lipschitz domination violated"
    out["D_lipschitz"] = {"Lip": lip, "ratio_by_K": ratios,
                          "monotone_to_1": all(ratios[a] <= ratios[b] + 1e-12
                                               for a, b in zip([3, 4, 8, 16], [4, 8, 16, 32]))}

    # E — de-compactification at fixed physical frequency omega = 2*pi
    panel_E = {}
    for T in (1.0, 2.0, 4.0, 8.0, 16.0):
        q = int(round(T))                            # omega_q = 2 pi q / T = 2 pi
        K = 2 * q + 2
        panel_E[T] = translation_seminorm(shift_matrix(K, q), K, T)
    out["E_decompactification"] = {"L_by_T": panel_E,
                                   "max_dev_from_2pi": max(abs(v - 2 * np.pi)
                                                           for v in panel_E.values())}
    return out


if __name__ == "__main__":
    res = run()
    (ROOT / "data" / "wh7_toeplitz_probe.json").write_text(
        json.dumps(res, indent=1, default=float), encoding="utf-8")
    a = res["A_single_mode"]["max_abs_err"]
    b = res["B_p45_control"]["max_L"]
    c0, c1 = res["C_kernel"]["const_L"], res["C_kernel"]["min_nonconst_L"]
    e = res["E_decompactification"]["max_dev_from_2pi"]
    print(f"A single-mode exactness : max |L - 2*pi*q/T| = {a:.2e}  -> {'PASS' if a < 1e-10 else 'FAIL'}")
    print(f"B P45 control (L == 0)  : max L = {b:.2e}              -> {'PASS' if b < 1e-12 else 'FAIL'}")
    print(f"C kernel condition      : const {c0:.1e}, min nonconst {c1:.3f} -> "
          f"{'PASS' if c0 < 1e-12 and c1 > 1e-6 else 'FAIL'}")
    print(f"D Lipschitz domination  : ratios {res['D_lipschitz']['ratio_by_K']}, "
          f"monotone {res['D_lipschitz']['monotone_to_1']}")
    print(f"E fixed-omega T-sweep   : max |L - 2pi| = {e:.2e}      -> {'PASS' if e < 1e-10 else 'FAIL'}")
