"""Frozen falsifiers for Paper 38's unconditional theorem (2026-06-10).

Guards the two load-bearing facts of the translation-seminorm
framework (Paper 38 §sec:named_gaps; sprint memo
debug/sprint_p38_g1g2_phaseA_memo.md):

  T1. Per-band injectivity (premise of the kernel proposition):
      on the chirality-doubled spinor multiplier system, every
      harmonic band N <= 2 n_max - 1 has full vec-rank N^2.
  T2. The lifted-state smoothing identity (scalar sector, exact):
      upsilon(P M_f P)(g) = sigma_J f(g) with sigma the unit-step
      Clebsch-Gordan multiplier symbol — checked on the SU(2)
      Peter-Weyl window j <= 1/2 (the n_max = 2 scalar window),
      where sigma = (1, sqrt(2)/3, 2/9), the published values.

If T1 ever fails the kernel proposition's premise is wrong; if T2
fails the dual map's structural identity (and with it the
almost-inverse bookkeeping of the unconditional theorem) is wrong.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction

import numpy as np
import pytest

from geovac.full_dirac_operator_system import FullDiracTruncatedOperatorSystem


@pytest.mark.parametrize("n_max", [2, 3])
def test_per_band_injectivity(n_max):
    osys = FullDiracTruncatedOperatorSystem(n_max)
    bands = defaultdict(list)
    for (N, L, M), mat in zip(osys.multiplier_labels, osys.multiplier_matrices):
        bands[N].append(mat)
    for N, mats in bands.items():
        stack = np.stack([m.ravel() for m in mats])
        sv = np.linalg.svd(stack, compute_uv=False)
        rank = int(np.sum(sv > 1e-10 * sv[0]))
        assert len(mats) == N * N, (n_max, N, len(mats))
        assert rank == N * N, (n_max, N, rank)


def _su2_generators(j):
    ms = [float(Fraction(int(round(2 * j)), 2)) - i
          for i in range(int(round(2 * j)) + 1)]
    d = len(ms)
    Jz = np.diag(ms).astype(complex)
    Jp = np.zeros((d, d), dtype=complex)
    for i in range(1, d):
        m = ms[i]
        Jp[i - 1, i] = np.sqrt(j * (j + 1) - m * (m + 1))
    Jy = (Jp - Jp.conj().T) / 2j
    return Jz, Jy


def _wigner_D(j, a, b, c):
    from scipy.linalg import expm
    Jz, Jy = _su2_generators(j)
    return expm(-1j * a * Jz) @ expm(-1j * b * Jy) @ expm(-1j * c * Jz)


def test_lifted_state_smoothing_identity_scalar_window():
    """Window W = {0, 1/2} (n=2 scalar levels), bands J = 1/2, 1.

    Exact CG data at this size (no sympy needed):
      basis: (0,0,0), (1/2,a,b) a,b in {+-1/2}  -> dim 5
      sigma_{1/2} = 2*sqrt(2)/6 = sqrt(2)/3,  sigma_1 = 2/9... note:
      at the n=2 window the published symbols are (1, sqrt2/3, 2/9)
      for J = (0, 1/2, 1).
    """
    from sympy import S
    from sympy.physics.quantum.cg import CG as _CG

    js = [Fraction(0), Fraction(1, 2)]
    Z = sum(int(2 * j) + 1 for j in js)  # = 3

    def jm(j):
        k = int(round(2 * j))
        return [Fraction(k, 2) - i for i in range(k + 1)]

    index, idx = {}, 0
    for j in js:
        for a in jm(j):
            for b in jm(j):
                index[(j, a, b)] = idx
                idx += 1
    dim = idx  # 5

    def cg(j1, m1, j2, m2, j3, m3):
        return float(_CG(S(j1), S(m1), S(j2), S(m2), S(j3), S(m3)).doit())

    def multiplier(J, A, B):
        M = np.zeros((dim, dim), dtype=complex)
        for j2 in js:
            for j1 in js:
                if abs(J - j2) > j1 or j1 > J + j2:
                    continue
                if (j1 + j2 - J) % 1 != 0:
                    # CG coefficients vanish here anyway; skip for speed
                    continue
                pref = np.sqrt((int(2 * j2) + 1) / (int(2 * j1) + 1))
                for a2 in jm(j2):
                    a1 = A + a2
                    if abs(a1) > j1:
                        continue
                    ca = cg(J, A, j2, a2, j1, a1)
                    if ca == 0.0:
                        continue
                    for b2 in jm(j2):
                        b1 = B + b2
                        if abs(b1) > j1:
                            continue
                        cb = cg(J, B, j2, b2, j1, b1)
                        if cb == 0.0:
                            continue
                        M[index[(j1, a1, b1)], index[(j2, a2, b2)]] = \
                            pref * ca * cb
        return M

    def U(alpha, beta, gamma_):
        from scipy.linalg import block_diag
        blocks = []
        for j in js:
            Dg = _wigner_D(float(j), alpha, beta, gamma_)
            blocks.append(np.kron(Dg.conj(), np.eye(Dg.shape[0])))
        return block_diag(*blocks)

    # Fejer vector: 1 on diagonal elements
    xi = np.zeros(dim, dtype=complex)
    for j in js:
        for a in jm(j):
            xi[index[(j, a, a)]] = 1.0
    xi /= np.linalg.norm(xi)

    rng = np.random.default_rng(11)
    grid = [tuple(rng.uniform(0, 2 * np.pi, 3)) for _ in range(8)]

    expected_sigma = {Fraction(1, 2): np.sqrt(2) / 3, Fraction(1): 2.0 / 9.0}
    for J, sig in expected_sigma.items():
        ms = jm(J)
        A, B = ms[0], ms[-1]
        M = multiplier(J, A, B)
        for g in grid:
            xv = U(*g) @ xi
            ups = np.vdot(xv, M @ xv)
            D = _wigner_D(float(J), *g)
            fg = D[0, -1]
            assert abs(ups - sig * fg) < 1e-12, (float(J), g)
