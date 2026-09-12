"""Backing for Paper 60's [PRIOR ART] / [MEASURED] / [SYMBOLIC] passage on the
Kac-Murdock-Szego attribution of eq:sigma_law (added 2026-09-11).

Written as a SEPARATE pass from the paper edit it protects, per CLAUDE.md
Sec. 9 ("guard-writing is a separate, separately-reviewed activity").  The
question asked of each guard below is *what wrong answer would this accept?*,
named explicitly in each docstring, and each is fire-tested against that
wrong answer via debug/qa/fire_test.py.

Claims backed:
  1. c_1 = pi^2 in the KMS asymptotic lam_min ~ (c_alpha / n^{2 alpha}) b(1).
  2. Our symbol is the alpha = 1 case with curvature b(1) = (kR)^2 / 24.
  3. (1) x (2) reproduces eq:sigma_law exactly -- so the law is KMS, not ours.
  4. The chi -> 0 chirp gives |c_j| ~ j^{-5/4}, hence sum_j j |c_j| DIVERGES,
     so the Boettcher-Widom smoothness hypothesis fails on our symbol.
  5. The restatement 1 - sigma_max = (1/6) (R / L_max)^2, L_max = 2n/(pi k).
  6. Proposition D: a block-diagonal congruence cannot orthogonalize a metric
     that is not itself block diagonal -- the l-selection loss is independent
     of conditioning.
"""
import numpy as np
import pytest
import sympy as sp

KR_VALUES = (1.0, 2.0, 5.0)


# ---------------------------------------------------------------- KMS constant
def test_kms_constant_c1_is_pi_squared():
    """c_1 = pi^2, checked in closed form on the canonical symbol.

    For a(theta) = |1 - e^{i theta}|^2 = 4 sin^2(theta/2) (so b == 1, b(1) = 1),
    T_n(a) is the standard tridiagonal (2, -1) matrix with eigenvalues known
    exactly: 4 sin^2(k pi / (2(n+1))).  KMS then says lam_min (n+1)^2 -> c_1.

    WRONG ANSWER REJECTED: any c_1 other than pi^2 -- in particular pi^2/2,
    2 pi^2, or 4 (the naive 'curvature' guess).  The tolerance below is 3e-4
    at n = 10000, far tighter than the gap to any of those.
    """
    n = 10_000
    lam_min = 4 * np.sin(np.pi / (2 * (n + 1))) ** 2
    c1 = lam_min * (n + 1) ** 2
    assert abs(c1 - np.pi**2) < 3e-4, f"c_1 = {c1}, expected pi^2 = {np.pi**2}"
    # and it must CONVERGE, not merely sit near pi^2 at one n
    prev = 4 * np.sin(np.pi / (2 * 101)) ** 2 * 101**2
    assert abs(c1 - np.pi**2) < abs(prev - np.pi**2) / 100


def test_symbol_curvature_b1_is_kR_squared_over_24():
    """b(1) = (kR)^2 / 24, exactly, symbolically in kR.

    Paper 60's symbol is a(chi) = j0(kR cot(chi/2)); the small quantity at the
    IR pole is theta = pi - chi, and 1 - a = b(1) theta^2 + O(theta^4).

    WRONG ANSWER REJECTED: any other rational multiple of (kR)^2 -- 1/6 (the
    bare sinc expansion, forgetting cot(chi/2) -> theta/2), 1/12, 1/48 -- and
    any wrong POWER of kR.  Checked as an exact symbolic equality, so there is
    no tolerance to hide in.
    """
    th, kR = sp.symbols("theta kR", positive=True)
    x = kR / sp.tan((sp.pi - th) / 2)
    one_minus_a = 1 - sp.sin(x) / x
    series = sp.series(one_minus_a, th, 0, 4).removeO()
    b1 = sp.simplify(series.coeff(th, 2))
    assert sp.simplify(b1 - kR**2 / 24) == 0, f"b(1) = {b1}, expected kR**2/24"
    # the leading order must be theta^2 (alpha = 1), not theta^0 or theta^4
    assert sp.simplify(series.coeff(th, 0)) == 0
    assert sp.simplify(b1) != 0


def test_kms_product_reproduces_sigma_law():
    """c_1 b(1) / n^2 IS eq:sigma_law, i.e. the law is KMS and not ours.

    WRONG ANSWER REJECTED: the claim that eq:sigma_law's constant differs from
    the KMS prediction by any factor != 1 (which is what a 'we derived this
    independently' reading would need).  Exact symbolic difference.
    """
    n, kR = sp.symbols("n kR", positive=True)
    kms = sp.pi**2 * (kR**2 / 24) / n**2          # c_1 * b(1) / n^2
    paper = sp.pi**2 * kR**2 / (24 * n**2)        # eq:sigma_law as printed
    assert sp.simplify(kms - paper) == 0


# ------------------------------------------------------------- the UV residue
def _symbol_cosine_coeff(j, kR, Smax=2000.0, gl_n=48):
    """c_j = (1/pi) int_0^pi cos(j chi) a(chi) dchi, via s = kR cot(chi/2).

    Fixed-order Gauss-Legendre on panels sized so BOTH phases advance <= pi per
    panel.  Deterministic: no adaptive solver that can silently not converge
    (mpmath quadosc does exactly that here past j ~ 512 and returns values that
    GROW with j, which is impossible for a continuous symbol).
    """
    xg, wg = np.polynomial.legendre.leggauss(gl_n)
    edges, s = [0.0], 0.0
    while s < Smax:
        s += min(np.pi / 2, np.pi / max(2 * j * kR / (kR * kR + s * s), 1e-300))
        edges.append(min(s, Smax))
    e = np.asarray(edges)
    mid, half = 0.5 * (e[:-1] + e[1:]), 0.5 * (e[1:] - e[:-1])
    sv = (mid[:, None] + half[:, None] * xg[None, :]).ravel()
    wv = (half[:, None] * wg[None, :]).ravel()
    integ = np.cos(2 * j * np.arctan2(kR, sv)) * np.sinc(sv / np.pi) \
        * 2 * kR / (kR * kR + sv * sv)
    return float(np.dot(wv, integ) / np.pi)


@pytest.mark.parametrize("kR", KR_VALUES)
def test_chirp_envelope_exponent_is_five_fourths(kR):
    """|c_j| ~ j^{-5/4}: stationary phase on the chi -> 0 chirp, parameter-free.

    The prediction carries a modulating sin(2 sqrt(2 kR j) + pi/4), so a
    pointwise test would fail at its zeros; the ENVELOPE is tested by taking a
    running maximum over a window of j, which the modulation cannot suppress.

    WRONG ANSWER REJECTED: exponents 1 (a plain jump/BV symbol), 3/2, or 2 (a
    Lipschitz symbol) -- each 0.25 or more away from 1.25, against a tolerance
    of 0.08.  Also rejects exponential decay, which would drive the fitted
    slope steeply negative and out of band.
    """
    js = np.unique(np.round(np.logspace(np.log10(64), np.log10(6000), 40)).astype(int))
    c = np.array([abs(_symbol_cosine_coeff(int(j), kR)) for j in js])
    # envelope: running max over a half-decade window kills the modulation
    env, jj = [], []
    for lo in np.logspace(np.log10(64), np.log10(2000), 8):
        m = (js >= lo) & (js < lo * 2.2)
        if m.sum() >= 2:
            env.append(c[m].max())
            jj.append(js[m][np.argmax(c[m])])
    slope = np.polyfit(np.log(jj), np.log(env), 1)[0]
    assert abs(slope + 1.25) < 0.08, f"envelope exponent {slope:.3f}, expected -1.25"


def test_bottcher_widom_smoothness_hypothesis_fails():
    """sum_j j |c_j| DIVERGES, so the hypothesis proving the constant fails.

    With |c_j| ~ j^{-5/4}, the partial sums of j|c_j| grow like K^{3/4}.  The
    test asserts growth, not a value: doubling the cut must raise the partial
    sum by a factor clearly above 1 (2^{3/4} = 1.68).

    WRONG ANSWER REJECTED: convergence -- i.e. the claim that our symbol DOES
    satisfy the hypothesis.  A convergent sum would give a ratio -> 1.0; the
    assertion below requires > 1.35, which no convergent series can sustain
    across two successive doublings.
    """
    kR = 2.0
    js = np.arange(1, 1601)
    w = np.array([abs(_symbol_cosine_coeff(int(j), kR)) for j in js]) * js
    partial = [w[:k].sum() for k in (200, 400, 800, 1600)]
    ratios = [partial[i + 1] / partial[i] for i in range(3)]
    assert all(r > 1.35 for r in ratios), f"partial-sum ratios {ratios} -- looks convergent"


# --------------------------------------------------- the wavelength restatement
def test_Lmax_restatement_matches_sigma_law():
    """1 - sigma_max = (1/6)(R/L_max)^2 with L_max = 2n/(pi k) IS eq:sigma_law.

    WRONG ANSWER REJECTED: any other prefactor (1/24, 1/2) or any other
    L_max convention (n/(pi k), 2n/k) -- all shift the identity by a factor
    that this exact symbolic difference catches.
    """
    n, k, R = sp.symbols("n k R", positive=True)
    L_max = 2 * n / (sp.pi * k)
    restated = sp.Rational(1, 6) * (R / L_max) ** 2
    law = sp.pi**2 * (k * R) ** 2 / (24 * n**2)
    assert sp.simplify(restated - law) == 0


# ------------------------------------------------------------- Proposition D
def test_prop_d_block_diagonal_congruence_preserves_block_structure():
    """If X is block diagonal and X^dag S X is block diagonal, so is S.

    Forward direction, constructively: S block diagonal => S^{-1/2} block
    diagonal (this is the m-selection half, which SURVIVES).  Contrapositive,
    constructively: for an S that is NOT block diagonal, no block-diagonal X
    makes X^dag S X block diagonal -- sampled over random block-diagonal X.

    WRONG ANSWER REJECTED: the reading that l-selection loss is a CONDITIONING
    effect that relaxes as cond(S) -> 1.  The second half of this test uses an
    S whose off-block coupling is scaled to eps = 1e-6, i.e. cond(S) = 1 + 2e-6,
    essentially perfectly conditioned -- and the block structure is destroyed
    just the same.  A guard that only tested an ill-conditioned S would accept
    that wrong reading.
    """
    rng = np.random.default_rng(6002)
    sizes = (3, 4)          # two "l blocks"
    N = sum(sizes)

    def is_block_diagonal(M, tol=1e-10):
        off = M[: sizes[0], sizes[0]:]
        return np.abs(off).max() < tol

    # forward: S block diagonal => S^{-1/2} block diagonal
    blocks = [rng.normal(size=(s, s)) for s in sizes]
    S_bd = np.zeros((N, N))
    i = 0
    for b, s in zip(blocks, sizes):
        S_bd[i:i + s, i:i + s] = b @ b.T + s * np.eye(s)
        i += s
    ev, U = np.linalg.eigh(S_bd)
    S_bd_invsqrt = U @ np.diag(ev**-0.5) @ U.T
    assert is_block_diagonal(S_bd_invsqrt)

    # contrapositive, at essentially PERFECT conditioning
    for eps in (1e-6, 1e-3, 0.3):
        S = np.eye(N)
        S[: sizes[0], sizes[0]:] = eps * rng.normal(size=(sizes[0], sizes[1]))
        S = 0.5 * (S + S.T) + np.eye(N) * 0.0
        assert not is_block_diagonal(S)
        cond = np.linalg.cond(S)
        ev, U = np.linalg.eigh(S)
        assert not is_block_diagonal(U @ np.diag(ev**-0.5) @ U.T), (
            f"S^{{-1/2}} came out block diagonal at cond = {cond}")
        for _ in range(25):
            X = np.zeros((N, N))
            i = 0
            for s in sizes:
                X[i:i + s, i:i + s] = rng.normal(size=(s, s))
                i += s
            if abs(np.linalg.det(X)) < 1e-8:
                continue
            assert not is_block_diagonal(X.T @ S @ X), (
                f"a block-diagonal congruence block-diagonalized S at cond = {cond}")
