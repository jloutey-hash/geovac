"""Paper 59 -- the (KW) factorization of the collinear observable T2 and its
high-precision value (CHANGELOG v4.104.0; drivers debug/beta2_t2_*.py).

The exact reduction: j0(z) = int_0^1 cos(zw) dw decouples the b = s+t phase and
P(s,k) = P(1-s,k) makes the inner integral real, so the (s,t) double integral
becomes the SQUARE of a one-dimensional integral:

    T2 = (8/pi) int_0^inf dk int_0^1 dw cos(kw) R(k,w)^2,
    R(k,w) = int_0^1 cos(kw(s - 1/2)) P(s,k) ds,
    P(x,k) = c e^{-D} (D^-3 + 3 D^-4 + 3 D^-5),  c = x(1-x),  D = sqrt(c k^2 + 1).

In this frame the complex off-axis (s,t) singularity that capped the old outer
quadrature at ~13-14 digits does not exist; the truncation tail is the whole
precision question and falls as K^-7.  The 66-digit certified value (six
parameter-disjoint runs, two independent parallel configurations to 1.7e-67,
identity verified against an independent 2D evaluation to 69-96 digits) is

    T2 = 0.395355765901713964325229296804847564260563977867082108935234265469

which corrects the previously frozen ~19-digit anchor in its 19th digit
(0.3953557659017139641 -> ...0171396432), consistent with the v4.97.0 finding
that the honest cross-validated ceiling of the old frame was 15-16 digits.
"""
import numpy as np
import pytest

T2_66 = "0.395355765901713964325229296804847564260563977867082108935234265469"
T2_OLD_ANCHOR = "0.3953557659017139641"


def _kw_value(K: float, n_k: int = 4000, n_w: int = 120, n_s: int = 240) -> float:
    """Plain truncated (KW) evaluation on [0, K] (float64, Gauss-Legendre)."""
    xk, wk = np.polynomial.legendre.leggauss(n_k)
    k = 0.5 * K * (xk + 1.0)
    wk = 0.5 * K * wk
    xw, ww = np.polynomial.legendre.leggauss(n_w)
    w = 0.5 * (xw + 1.0)
    ww = 0.5 * ww
    xs, ws = np.polynomial.legendre.leggauss(n_s)
    s = 0.5 * (xs + 1.0)
    ws = 0.5 * ws
    c = s * (1.0 - s)
    total = 0.0
    for ki, wki in zip(k, wk):
        D = np.sqrt(c * ki * ki + 1.0)
        P = c * np.exp(-D) * (D ** -3 + 3.0 * D ** -4 + 3.0 * D ** -5)
        # R(k, w) for all w at once: cos(k w (s-1/2)) @ (ws * P)
        phase = np.cos(np.outer(ki * w, s - 0.5))
        R = phase @ (ws * P)
        total += wki * np.dot(ww * np.cos(ki * w), R * R)
    return (8.0 / np.pi) * total


def test_kw_identity_converges_to_value():
    """Truncated (KW) approaches the certified value; the K^-7 tail dominates
    until the float64 roundoff floor (~1e-12 at this resolution)."""
    target = float(T2_66)
    e20 = abs(_kw_value(20.0, n_k=2000) - target)
    e40 = abs(_kw_value(40.0) - target)
    e80 = abs(_kw_value(80.0, n_k=8000) - target)
    assert e20 < 2e-9           # tail ~6e-10 at K=20
    assert e40 < 1e-11          # tail ~5e-12 at K=40
    assert e80 < 1e-11          # float64 floor; certified digits live in the
    assert e20 > 50.0 * e40     # mpmath drivers (K^-7: 2^7 = 128 nominal)


def test_kw_matches_independent_2d_at_fixed_k():
    """The factorization identity at fixed k: the (s,t) double integral of
    cos(kw(s+t)) P(s)P(t) equals R(k,w)^2 (reality via s -> 1-s symmetry)."""
    rng_k, wv = 3.7, 0.6
    n = 400
    xs, ws = np.polynomial.legendre.leggauss(n)
    s = 0.5 * (xs + 1.0)
    wgt = 0.5 * ws
    c = s * (1.0 - s)
    D = np.sqrt(c * rng_k ** 2 + 1.0)
    P = c * np.exp(-D) * (D ** -3 + 3.0 * D ** -4 + 3.0 * D ** -5)
    # 2D side: int int cos(k w (s+t)) P(s) P(t) ds dt
    ck = np.cos(rng_k * wv * (s[:, None] + s[None, :]))
    lhs = np.einsum("i,j,ij->", wgt * P, wgt * P, ck)
    # 1D side: Re[Q^2] with Q = e^{i k w/2} R  =>  cos(kw) R^2
    R = np.dot(wgt * P, np.cos(rng_k * wv * (s - 0.5)))
    assert abs(lhs - np.cos(rng_k * wv) * R * R) < 1e-14


def test_old_anchor_corrected_in_digit_19():
    """Regression lock of the two recorded strings (18 shared digits, the 19th
    differing); the EXECUTABLE witness of the correction is
    test_kw_mpmath_witness below."""
    assert T2_66[:20] == T2_OLD_ANCHOR[:20]          # first 18 significant digits
    assert T2_66[:21] != T2_OLD_ANCHOR[:21]          # the 19th digit differs
    diff = abs(float(T2_66[:25]) - float(T2_OLD_ANCHOR))
    assert diff < 5e-19                               # a last-digit correction, not a discrepancy


@pytest.mark.slow
def test_kw_mpmath_witness():
    """Executable witness from the tracked tree: geovac.t2_kw evaluates the (KW)
    representation with the analytic Watson tail at K=60 and reproduces the
    certified value to ~21 digits -- independently confirming the digit-19
    correction of the old anchor (no reliance on stored literals)."""
    import mpmath as mp
    from geovac.t2_kw import F_buckets, TAIL, int_F

    mp.mp.dps = 30
    main = int_F(mp.mpf(0), mp.mpf(60), 15, 48, 220, p=6)
    tail = TAIL(mp.mpf(60), F_buckets(60))
    val = (8 / mp.pi) * (main + tail)
    assert abs(val - mp.mpf(T2_66)) < mp.mpf("1e-20")
    assert abs(val - mp.mpf(T2_OLD_ANCHOR)) > mp.mpf("1e-19")   # old anchor excluded
