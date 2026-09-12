"""Backing for Paper 60 sec:resource's direct-encoding construction (2026-09-12).

The open item after v5.11.1 was whether G admits a DIRECT block-encoding rather
than a composed one.  Composing P^-1/2 with (I-C) makes the subnormalization
inherit ||P^-1/2||^2 ~ n^2, so the total scales as n^2 where eq:amplitude_floor
allows n -- a factor ||G||-vs-alpha of 8.6e3 wasted at n = 160 on operation
order alone.

The construction.  G's symbol is the RATIO of the two symbols, and it is
bounded: numerator and denominator both vanish quadratically at chi = pi so the
quotient tends to (kR)^2/24, and at chi -> 0 it tends to 1/4.  In the
s = kR cot(chi/2) variable,

    ratio(s) = (1 - sinc(s)) (s^2 + (kR)^2) / (4 s^2),

so the Toeplitz-minus-Hankel matrix B built from the ratio's OWN cosine
coefficients is directly constructible, and a circulant-embedded encoding of it
carries alpha = ||ratio||_inf = O(1).

B is not G -- it differs by the finite-section commutator, ~12% in operator
norm -- and the point of this test is that the difference does not matter:
whitening with X = P^-1/2 B^-1/2 gives X^T (I-C) X = B^-1/2 G B^-1/2, whose
conditioning is bounded near 1.23, and ||X|| lands on the invariant floor.

Self-contained: the coefficient quadrature is inlined rather than imported from
the prunable debug/ tree (C22 check D).
"""
import numpy as np
import pytest

from geovac.sturmian_sigma_law import sw_cross_block

GL_N = 48
_xg, _wg = np.polynomial.legendre.leggauss(GL_N)
KR = 2.0
M_QUAD = 200_001


def _ratio_less_quarter(s, kR):
    with np.errstate(divide="ignore", invalid="ignore"):
        v = (kR**2 - np.sinc(s / np.pi) * (s**2 + kR**2)) / (4 * s**2)
    return np.where(s < 1e-8, (kR**2 / 6.0 - 1.0) / 4.0, v)


def ratio_coeff(j, kR, Smax=3000.0):
    """(1/pi) int_0^pi cos(j chi) ratio(chi) dchi.

    The s -> infinity limit 1/4 is subtracted (its cosine coefficients vanish for
    j >= 1, so this is exact) to make the integrand decay like 1/s^3; it is added
    back into r_0 analytically.
    """
    edges, s = [0.0], 0.0
    while s < Smax:
        s += min(np.pi / 2, np.pi / max(2 * j * kR / (kR * kR + s * s), 1e-300))
        edges.append(min(s, Smax))
    e = np.asarray(edges)
    mid, half = 0.5 * (e[:-1] + e[1:]), 0.5 * (e[1:] - e[:-1])
    sv = (mid[:, None] + half[:, None] * _xg[None, :]).ravel()
    wv = (half[:, None] * _wg[None, :]).ravel()
    integ = (np.cos(2 * j * np.arctan2(kR, sv)) * _ratio_less_quarter(sv, kR)
             * 2 * kR / (kR * kR + sv * sv))
    val = float(np.dot(wv, integ) / np.pi)
    return val + 0.25 if j == 0 else val


def tmh(coeffs, n):
    g = lambda j: coeffs.get(j, 0.0)
    return np.array([[g(abs(a - b)) - g(a + b) for b in range(1, n + 1)]
                     for a in range(1, n + 1)])


def tridiag(n):
    return (np.diag(2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1)
            + np.diag(np.ones(n - 1), -1))


def inv_sqrt(M):
    ev, U = np.linalg.eigh(M)
    assert ev.min() > 0, f"not PD, lam_min={ev.min():.3e}"
    return U @ np.diag(ev ** -0.5) @ U.T


def _build(n):
    coeffs = {j: ratio_coeff(j, KR) for j in range(0, 2 * n + 2)}
    B = tmh(coeffs, n)
    A = np.eye(n) - sw_cross_block(KR, n, M=M_QUAD)
    P_is = inv_sqrt(tridiag(n))
    return A, B, P_is, P_is @ A @ P_is


def test_ratio_symbol_sup_equals_the_norm_of_G():
    """||ratio||_inf = ||G||, so a direct encoding carries no excess alpha.

    WRONG ANSWER REJECTED: that a direct encoding would still pay the composed
    penalty.  The composed alpha is ||P^-1/2||^2 ||A|| ~ 3.2e3 at n=160; the sup
    of the ratio symbol is ~0.37.  The test requires the sup to match ||G|| to
    1% AND to be smaller than the composed alpha by more than three orders of
    magnitude, so a value anywhere near the composed penalty fails.
    """
    chi = np.linspace(1e-9, np.pi - 1e-12, 200_000)
    sup = float(np.max(np.abs(_ratio_less_quarter(KR / np.tan(chi / 2), KR) + 0.25)))

    n = 80
    A, B, P_is, G = _build(n)
    assert abs(sup / np.linalg.norm(G, 2) - 1.0) < 0.01, (
        f"sup(ratio)={sup:.4f} vs ||G||={np.linalg.norm(G,2):.4f}")
    composed = np.linalg.norm(P_is, 2) ** 2 * np.linalg.norm(A, 2)
    assert composed / sup > 1e3, f"composed penalty {composed:.1f} not >> sup {sup:.4f}"


def test_direct_object_whitens_to_bounded_residual_conditioning():
    """The decisive claim: cond(B^-1/2 G B^-1/2) is BOUNDED, near 1.23.

    B differs from G by ~12% in operator norm, and the test asserts BOTH that
    the discrepancy is real and that it is harmless -- otherwise a reader could
    conclude either that B equals G (it does not) or that a 12% error disqualifies
    it (it does not).

    WRONG ANSWER REJECTED: that the finite-section discrepancy grows with n, or
    that it degrades the conditioning.  The residual conditioning must stay under
    1.3 at every n and its increments must shrink.
    """
    resids, diffs = [], []
    for n in (20, 40, 80):
        A, B, P_is, G = _build(n)
        Bis = inv_sqrt(B)
        resids.append(np.linalg.cond(Bis @ G @ Bis))
        diffs.append(np.linalg.norm(G - B, 2) / np.linalg.norm(G, 2))

    assert all(0.08 < d < 0.16 for d in diffs), f"discrepancy not ~12%: {diffs}"
    assert max(diffs) - min(diffs) < 0.02, f"discrepancy grows with n: {diffs}"
    assert all(r < 1.3 for r in resids), f"residual conditioning unbounded: {resids}"
    incs = [b - a for a, b in zip(resids, resids[1:])]
    assert all(abs(b) < abs(a) for a, b in zip(incs, incs[1:])), f"not settling: {resids}"


def test_direct_whitening_attains_the_amplitude_floor():
    """||P^-1/2 B^-1/2|| -> ||A^-1/2||: the floor of eq:amplitude_floor, reached.

    This is what turns the O(n^2) composed total into O(n).

    WRONG ANSWER REJECTED: that only the exact G^-1/2 reaches the floor and the
    directly-built B^-1/2 does not.  The ratio must approach 1 from above AND
    tighten with n -- a fixed offset, however small, would fail the second.
    """
    ratios = []
    for n in (20, 40, 80):
        A, B, P_is, _ = _build(n)
        X = P_is @ inv_sqrt(B)
        ratios.append(np.linalg.norm(X, 2) / np.linalg.norm(inv_sqrt(A), 2))

    assert all(1.0 <= r < 1.02 for r in ratios), f"floor not attained: {ratios}"
    for a, b in zip(ratios, ratios[1:]):
        assert b - 1.0 < (a - 1.0), f"not tightening toward the floor: {ratios}"


@pytest.mark.slow
def test_direct_encoding_still_holds_one_cutoff_further():
    """n = 160, outside the window the paper tabulates (guard-asymptotics rule)."""
    A, B, P_is, G = _build(160)
    Bis = inv_sqrt(B)
    assert np.linalg.cond(Bis @ G @ Bis) < 1.3
    assert np.linalg.norm(P_is @ Bis, 2) / np.linalg.norm(inv_sqrt(A), 2) < 1.01
