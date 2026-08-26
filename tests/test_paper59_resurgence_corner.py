"""Backing tests for Paper 59 sec:modular reconciliation (v4.100.0).

Two [MEASURED] facts the reconciled sec:modular asserts, both self-contained
(no debug/ import, per the permanent-record policy):

  (A) Corner leading orders of the fibre J(s,t) [backs the outer-wall reframe]:
      only the (0,0) corner is fractionally singular (rho^{3/2}); the three
      oscillatory corners (b=s+t != 0) are integer-leading (rho^2).  Hence no
      fractional-power "corner subtraction" closes the outer integral there --
      the >=40-digit route is the modular/q-series representation, not more (s,t)
      quadrature.

  (B) Fibre resurgence [backs "T2 is a resurgent Lambert series, the finite
      Gamma(2) MMV is its D->0 shadow"]:  the one-mass fibre N(D)'s large-D
      asymptotic series is Gevrey-1 (factorially divergent) with Borel radius 2 --
      i.e. a period of the IRREGULAR rank-4 connection eq:pf, not a Fuchsian
      modular period.  (The structural fact that eq:pf is irregular is the
      already-tracked test_L4_irregular_at_infinity; this is the resurgence
      signature of the physical master period itself.)
"""
import mpmath as mp
import pytest


# --- fibre J(s,t) = int_0^inf j0(k(s+t)) P(s,k) P(t,k) dk (Paper 59 collinear) ---
def _P(x, k):
    c = x * (1 - x)
    D = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-D) * (D ** -3 + 3 * D ** -4 + 3 * D ** -5)


def _j0(z):
    return mp.mpf(1) if z == 0 else mp.sin(z) / z


def _J(s, t):
    b = s + t
    f = lambda k: _j0(b * k) * _P(s, k) * _P(t, k)
    if b < mp.mpf('0.2'):                       # (0,0) corner: non-oscillatory
        return mp.quad(f, [0, 1, 2, 4, 8, 16, mp.inf])
    return mp.quadosc(f, [0, mp.inf], period=2 * mp.pi / b)   # oscillatory corners


def _corner_exponent(which, rhos):
    """Local exponent alpha of J ~ rho^alpha along a ray into the named corner,
    measured as log2(J(rho)/J(rho/2)) for the finest pair."""
    Js = []
    for r in rhos:
        s = t = r if which == '00' else 1 - r
        Js.append(_J(s, t))
    return mp.log(Js[-2] / Js[-1]) / mp.log(2)   # finest halving = closest to the limit


@pytest.mark.slow
def test_paper59_corner_orders_only_00_is_fractional():
    """(A) (0,0) is rho^{3/2} (fractional); (1,1) is rho^2 (integer) -- distinct
    classes, so only (0,0) needs the sigma^2 Duffy; the oscillatory corners carry
    no fractional-power subtraction."""
    mp.mp.dps = 30
    rhos = [mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005'), mp.mpf('0.0025')]
    a00 = _corner_exponent('00', rhos)
    a11 = _corner_exponent('11', rhos)
    # (0,0): fractional, converging to 3/2 from below (subleading rho^{5/2})
    assert mp.mpf('1.40') < a00 < mp.mpf('1.62'), a00
    # (1,1): integer leading, ~2
    assert mp.mpf('1.90') < a11 < mp.mpf('2.05'), a11
    # the two corners are genuinely different classes (not both ~integer, not both ~3/2)
    assert a11 - a00 > mp.mpf('0.35'), (a00, a11)


def _watson_coeffs(rho, kmax):
    """Large-D asymptotic N(D) ~ e^{-D} sum_k b_k D^{-(k+1/2)} via Watson's lemma
    at the branch point x=1 (x=1+u; 1/sqrt(Q)=u^{-1/2} g(u), b_k = g_k Gamma(k+1/2))."""
    g = mp.taylor(lambda u: 1 / mp.sqrt((2 + u) * (rho * (1 + u) ** 2 + 1 - rho)), 0, kmax + 2)
    return [g[k] * mp.gamma(k + mp.mpf('0.5')) for k in range(len(g))]


def _N_exact(D, rho):
    f = lambda x: mp.e ** (-D * x) / mp.sqrt((x * x - 1) * (rho * x * x + 1 - rho))
    return mp.quad(f, [1, mp.inf])


def test_paper59_fibre_resurgence_gevrey1_borel2():
    """(B) N(D)'s large-D series is Gevrey-1 (optimal truncation k* exists, error
    rises past it, k* grows with D) with Borel radius 2 (|b_{k+1}/b_k|/(k+1/2)->1/2).
    => the fibre is an irregular (resurgent) period, not a Fuchsian modular period."""
    mp.mp.dps = 60
    rho = mp.mpf(1) / 5
    b = _watson_coeffs(rho, 40)

    kstars = {}
    for D in (mp.mpf(8), mp.mpf(12)):
        Nex = _N_exact(D, rho)
        partial = mp.mpf(0)
        errs = []
        for k in range(len(b)):
            partial += b[k] * D ** (-(k + mp.mpf('0.5')))
            errs.append(abs(mp.e ** (-D) * partial - Nex))
        kstar = min(range(len(errs)), key=lambda k: errs[k])
        kstars[D] = kstar
        # divergent (asymptotic) signature: error turns around at k* then RISES
        assert errs[kstar + 3] > errs[kstar], (D, kstar)
        assert errs[kstar + 6] > errs[kstar + 3], (D, kstar)
        # and the optimal truncation actually resolves the value (not a trivial k*=0)
        assert kstar >= 5, (D, kstar)

    # optimal truncation grows with D (k* ~ S*D/2) -- the hallmark of a fixed Borel radius
    assert kstars[mp.mpf(12)] > kstars[mp.mpf(8)], kstars

    # Borel radius S=2: |b_{k+1}/b_k| / (k+1/2) -> 1/S = 1/2
    ratio = abs(b[36] / b[35]) / (35 + mp.mpf('0.5'))
    assert abs(ratio - mp.mpf('0.5')) < mp.mpf('0.05'), ratio
