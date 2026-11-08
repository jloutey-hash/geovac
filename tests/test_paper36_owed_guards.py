"""Owed guards for Paper 36 after the 2026-11-07 narrowing.

Paper: papers/group6_precision_observations/paper_36_bound_state_qed.tex.

Each test names the wrong answer it rejects ("Rejects:") and gets its
expected value from a route that does not share the code under test:
exact / high-precision closed forms, sympy, Wigner 3j sums, a literature
value, or an independent Fraction re-implementation of the LS-8a sum.

Loci backed:
  * abstract + sec:result + sec:residual: residuals -5.65 MHz (-0.534%)
    vs Lundeen-Pipkin and -5.64 MHz vs Bezginov 2019;
  * sec:bethe "Acceleration form": ln k0(2S) = 2.7858 at N = 40 (-0.92%)
    and 2.8955 at N = 50 (+2.98%), the N = 40 agreement a crossing;
  * sec:bethe "Velocity form": ln k0(2P) +2.40% at N = 40, +1.24% at N = 55;
  * sec:bethe "Definition and normalization": the (2l+1) bookkeeping and
    the s-state sum rule I_v(nS) = 2 Z^4 / n^3;
  * sec:twoloop, tab:ls8a_native: model totals and the 5 -> 6 slope 2.8
    under the corrected channel count (old count: 44.8 ... 3052.8, 3.3).
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent / "paper36_lamb_support"))
from bethe_log_sturmian import (  # noqa: E402
    bethe_log_acceleration,
    bethe_log_drake,
    build_velocity_spectrum,
)

# Drake & Swainson, Phys. Rev. A 41, 1243 (1990), Table I -- typed here
# from the paper's sec:result / sec:bethe, NOT imported from the support
# module, so a corrupted support-module reference cannot pass silently.
DRAKE_2S = 2.8117698931
DRAKE_2P = -0.0300167089

LAMB_LP81_MHZ = "1057.845"      # Lundeen & Pipkin, PRL 46, 232 (1981)
LAMB_BEZ19_MHZ = "1057.8298"    # Bezginov et al., Science 365, 1007 (2019)


def _pct(v: float, ref: float) -> float:
    """Signed relative error in percent, sign as printed in Paper 36
    (positive = larger magnitude than the reference)."""
    return 100.0 * (abs(v) - abs(ref)) / abs(ref)


# ===========================================================================
# 1. Residuals of the 1052.19 MHz assembly (two independent unit routes)
# ===========================================================================

def _lamb_assembly_mp(route: str) -> mp.mpf:
    """Paper 36 eqs. (eq:eides_se), (eq:eides_se_2p), (eq:uehling_vp),
    n = 2, Z = 1, m_r = m, at 40 digits with exact rational constants.

    route 'hartree': common factor alpha^3 Z^4 Ha, 1 Ha = 6 579 683 920.502 MHz.
    route 'mc2':     common factor alpha (Z alpha)^4 m c^2 from
                     m_e c^2 = 510 998.950 00 eV and 1 eV = 2.417 989 242e8 MHz.
    The two routes share no unit constant.
    """
    with mp.workdps(40):
        a = 1 / mp.mpf("137.035999084")
        Z, n = 1, 2
        if route == "hartree":
            unit = a ** 3 * Z ** 4 * mp.mpf("6579683920.502")
        else:
            unit = (a * (Z * a) ** 4 * mp.mpf("510998.95000")
                    * mp.mpf("2.417989242e8"))
        c = lambda fr: mp.mpf(fr.numerator) / fr.denominator  # noqa: E731
        se_2s = (c(Fraction(4, 3)) * (mp.log(1 / (Z * a) ** 2)
                                      - mp.mpf("2.8117698931"))
                 + c(Fraction(10, 9)))
        vp_2s = -c(Fraction(4, 15))
        se_2p = -c(Fraction(4, 3)) * mp.mpf("-0.0300167089") - c(Fraction(1, 6))
        return unit / (mp.pi * n ** 3) * (se_2s + vp_2s - se_2p)


def test_paper36_lamb_assembly_two_unit_routes_agree() -> None:
    """sec:result: 1052.19 MHz from two unit routes (Ha->MHz vs m c^2 in eV).

    Rejects: an assembly whose value depends on the unit route (a wrong
    Hartree/MHz constant or a dropped alpha power moves one route only),
    and the retired 38/45 constant (1025.06 MHz)."""
    h = _lamb_assembly_mp("hartree")
    e = _lamb_assembly_mp("mc2")
    assert abs(h - e) / h < 2e-8, (h, e)
    assert abs(h - mp.mpf("1052.19")) < mp.mpf("0.005"), h
    assert abs(h - mp.mpf("1025.06")) > 20


def test_paper36_residual_vs_lundeen_pipkin() -> None:
    """Abstract, sec:result table, conclusion: residual -5.65 MHz (-0.534%)
    against Lundeen--Pipkin 1057.845(9) MHz.

    Rejects: the retired residual -32.78 MHz (-3.10%) of the 38/45
    assembly, and a residual of the wrong sign."""
    L = _lamb_assembly_mp("hartree")
    r = L - mp.mpf(LAMB_LP81_MHZ)
    assert abs(r - mp.mpf("-5.65")) < mp.mpf("0.005"), r
    pct = 100 * r / mp.mpf(LAMB_LP81_MHZ)
    assert abs(pct - mp.mpf("-0.534")) < mp.mpf("0.0005"), pct


def test_paper36_residual_vs_bezginov_2019() -> None:
    """Abstract + sec:result: residual -5.64 MHz against Bezginov et al.
    2019, 1057.8298(32) MHz.

    Rejects: quoting the Lundeen--Pipkin residual (-5.65) against Bezginov;
    the two differ by 0.0152 MHz, which the 2-decimal print resolves."""
    L = _lamb_assembly_mp("hartree")
    r = L - mp.mpf(LAMB_BEZ19_MHZ)
    assert abs(r - mp.mpf("-5.64")) < mp.mpf("0.005"), r
    assert abs(r - mp.mpf("-5.65")) > mp.mpf("0.005"), r


# ===========================================================================
# 2. 2S acceleration-form Bethe log: not converged; N = 40 is a crossing
# ===========================================================================

def test_paper36_2s_acceleration_drift_predicts_crossing_fast() -> None:
    """sec:bethe acceleration form (fast sibling, N = 25, 30, 35, ~10 s).

    The error against Drake--Swainson is negative and still rising by more
    than 2.5 percentage points per five basis functions, and the linear
    continuation of the last step crosses zero between N = 40 and N = 50.

    Rejects: "converged at N = 40" -- a sequence converged near N = 40
    would already be moving by well under a point per five functions."""
    errs = [_pct(bethe_log_acceleration(2, 0, N=N, lam=0.5)["ln_k0"], DRAKE_2S)
            for N in (25, 30, 35)]
    assert all(e < 0 for e in errs), errs
    steps = [errs[1] - errs[0], errs[2] - errs[1]]
    assert all(s > 2.5 for s in steps), (errs, steps)
    slope = steps[-1] / 5.0
    n_zero = 35 + (-errs[-1]) / slope
    assert 40.0 < n_zero < 50.0, (errs, n_zero)


@pytest.mark.slow
def test_paper36_2s_acceleration_n40_n50_sign_change() -> None:
    """sec:bethe acceleration form: ln k0(2S) = 2.7858 at N = 40 (-0.92%)
    and 2.8955 at N = 50 (+2.98%); one cutoff past (N = 55) the drift
    continues upward (measured +3.48%).  ~50 s.

    Rejects: "converged at N = 40" (the error changes sign and keeps
    growing past the data the paper quotes)."""
    vals = {N: bethe_log_acceleration(2, 0, N=N, lam=0.5)["ln_k0"]
            for N in (40, 50, 55)}
    e = {N: _pct(v, DRAKE_2S) for N, v in vals.items()}
    assert abs(vals[40] - 2.7858) < 5e-4 and abs(e[40] - (-0.92)) < 0.01, e
    assert abs(vals[50] - 2.8955) < 5e-4 and abs(e[50] - 2.98) < 0.01, e
    assert e[40] < 0 < e[50]
    assert e[55] > e[50] + 0.2, e  # still drifting one cutoff past


# ===========================================================================
# 3. 2P velocity form: 2.4% at N = 40, 1.24% at N = 55
# ===========================================================================

def test_paper36_2p_velocity_monotone_from_above_fast() -> None:
    """sec:bethe velocity form (fast sibling, N = 20, 25, 30, ~6 s): the
    computed 2P logarithm approaches Drake--Swainson from the larger-
    magnitude side, its error falling monotonically (11.7, 7.2, 4.7 %).

    Rejects: a 2P value that undershoots the tabulated magnitude or has
    stalled (a stall above 2.4% would contradict the N = 40 and 55 values)."""
    errs = [_pct(bethe_log_drake(2, 1, N=N, lam=0.5)["ln_k0_drake"], DRAKE_2P)
            for N in (20, 25, 30)]
    assert all(e > 0 for e in errs), errs
    assert errs[0] > errs[1] > errs[2] > 2.4, errs
    assert abs(errs[2] - 4.75) < 0.05, errs


@pytest.mark.slow
def test_paper36_2p_velocity_n55() -> None:
    """sec:bethe velocity form: ln k0(2P) = -0.03039 at N = 55, +1.24% from
    -0.0300167 (~30 s).  The N = 40 value (+2.40%) is pinned in
    test_paper36_lamb_chain.py.

    Rejects: "2.4% is the converged 2P accuracy" and any 2P value at N = 55
    that is not between Drake and the N = 40 value."""
    v = bethe_log_drake(2, 1, N=55, lam=0.5)["ln_k0_drake"]
    assert abs(v - (-0.03039)) < 5e-6, v
    assert abs(_pct(v, DRAKE_2P) - 1.24) < 0.01, _pct(v, DRAKE_2P)


# ===========================================================================
# 4. The (2l+1) normalization and the s-state sum rule
# ===========================================================================

def _msum_dipole_angular(l: int, lp: int) -> sp.Expr:
    """sum over m, m', q of |<l m | C^1_q | l' m'>|^2 by brute Wigner-3j
    summation (Wigner-Eckart, reduced element <l||C^1||l'>)."""
    from sympy.physics.wigner import wigner_3j
    red2 = (2 * l + 1) * (2 * lp + 1) * wigner_3j(l, 1, lp, 0, 0, 0) ** 2
    tot = sp.Integer(0)
    for m in range(-l, l + 1):
        for mp_ in range(-lp, lp + 1):
            for q in (-1, 0, 1):
                tot += wigner_3j(l, 1, lp, -m, q, mp_) ** 2 * red2
    return sp.nsimplify(tot)


def test_paper36_angular_weight_is_double_m_sum() -> None:
    """sec:bethe normalization: the code's angular weight l_> is the sum
    over BOTH initial and final magnetic substates; per initial substate it
    is l_> / (2l + 1).  Exact Wigner-3j route.

    Rejects: reading l_> as the per-substate weight (which would make the
    uniform denominator 2 Z^4/n^3 correct without the 2l+1 factor)."""
    for l, lp in ((1, 0), (1, 2), (2, 1), (2, 3), (0, 1)):
        tot = _msum_dipole_angular(l, lp)
        assert tot == max(l, lp), (l, lp, tot)
        assert tot / (2 * l + 1) == sp.Rational(max(l, lp), 2 * l + 1)


def test_paper36_2p_needs_2l_plus_1_factor() -> None:
    """sec:bethe normalization: 'without it the computed 2P logarithm would
    be three times the tabulated one'.  N = 20 (fast); the (2l+1)-divided
    normalization lands near Drake--Swainson and the undivided one near
    three times it.

    Rejects: the denominator 2 Z^4/n^3 applied to the m-summed numerator
    (gives ~ -0.10, 3.3x the table at N = 20)."""
    res = bethe_log_drake(2, 1, N=20, lam=0.5)
    without = res["J_v"] / (2.0 / 2 ** 3)          # 2 Z^4 / n^3, no 2l+1
    with_f = res["ln_k0_drake"]
    assert abs(with_f / DRAKE_2P - 1.0) < 0.15, with_f
    assert abs(without / DRAKE_2P - 3.0) < 0.5, without


def test_paper36_s_state_sum_rule_closed_form() -> None:
    """sec:bethe normalization: I_v(nS) = (1/2)<nabla^2 V> = 2 Z^4 / n^3,
    from sympy's hydrogen radial function: (1/2) 4 pi Z |psi(0)|^2 =
    Z R_n0(0)^2 / 2.

    Rejects: the old 'Z^4/n^3 including spin factor 2' (1/8 at n = 2)."""
    from sympy.physics.hydrogen import R_nl
    r, Z = sp.symbols("r Z", positive=True)
    for n in (1, 2, 3, 4):
        R0 = sp.limit(R_nl(n, 0, r, Z), r, 0)
        I = sp.simplify(Z * R0 ** 2 / 2)
        assert sp.simplify(I - 2 * Z ** 4 / n ** 3) == 0, (n, I)
    R0 = sp.limit(R_nl(2, 0, r, 1), r, 0)
    assert sp.Rational(1, 1) * R0 ** 2 / 2 == sp.Rational(1, 4)


def test_paper36_s_state_sum_rule_pseudostates_fast() -> None:
    """sec:bethe normalization (numerical side, N = 10, 15, 20): the
    pseudostate sum I_v(2S) rises monotonically (nested Laguerre bases at
    fixed lambda) and already exceeds Z^4/n^3 = 1/8 at N = 15 while staying
    below 2 Z^4/n^3 = 1/4.

    Rejects: the limit 1/8 (a monotone-from-below sequence cannot pass its
    limit)."""
    I = [build_velocity_spectrum(2, 0, N, 0.5)["I_v"] for N in (10, 15, 20)]
    assert I[0] < I[1] < I[2] < 0.25, I
    assert I[1] > 0.125 * 1.2, I


@pytest.mark.slow
def test_paper36_s_state_sum_rule_n40() -> None:
    """The 0.2096 quoted with the sum rule (task record; not printed in the
    paper): I_v(2S) at N = 40, below 1/4 and above 1/8 (~11 s).

    Rejects: 1/8 as the limit; a value at N = 40 above the exact 1/4."""
    I40 = build_velocity_spectrum(2, 0, 40, 0.5)["I_v"]
    assert abs(I40 - 0.2096) < 5e-4, I40
    assert 0.125 < I40 < 0.25


# ===========================================================================
# 5. LS-8a model totals (independent exact Fraction route)
# ===========================================================================

def _w_independent(n1: int, n2: int, q: int, old: bool = False) -> int:
    """Independent channel count on 2j-integers.  Corrected rule: both
    spinors ((n+1)/2, n/2), vector components ((q+1)/2,(q-1)/2) and
    ((q-1)/2,(q+1)/2), SU(2) triangle AND integrality in L and R, inside
    the window |n1-n2| <= q <= n1+n2, n1+n2+q odd.  old=True reproduces
    the retired count (n2's content swapped, no integrality)."""
    if q < 1 or not (abs(n1 - n2) <= q <= n1 + n2) or (n1 + n2 + q) % 2 == 0:
        return 0
    a = (n1 + 1, n1)
    b = (n2, n2 + 1) if old else (n2 + 1, n2)

    def tri(x: int, y: int, z: int) -> bool:
        if not abs(x - y) <= z <= x + y:
            return False
        return True if old else (x + y + z) % 2 == 0

    return sum(1 for gL, gR in ((q + 1, q - 1), (q - 1, q + 1))
               if tri(a[0], gL, b[0]) and tri(a[1], gR, b[1]))


def _ls8a_fraction(n_max: int, old: bool = False, n_ext: int = 1):
    """Rainbow and crossed sums of Paper 36 sec:twoloop in exact Fractions:
    lambda_n^2 = (n+3/2)^2, g_n = 2(n+1)(n+2), d_q / mu_q = 1."""
    lam2 = lambda n: Fraction(2 * n + 3, 2) ** 2  # noqa: E731
    g = lambda n: 2 * (n + 1) * (n + 2)          # noqa: E731
    w = lambda x, y, q: _w_independent(x, y, q, old)  # noqa: E731
    rain = cr = Fraction(0)
    rng = range(n_max + 1)
    qs = range(1, 2 * n_max + 2)
    for n1 in rng:
        for n2 in rng:
            for n3 in rng:
                pre = Fraction(g(n1) * g(n2) * g(n3)) / (lam2(n1) * lam2(n2) * lam2(n3))
                for qa in qs:
                    for qb in qs:
                        wr = w(n_ext, n1, qa) * w(n3, n_ext, qa) * w(n1, n2, qb) * w(n2, n3, qb)
                        wc = w(n_ext, n1, qa) * w(n2, n3, qa) * w(n1, n2, qb) * w(n3, n_ext, qb)
                        # mu_q = q(q+2) and d_q = q(q+2) cancel per photon line
                        rain += pre * wr
                        cr += pre * wc
    return rain, cr


# tab:ls8a_native, printed (rainbow, crossed, total, C)
_TAB = {2: (587.4, 467.0, 1054.4, 6.68), 3: (1680.9, 1313.8, 2994.8, 18.96),
        4: (3653.0, 2786.3, 6439.3, 40.78), 5: (6761.9, 5140.6, 11902.5, 75.37),
        6: (11265.7, 8507.5, 19773.2, 125.22)}
_OLD_TOTALS = {2: 44.8, 3: 266.3, 4: 761.3, 5: 1657.8, 6: 3052.8}


def test_paper36_ls8a_table_exact_route() -> None:
    """tab:ls8a_native: every printed cell, plus totals 1054.36 (n_max = 2)
    and 19773.16 (n_max = 6), from the independent exact Fraction sum;
    C = total / (4 pi)^2.

    Rejects: the pre-2026-11-07 table (totals 44.8 ... 3052.8)."""
    for nm, (r_p, c_p, t_p, C_p) in _TAB.items():
        r, c = _ls8a_fraction(nm)
        tot = float(r + c)
        assert abs(float(r) - r_p) < 0.05 and abs(float(c) - c_p) < 0.05, nm
        assert abs(tot - t_p) < 0.05, (nm, tot)
        assert abs(tot / (16 * math.pi ** 2) - C_p) < 0.005, nm
        assert abs(tot - _OLD_TOTALS[nm]) > 40, nm
    assert abs(float(sum(_ls8a_fraction(2))) - 1054.36) < 0.005
    assert abs(float(sum(_ls8a_fraction(6))) - 19773.16) < 0.005


def test_paper36_ls8a_slope_and_old_count() -> None:
    """sec:twoloop: local log-log slope of the total between n_max = 5 and 6
    is 2.8 (it was 3.3, i.e. 3.349, under the old count, whose totals were
    44.8, 266.3, 761.3, 1657.8, 3052.8 -- reproduced here by the old rule,
    which shows the 'old' branch is the count the paper describes).

    Rejects: slope 3.3 / the old totals as the current table, and 3.4 as the
    old-count slope (the previously printed rounding)."""
    t5, t6 = (float(sum(_ls8a_fraction(n))) for n in (5, 6))
    assert abs(math.log(t6 / t5) / math.log(6 / 5) - 2.8) < 0.05
    old = {n: float(sum(_ls8a_fraction(n, old=True))) for n in _OLD_TOTALS}
    for n, v in _OLD_TOTALS.items():
        assert abs(old[n] - v) < 0.05, (n, old[n])
    # The old-count slope is ln(3052.8/1657.8)/ln 1.2 = 3.349 by this route;
    # the paper prints "3.3".  Pinned at 3.35 +- 0.005 (rejects 3.4 and 2.8).
    assert abs(math.log(old[6] / old[5]) / math.log(6 / 5) - 3.35) < 0.005


def test_paper36_ls8a_production_matches_exact_route() -> None:
    """sec:twoloop: geovac.two_loop_self_energy (mpmath, channel count from
    geovac.qed_vertex.so4_channel_count) equals the independent Fraction
    route at n_max = 2..6 to 1e-12 relative, and its window W lies in
    {0, 2} as the paper states.

    Rejects: the old channel count in geovac/qed_vertex.py (n2's L/R
    content swapped and no integrality), which gives 44.8 at n_max = 2."""
    from geovac.two_loop_self_energy import iterated_se_dimensionless
    from geovac.qed_vertex import so4_channel_count
    for nm in range(2, 7):
        r, c = _ls8a_fraction(nm)
        prod = iterated_se_dimensionless(n_ext=1, n_max=nm)
        assert abs(float(prod["rainbow"]) / float(r) - 1) < 1e-12, nm
        assert abs(float(prod["crossed"]) / float(c) - 1) < 1e-12, nm
    vals = {so4_channel_count(a, b, q) for a in range(7) for b in range(7)
            for q in range(1, 14)}
    assert vals <= {0, 2} and 2 in vals, vals
