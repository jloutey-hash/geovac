"""Paper 1 abstract claim (iii): Berry phase = 0; log-holonomy = -2 ln((n+1)/n).

WHY THIS FILE EXISTS (2026-09-02, trunk-FULL remediation, carryforward B6).
The 2026-09-01 /qa trunk run found Paper 1's third abstract headline had ZERO
live backing: the only artifact was tests/_archive/dead_ends/test_berry_phase.py
-- not collected by default, tied to the RETRACTED k = 2.113 claim, and testing
a different quantity (summed multi-plaquette totals).  Section 13.4a requires a
computational verification for every paper equation; this file supplies it for
eq:berry_arg / eq:berry_zero / eq:log_holonomy / eq:log_holonomy_exact.

It rebuilds the quantities FROM THE PAPER'S OWN PRINTED FORMULAS -- the
Biedenharn-Louck matrix elements and the topological edge weights w = 1/(n1 n2)
-- not from any solver code, so it verifies the paper's claim, not another
test's implementation.

NON-TAUTOLOGY CONTROLS.  theta = 0 is exactly the "suspiciously clean" shape
QA doctrine warns about, so each zero carries a control proving the machinery
CAN detect a nonzero: a synthetically phased plaquette must yield arg != 0,
and a synthetically warped weight function must break the closed form.
"""
from __future__ import annotations

import math
import cmath

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# The paper's printed matrix elements (eq:berry_arg context):
#   <l, m+1 | L+ | l, m>   = sqrt((l - m)(l + m + 1))          (SU(2))
#   <n+1, l | T+ | n, l>   = sqrt((n - l)(n + l + 1) / 4)      (SU(1,1))
# Valid plaquette at (n, l, m) needs: l <= n - 1 is NOT required by the
# operators themselves -- T+ needs n > l (else the element vanishes and the
# plaquette degenerates), L+ needs m < l.
# ---------------------------------------------------------------------------

def t_plus(n: int, l: int) -> float:
    return math.sqrt((n - l) * (n + l + 1) / 4.0)


def l_plus(l: int, m: int) -> float:
    return math.sqrt((l - m) * (l + m + 1))


def _states(n_max: int):
    return [(n, l, m) for n in range(1, n_max + 1) for l in range(n) for m in range(-l, l + 1)]


def build_operators(n_max: int, phase=None):
    """T+ and L+ as matrices from the paper's printed elements.  `phase`, if
    given, is a function (n, l, m) -> complex unit multiplying the L+ element
    that LEAVES state (n, l, m) -- a state-dependent U(1) dressing (Paper 1
    Sec. V.C item 2), used by the non-tautology control below."""
    st = _states(n_max); idx = {s: i for i, s in enumerate(st)}; N = len(st)
    T = np.zeros((N, N), dtype=complex); L = np.zeros((N, N), dtype=complex)
    for (n, l, m) in st:
        i = idx[(n, l, m)]
        if (n + 1, l, m) in idx:
            T[idx[(n + 1, l, m)], i] = t_plus(n, l)
        if m < l:
            ph = 1.0 if phase is None else phase(n, l, m)
            L[idx[(n, l, m + 1)], i] = l_plus(l, m) * ph
    return st, idx, T, L


def plaquette_holonomy(n: int, l: int, m: int, ops=None) -> complex:
    """The eq:berry_arg product around |n,l,m> -> T+ -> L+ -> T- -> L-, read
    off the OPERATOR MATRICES with true adjoints (T- = T+^dagger,
    L- = L+^dagger).  Rebuilt 2026-09-03: the earlier version re-used the
    same two scalars for the return legs, making the product a perfect
    square for ANY phase assignment -- a guard that could not fail."""
    st, idx, T, L = ops if ops is not None else build_operators(max(n + 1, 3))
    i0 = idx[(n, l, m)]; i1 = idx[(n + 1, l, m)]; i2 = idx[(n + 1, l, m + 1)]; i3 = idx[(n, l, m + 1)]
    Tm = T.conj().T; Lm = L.conj().T
    return complex(Lm[i0, i3] * Tm[i3, i2] * L[i2, i1] * T[i1, i0])


def valid_plaquettes(n_max: int):
    for n in range(1, n_max):          # T+ must land inside the lattice
        for l in range(0, n):          # n > l => t_plus nonzero
            for m in range(-l, l):     # m < l => l_plus nonzero
                yield n, l, m


# ---------------------------------------------------------------------------
# eq:berry_zero -- theta = 0 for EVERY plaquette
# ---------------------------------------------------------------------------

def test_berry_phase_vanishes_on_every_plaquette():
    checked = 0
    ops = build_operators(13)
    for n, l, m in valid_plaquettes(12):
        h = plaquette_holonomy(n, l, m, ops)
        assert h.real > 0.0, (n, l, m, h)          # nondegenerate
        assert abs(cmath.phase(h)) == 0.0, (n, l, m, h)
        checked += 1
    assert checked > 200      # 12-shell lattice: many plaquettes, not a corner


def test_berry_phase_detector_actually_detects():
    """NON-TAUTOLOGY CONTROL (rebuilt 2026-09-03): dress the L+ elements with
    a STATE-DEPENDENT U(1) phase exp(i alpha n) (Paper 1 Sec. V.C item 2).
    Around the plaquette the two L legs sit at shells n+1 and n, so the
    holonomy picks up exp(i alpha) -- a genuine nonzero Berry phase that the
    matrix-built detector must see.  A sign-only convention (Condon-Shortley
    (-1)^m) cancels on the adjoint leg and gives 0, as the paper says."""
    alpha = 0.7
    ops = build_operators(5, phase=lambda n, l, m: cmath.exp(1j * alpha * n))
    h = plaquette_holonomy(3, 1, 0, ops)
    assert abs(cmath.phase(h) - alpha) < 1e-12
    ops_cs = build_operators(5, phase=lambda n, l, m: (-1.0) ** m)
    assert abs(cmath.phase(plaquette_holonomy(3, 1, 0, ops_cs))) < 1e-12


def test_berry_zero_traces_to_realness_not_magnitude():
    """The paper's mechanism claim: the phase vanishes BECAUSE the four
    Biedenharn-Louck factors are real positive -- not because they are
    small, unit, or mutually cancelling.  Verify the factors genuinely vary
    (no hidden normalisation) while the phase stays pinned at zero."""
    mags = set()
    ops = build_operators(9)
    for n, l, m in valid_plaquettes(8):
        h = plaquette_holonomy(n, l, m, ops)
        mags.add(round(abs(h), 9))
        assert cmath.phase(h) == 0.0
    assert len(mags) > 20     # magnitudes spread; only the phase is rigid


# ---------------------------------------------------------------------------
# eq:log_holonomy / eq:log_holonomy_exact -- Theta(n) = -2 ln((n+1)/n)
# ---------------------------------------------------------------------------

def edge_weight(n1: int, n2: int) -> float:
    """The paper's topological weight w = 1/(n1 n2) on an edge between
    shells n1 and n2 (azimuthal edges have n1 = n2)."""
    return 1.0 / (n1 * n2)


def log_holonomy(n: int) -> float:
    """eq:log_holonomy around the standard plaquette based at shell n:
    ln w12 + ln w23 - ln w34 - ln w41 with
      w12 = T+ edge (n, n+1),  w23 = L+ edge at shell n+1,
      w34 = T- edge (n+1, n),  w41 = L- edge at shell n."""
    w12 = edge_weight(n, n + 1)
    w23 = edge_weight(n + 1, n + 1)
    w34 = edge_weight(n + 1, n)
    w41 = edge_weight(n, n)
    return (math.log(w12) + math.log(w23)
            - math.log(w34) - math.log(w41))


def test_log_holonomy_closed_form_exact():
    """Theta(n) = -2 ln((n+1)/n), exactly, from the edge weights."""
    for n in range(1, 60):
        assert abs(log_holonomy(n) - (-2.0 * math.log((n + 1) / n))) < 1e-13


def test_log_holonomy_sign_and_asymptote():
    """Theta(n) is NEGATIVE and Theta(n) ~ -2/n: n * Theta(n) -> -2.

    (The v1 paper printed the asymptote with a dropped sign; the quantity
    -2 ln((n+1)/n) is manifestly negative.)"""
    for n in (1, 5, 20, 100):
        assert log_holonomy(n) < 0.0
    products = [n * log_holonomy(n) for n in (100, 400, 1600)]
    for p in products:
        assert p < 0.0
    # residual decays as 1/n (next Taylor order): at n=1600 it is ~6.3e-4
    assert abs(products[-1] - (-2.0)) < 1e-3
    # monotone approach to -2 from below in magnitude terms
    assert abs(products[0] + 2.0) > abs(products[-1] + 2.0)


def test_log_holonomy_decay_exponent_is_one():
    """The n^{-1} decay claim (k = 1.0): log-log slope of |Theta| vs n."""
    ns = np.array([10, 20, 40, 80, 160, 320], dtype=float)
    th = np.array([abs(log_holonomy(int(n))) for n in ns])
    slope = np.polyfit(np.log(ns), np.log(th), 1)[0]
    assert abs(slope - (-1.0)) < 2e-2, slope


def test_log_holonomy_detector_actually_detects():
    """NON-TAUTOLOGY CONTROL: a warped weight function (w = 1/(n1 n2)^2)
    must BREAK the -2 ln((n+1)/n) closed form -- proving the closed-form
    assertion is sensitive to the weight law, not an algebraic identity of
    any four logs."""
    def warped(n1, n2):
        return 1.0 / (n1 * n2) ** 2

    n = 7
    w12, w23 = warped(n, n + 1), warped(n + 1, n + 1)
    w34, w41 = warped(n + 1, n), warped(n, n)
    theta_warped = (math.log(w12) + math.log(w23)
                    - math.log(w34) - math.log(w41))
    assert abs(theta_warped - (-2.0 * math.log((n + 1) / n))) > 0.05
    # (it lands on -4 ln((n+1)/n): the exponent doubles the circulation)
    assert abs(theta_warped - (-4.0 * math.log((n + 1) / n))) < 1e-13
