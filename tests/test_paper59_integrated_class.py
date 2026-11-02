"""Paper 59, sec:scope, the integrated class (record IC, docs/predictions/2026-11-01_integrated_class.md, GO, with
its Caps erratum): regression guard.  Written 2026-11-02 as a separate step after the paper (CLAUDE.md Sec. 9,
guard-writing rule).

Claim guarded: as a function of the bond angle, the integrated (XY|XZ) is not in the triangular class (every
composition factor of order one), so it leaves the two-centre class; the evidence is constant components whose
minimal operator, pulled back to w = |Y - Z|, has no hyperexponential solution and stops at order 6.

Instrument: verbatim copies of the frozen IC instrument in tests/paper59_intclass_support/ (hashes checked
below).  Data: tests/data/paper59_intclass_2026_11_02.json, built by debug/build_intclass_guard_data_2026_11_02.py
(case W1, prime 1073741833, L = 600, the two ln components; nothing is imported from debug/).

  * the copies are the frozen instrument (sha256);
  * component (-231/50, ln 5/2) of W1: minimal x-operator of order 6, both operator sets peel to a clean STUCK
    residual of order 6, no order-one left factor either; every exponent at every finite singular point is an
    integer, and every genuine (non-apparent) singular factor has degree one, so the search is complete there
    (the Caps erratum's condition: a miss needs a factor of degree >= 2 with a non-integer exponent class).  The
    residual also has apparent factors of degree 2 and 12 (exponents 0-4, 6); integer exponents there are
    absorbed by the polynomial part, so they cannot hide a solution;
  * controls run here, at the controls' geometry (9/5, 13/5): rho^L/(L+1) (a logarithm, triangular) peels
    to TRIANGULAR; the central-binomial series P1 (an elliptic-type period) is STUCK;
  * the erratum's planted control 4(c - w^2) y'' - 4 w y' + y (solutions (w -+ s)^(1/2), triangular) peels
    at c = 2, a square mod p, and is returned STUCK at c = 5, a non-square: the gap is pinned, so the scope of
    the claim stays honest (it rests on degree-one factors).
Scope: one component, one case, one prime, L = 600 (the record used L = 2000, two cases, two primes, six
components each); Qbar-independence of the 31 constants and minimality of the operators are assumed, as in the
record.  About 3 s.

Fire-tested 2026-11-02 with debug/qa/fire_test.py; plants 1 and 3-5 FIRED:
  1. the component replaced by the triangular control rho^L/(L+1) (planted triangular operator);
  2. (stayed green, correctly) the pullback taken with the sign of a2 flipped: that is w -> i w, a linear change
     of variable over Qbar, which cannot change whether every composition factor has order one;
  3. the conjugate-root control's expected verdicts swapped (c = 5 expected to peel);
  4. the support copy of intclass_lib.py edited by one character;
  5. the genuine-factor check widened to every factor (it then fails on the apparent degree-2 factor, which is
     how the first draft of this guard read the erratum's looser wording; see the record's addendum).
"""
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import flint
import pytest

HERE = Path(__file__).parent
SUP = HERE / 'paper59_intclass_support'
# sha256 after CRLF -> LF (git may rewrite line endings on checkout).  intclass_lib.py was frozen with CRLF
# endings, raw sha256 b1655d64...67a6 as in commit 40093e7c; the other three were frozen with LF.
FROZEN = {'intclass_lib.py': '2a142e5b5667b432c7de4be885f5b6767daf3d209c86261076dc6ae5e4f641d4',
          'intclass_score.py': 'aaf700e482f53c0a204939c08e6e12c9786d544e48c77e176503f21e5add91ca',
          'na_modguess.py': '84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d',
          'na_modanalyze.py': 'de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005'}

D = json.loads((HERE / 'data' / 'paper59_intclass_2026_11_02.json').read_text())
P = D['prime']
COMP = "(Fraction(-231, 50), (('ln', Fraction(5, 2)),))"


def _load():
    sys.path.insert(0, str(SUP))
    try:
        spec = importlib.util.spec_from_file_location('p59_intclass_score', SUP / 'intclass_score.py')
        sc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sc)
    finally:
        sys.path.remove(str(SUP))
    return sc, sc.I


SC, I = _load()


def test_support_copies_are_the_frozen_instrument():
    """Wrong answer rejected: an instrument edited after the record was scored."""
    for name, h in FROZEN.items():
        assert hashlib.sha256((SUP / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == h, name


def pulled_back_peel(c_x, D1, D2):
    a0 = I.to_fp((D1 ** 2 + D2 ** 2) / (2 * D1 * D2), P)
    a2 = I.to_fp(Fraction(-1) / (2 * D1 * D2), P)
    cw = I.pullback_quadratic(c_x, a0, a2, P)
    return I.peel(cw, P, verbose=False)


def test_ln_component_is_stuck_at_order_six():
    """Wrong answers rejected: a component whose operator peels to order 0 (triangular); a residual of a
    different order; a residual with an order-one left factor; a stuck verdict that rests on a singular factor
    of degree two or more that is genuine, or on non-integer exponents (where the search can miss without a
    flag)."""
    D1, D2 = Fraction(D['case']['D1']), Fraction(D['case']['D2'])
    seq = D['components'][COMP]
    ops = SC.operator_sets(seq, P, 12)
    assert ops is not None and ops['R0'] == 6, ops and ops['R0']
    for key in ('A', 'B'):
        r = pulled_back_peel(ops[key], D1, D2)
        assert r['verdict'] == 'STUCK', (key, r['verdict'])
        assert r['residual_order'] == 6 and r['adjoint_has_order1'] is False, key
        cl = I.classify_residual(r, P)
        assert cl['clean'] and cl['order'] == 6, (key, cl['flags'])
        pts = I.finite_points(r['residual'], P)
        for pt in pts:
            assert all(k == 'int' for (k, v), m in pt['croots']), (key, pt['croots'])
            if not pt['apparent_type']:
                assert pt['fdeg'] == 1, (key, pt['fdeg'])
        assert sum(1 for pt in pts if not pt['apparent_type']) == 5, key


CTRL_D1, CTRL_D2 = Fraction(9, 5), Fraction(13, 5)
RHO = CTRL_D1 / CTRL_D2


def seq_log(n):
    return [I.to_fp(RHO ** L / (L + 1), P) for L in range(n)]


def seq_p1(n):
    out, cb = [], Fraction(1)
    for L in range(n):
        out.append(I.to_fp(RHO ** L * cb, P))
        cb = cb * Fraction((2 * L + 1) * (2 * L + 2), (L + 1) ** 2 * 4)
    return out


def test_controls_triangular_and_stuck():
    """Wrong answers rejected: an instrument that calls everything STUCK (the logarithm control would fail) or
    everything TRIANGULAR (the elliptic-type control would fail)."""
    for seq, expect in ((seq_log(300), 'TRIANGULAR'), (seq_p1(300), 'STUCK')):
        m = I.minimal_operator(seq, P, 16, 30)
        assert m is not None
        assert pulled_back_peel(m['c'], CTRL_D1, CTRL_D2)['verdict'] == expect, expect


CONJ_EXPECT = {2: 'TRIANGULAR', 5: 'STUCK'}


def test_conjugate_root_gap_is_pinned():
    """4(c - w^2) y'' - 4 w y' + y: triangular; peels when c is a square mod p, missed (STUCK) when it is not.
    Wrong answer rejected: a scope statement that the search is complete for factors of any degree."""
    assert pow(2, (P - 1) // 2, P) == 1 and pow(5, (P - 1) // 2, P) == P - 1
    for c, expect in CONJ_EXPECT.items():
        op = [flint.nmod_poly([1], P), flint.nmod_poly([0, -4 % P], P), flint.nmod_poly([4 * c % P, 0, -4 % P], P)]
        assert I.peel(op, P, verbose=False)['verdict'] == expect, (c, expect)
