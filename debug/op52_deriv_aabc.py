"""OP52 derivation memo: the endpoint-germ predictions on (AA|BC) itself (record MO's f = NA - Y_mu - zA W_mu).

Builds the per-piece minimal operators of f's five endpoint pieces from MO's cached sequences (debug/data/s1_coeffs_*,
s1_controls.piece_operator from the frozen support copy, unchanged) and compares their local branch series at x*
with the pinch-germ endpoint parts (debug/op52_deriv_nextorder.py's germ formula, extended to |a|^9):
  D1:  B_S - B_N = -G^{D1}_{V_A} / G^{D1}_{1/s}.lead   -> [t^{5/2}] = 0 (v_1(V_A) = 0), [t^{7/2}], [t^{9/2}] predicted;
  D2:  normalized V_A germ, rho_D2 = [t^{9/2}]/[t^{7/2}];  collision residual rho_D2 - rho_D1 (9/2 iff nonzero).
  P2:  the two Y/dY halves carry the same normalized series; so do the two D2 halves.
Usage: python debug/op52_deriv_aabc.py <case> <p>    -> debug/data/op52_deriv_aabc_<case>_<p>.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
SUP = ROOT / 'tests' / 'paper59_intclass_support'
sys.path.insert(0, str(SUP))
import s1_controls as S1          # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'debug' / f'{name}.py')
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


loc = _load('op52_deriv_local')
nxt = _load('op52_deriv_nextorder')
CASES = {'A1': ('9/5', '13/5', '13/10', '1', '9/10'), 'A2': ('3/2', '12/5', '1', '6/5', '4/5')}
XSTAR = {'A1': Fraction(125, 117), 'A2': Fraction(89, 80)}


def germ(V, D1, D2, zB, zC, nmax=9):
    """coefficients [t^{n/2}] (n = 3..nmax odd) at both endpoints, through relative |a|^2 (N dropped)."""
    K = 2 * D1 * D2
    out = {E: {n: sp.Integer(0) for n in range(3, nmax + 1, 2)} for E in ('D1', 'D2')}
    for p, vp in V.items():
        for q in (1, 3, 5):
            n = p + q + 3
            if n > nmax:
                continue
            gq = (-zC) ** q / sp.factorial(q)
            a, b = sp.Rational(p + 1, 2), sp.Rational(q + 1, 2)
            lead = nxt.J(lambda u: u ** a * (1 - u) ** b, D1, D2, zB)
            for E in out:
                out[E][n] += vp * nxt.Cpq(p, q) * gq * K ** sp.Rational(n, 2) * lead[E]
            if n + 2 <= nmax:
                corr = nxt.J(lambda u: u ** a * (1 - u) ** b * u * (1 - u) *
                             (zB / (2 * nxt.r) - (zB ** 2 - 2 * zB / nxt.r) / (2 * (p + q + 5))), D1, D2, zB)
                for E in out:
                    out[E][n + 2] += vp * nxt.Cpq(p, q) * gq * K ** sp.Rational(n + 2, 2) * corr[E]
    return {E: {k: sp.simplify(v) for k, v in d.items()} for E, d in out.items()}


def main(case, p):
    t0 = time.time()
    D1, D2, zB, zC, zA = (sp.Rational(v) for v in CASES[case])
    mu = 2 * zA
    one = lambda e: repr((Fraction(e), '1', None))                                         # noqa: E731
    keys = {'NA': one(-zB * D1), 'S-': one(-zB * D1 - mu * (D2 - D1)), 'S+': one(-zB * D1 - mu * (D2 + D1)),
            'D2+': one(zC * D1 - (zB + zC) * D2), 'D2-': one(-zC * D1 - (zB + zC) * D2)}
    seqs = json.loads((ROOT / 'debug' / 'data' / f's1_coeffs_{case}_{p}.json').read_text())['seqs']
    ser, orders = {}, {}
    for nm, k in keys.items():
        op = S1.piece_operator(seqs[k], p)
        orders[nm] = op['order']
        E, logs = loc.local_halfint_space(op['c'], XSTAR[case], p, N=7)
        assert len(E) == 1 and not logs, (nm, len(E), logs)
        ser[nm] = E[0]
        print(nm, op['order'], [str(x) for x in loc.vals(E)], round(time.time() - t0), flush=True)
    VA = {3: 2 * zA ** 4 / 3, 5: 8 * zA ** 6 / 45}                     # 1/s - S, odd terms (no s^-1, s^1)
    gV = germ(VA, D1, D2, zB, zC)
    gN = germ({-1: sp.Integer(1)}, D1, D2, zB, zC, nmax=3)
    lead = gN['D1'][3]
    pred = {'d5': sp.Integer(0), 'd7': -gV['D1'][7] / lead, 'd9': -gV['D1'][9] / lead,
            'rho_D2': gV['D2'][9] / gV['D2'][7], 'rho_D1': gV['D1'][9] / gV['D1'][7]}
    pred['resid'] = pred['rho_D2'] - pred['rho_D1']
    pred = {k: nxt.frac(v) for k, v in pred.items()}
    d = [(a - b) % p for a, b in zip(ser['S-'], ser['NA'])]
    inv7 = pow(d[3], p - 2, p)
    meas = {'d5': d[2], 'd7': d[3], 'd9': d[4], 'rho_D2': ser['D2+'][4] * pow(ser['D2+'][3], p - 2, p) % p,
            'rho_D1': d[4] * inv7 % p}
    meas['resid'] = (meas['rho_D2'] - meas['rho_D1']) % p
    row = {'orders': orders, 'val_D2': str(loc.vals([ser['D2+']])[0]),
           'P2_S_halves_equal': ser['S-'] == ser['S+'], 'P2_D2_halves_equal': ser['D2+'] == ser['D2-'],
           'pred': {k: str(v) for k, v in pred.items()},
           'match': {k: meas[k] == loc.fp(pred[k], p) for k in pred}, 'time': round(time.time() - t0)}
    print(json.dumps(row, indent=1))
    (ROOT / 'debug' / 'data' / f'op52_deriv_aabc_{case}_{p}.json').write_text(json.dumps(row, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
