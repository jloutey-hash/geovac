"""Scoring driver for PY, PXp, DA (engine debug/aaa_exact.py).  Written before the record; NOT run on the new
integrals by the engine author (pre-registration).  'PA' is accepted as a positive control (record AP, scored).

    python debug/aaa_run.py p2 <prime> [Lmax=1000] [integrals=PY,PXp,DA]
        constant sets mod <prime>, L = 0..Lmax, per integral and per channel (PY/PXp: a = P_L part, b = P_L^2 part,
        and their union), L-independence from L = 2 (b: from L = 3, since b_2 is the first nonzero), inside AB's
        candidate set C, gamma, equality with AB's 21 -> data/aaa_p2_<prime>.json
    python debug/aaa_run.py p1p3 [Lmax=240] [dps=2000] [integrals=PY,PXp,DA]
        exact coefficients over Q to Lmax; summed values at the referee points (R2 n = 34 from
        data/aaa_controls.json) with tail terms; beta_eff(L) = log|c_{L+1}/c_{L-1}| / log((L+1)/(L-1)) on
        c_L / rho^L for every channel series at L = 20, 40, ..., Lmax - 20; signs -> data/aaa_p1p3_<dps>.json
No pass/fail bands are coded here: the record sets them.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aaa_exact as A  # noqa: E402
import aabc_exact as X  # noqa: E402
import aabc_run as ABR  # noqa: E402
import na_exact as E  # noqa: E402

DATA = HERE / 'data'
POINTS = {'A1': ('1/2', '-1/3'), 'A2': ('0', '4/5')}


def channels(case: str, Lmax: int, K, integral: str, eng=None) -> dict:
    if integral in ('PY', 'PXp'):
        a, b = A.coefficients(case, Lmax, K, integral, eng=eng)
        return {'a': a, 'b': b}
    return {'F': A.coefficients(case, Lmax, K, integral, eng=eng)}


def p2(p: int, Lmax: int = 1000, integrals=('PY', 'PXp', 'DA')) -> None:
    out = {}
    K = E.Field(p)
    for case in A.CASES:
        na_keys = set(k for q in X.pieces(case, 40, K, which=('NA',))['NA'] for k in q)
        C = ABR.candidate_set(case, na_keys)
        abf = DATA / f'aabc_p2_{p}.json'
        ab21 = set(json.loads(abf.read_text())[case]['constants']) if abf.exists() else None
        eng = A.Engine(case, K, Lmax + 2)
        for integral in integrals:
            t0 = time.time()
            ch = channels(case, Lmax, K, integral, eng)
            sec = round(time.time() - t0)
            recs = {}
            series = dict(ch)
            if len(ch) == 2:
                series['a+b'] = [set(x) | set(y) for x, y in zip(ch['a'], ch['b'])]
            for name, f in series.items():
                sets = [frozenset(q) for q in f]
                L0 = 3 if name == 'b' else 2
                S0 = sets[L0]
                changes = [L for L in range(L0 + 1, len(sets)) if sets[L] != S0]
                union = set().union(*sets[L0:])
                first_fixed = next((L for L in range(len(sets)) if all(sets[M] == S0 for M in range(L, len(sets)))),
                                   None)
                rec = dict(L0=L0, n_constants_L0=len(S0), L_where_set_changes=changes[:20], fixed=not changes,
                           first_L_fixed=first_fixed, sizes_L0_to_3=[len(sets[L]) for L in range(4)],
                           n_union=len(union), inside_C=not (union - C), outside_C=sorted(map(repr, union - C)),
                           gamma_present=[repr(k) for k in union if k[1] == 'g'], constants=sorted(map(repr, union)))
                if ab21 is not None:
                    u = set(map(repr, union))
                    rec['equals_AB21'] = u == ab21
                    rec['vs_AB21'] = dict(missing=sorted(ab21 - u), extra=sorted(u - ab21))
                recs[name] = rec
            out[f'{case}_{integral}'] = dict(seconds=sec, Lmax=Lmax, series=recs)
            print(case, integral, p, sec, 's',
                  json.dumps({n: {k: v for k, v in r.items() if k != 'constants'} for n, r in recs.items()}, indent=1),
                  flush=True)
            (DATA / f'aaa_p2_{p}.json').write_text(json.dumps(out, indent=1))


def legendre_pair(Lmax: int, x):
    """P_L(x) and P_L^2(x) = (1 - x^2) P_L''(x) (no Condon-Shortley phase), L = 0..Lmax, |x| < 1."""
    P = [mp.mpf(1), x]
    for L in range(1, Lmax):
        P.append(((2 * L + 1) * x * P[L] - L * P[L - 1]) / (L + 1))
    om = 1 - x * x
    P2 = [mp.mpf(0), mp.mpf(0)] + [2 * x * L * (P[L - 1] - x * P[L]) / om - L * (L + 1) * P[L]
                                   for L in range(2, Lmax + 1)]
    return P[:Lmax + 1], P2


def referee_values() -> dict:
    d = json.loads((DATA / 'aaa_controls.json').read_text())['referee']
    out = {}
    for k, v in d.items():
        if k == 'floors':
            continue
        case, integral, xs = k.split('_', 2)
        out[(case, integral, v['x'])] = dict(R2_34=v['R2']['34'], R2_30=v['R2']['30'], R1=v['R1']['120,120,144'])
    return out


def p1p3(Lmax: int = 240, dps: int = 2000, integrals=('PY', 'PXp', 'DA')) -> None:
    out = {}
    K = E.Field(None)
    ref = referee_values() if (DATA / 'aaa_controls.json').exists() else {}
    for case in A.CASES:
        D1, D2, zB, zC, zA = A.params(case)
        eng = A.Engine(case, K, Lmax + 2)
        xsq = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
        for integral in integrals:
            t0 = time.time()
            ch = channels(case, Lmax, K, integral, eng)
            sec = round(time.time() - t0)
            mp.mp.dps = dps
            keys = set(k for f in ch.values() for q in f for k in q)
            consts = {k: E.value(k, mp) for k in keys}
            tw = A.twopiN(case)
            V = {n: [tw * sum(mp.mpf(int(v.p)) / int(v.q) * consts[k] for k, v in q.items()) for q in f]
                 for n, f in ch.items()}
            if 'b' in V:   # record AS (added by the PM before registration): the total as a P_L series,
                # P_L^2 = (1 - x^2) P_L'' = -L(L-1) P_L + 2 sum_{k < L, k = L mod 2} (2k+1) P_k  (no CS phase)
                b_ = V['b']
                tail = [mp.mpf(0), mp.mpf(0)]          # running sums of b over L > j, by parity
                cc = [mp.mpf(0)] * (Lmax + 1)
                for j in range(Lmax, -1, -1):
                    cc[j] = V['a'][j] - j * (j - 1) * b_[j] + 2 * (2 * j + 1) * tail[j % 2]
                    tail[j % 2] += b_[j]
                V_total = cc
            else:
                V_total = None
            xsm = mp.mpf(xsq.numerator) / xsq.denominator
            rho = 1 / (xsm + mp.sqrt(xsm * xsm - 1))
            rec = dict(engine_s=sec, dps=dps, Lmax=Lmax, values={}, diffs={}, beta={}, signs={})
            for xs_ in POINTS[case]:
                xr = mp.mpf(Fraction(xs_).numerator) / Fraction(xs_).denominator
                Pl, P2 = legendre_pair(Lmax, xr)
                basis = {'a': Pl, 'b': P2, 'F': Pl}
                val = sum(V[n][L] * basis[n][L] for n in V for L in range(Lmax + 1))
                tail = max(abs(V[n][Lmax] * basis[n][Lmax]) for n in V)
                rec['values'][xs_] = mp.nstr(val, 30)
                if V_total is not None:     # exact identity at truncation: checks the re-expansion formula
                    rec.setdefault('total_reexpansion_check', {})[xs_] = mp.nstr(
                        abs(sum(V_total[L] * Pl[L] for L in range(Lmax + 1)) - val), 3)
                r = ref.get((case, integral, xs_))
                rec['diffs'][xs_] = dict(tail=mp.nstr(tail, 5))
                if r:
                    rec['diffs'][xs_].update(R2_34=float(val - mp.mpf(r['R2_34'])),
                                             R2_30=float(val - mp.mpf(r['R2_30'])), R1=float(val - mp.mpf(r['R1'])))
            series = dict(V)
            if V_total is not None:
                series['c'] = V_total
                amp = json.loads((DATA / 'aaa_controls.json').read_text()).get('amp_toys', {}).get(case)
                if amp:
                    R = {}
                    for L in (219, 220, 221):
                        if L <= Lmax:
                            ta, tb = mp.mpf(amp['ta'][str(L)]), mp.mpf(amp['tb'][str(L)])
                            R[L] = float((V['a'][L] / ta) / (V['b'][L] / tb))
                    rec['P4_ratio_a_over_b'] = R
                    rec['P4_predicted_x_star_sq_minus_1'] = amp['x_star_sq_minus_1']
                    if Lmax >= 221:
                        rec['c_over_tc'] = {L: float(V_total[L] / mp.mpf(amp['tc'][str(L)])) for L in (219, 220, 221)}
            for n, F in series.items():
                G = [F[L] / rho ** L for L in range(Lmax + 1)]
                lo = 3 if n == 'b' else 1
                top = Lmax if n != 'c' else Lmax - 160     # c_j needs the b tail beyond j (geometric, rho^160)
                rec['beta'][n] = {L: float(mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1)))
                                  for L in sorted(set(list(range(20, top, 20)) + [220] + [Lmax - 20]))
                                  if L - 1 >= lo and L + 1 <= (Lmax if n != 'c' else top)}
                rec['signs'][n] = dict(L_ge_20=sorted(set(int(mp.sign(v)) for v in F[20:])),
                                       L_ge_5=sorted(set(int(mp.sign(v)) for v in F[5:])))
            rec['F_tail'] = {n: {L: mp.nstr(F[L], 25) for L in (219, 220, 221, Lmax - 2, Lmax - 1, Lmax) if 0 <= L <= Lmax}
                             for n, F in series.items()}   # record AS, P0 precision gate
            out[f'{case}_{integral}'] = rec
            print(case, integral, json.dumps(rec, indent=1), flush=True)
            (DATA / f'aaa_p1p3_{dps}.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    cmd = sys.argv[1]
    args = sys.argv[2:]
    ints = tuple(args[-1].split(',')) if args and not args[-1].isdigit() else ('PY', 'PXp', 'DA')
    nums = [int(a) for a in args if a.isdigit()]
    if cmd == 'p2':
        p2(nums[0], *(nums[1:2] or [1000]), integrals=ints)
    elif cmd == 'p1p3':
        p1p3(*(nums[:2] or [240, 2000]), integrals=ints)
