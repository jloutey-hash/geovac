"""Record AP scoring driver (docs/predictions/2026-11-14_aa_bc_p_orbitals.md), Phase 2.
    python debug/aapc_run.py p1p3 [Lmax] [dps]   : P1 (value at the referee points) + P3 (decay), PC and PA, both
                                                   cases, exact split at dps digits -> data/aapc_p1p3_<dps>.json
    python debug/aapc_run.py p2 <prime>          : constant sets L = 0..1000 mod prime -> data/aapc_p2_<p>.json
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aapc_exact as P  # noqa: E402
import aabc_exact as X  # noqa: E402
import aabc_run as ABR  # noqa: E402
import na_exact as E  # noqa: E402

DATA = HERE / 'data'
REF = {('A1', 'PC'): dict(R2={'22': '-0.2089629965563952', '26': '-0.20896299648621663'},
                          R2_stated='-0.20896299649', R1_finest='-0.20896299577928962', tolR2=3e-10, tolR1=5e-9),
       ('A2', 'PC'): dict(R2={'22': '-0.21839796649545454', '26': '-0.2183979664314794'},
                          R2_stated='-0.21839796643', R1_finest='-0.21839796926538216', tolR2=3e-10, tolR1=5e-9),
       ('A1', 'PA'): dict(R2={'22': '0.203231256920345', '26': '0.20323125696311292'},
                          R2_stated='0.20323125696', R1_finest='0.20323125699555303', tolR2=2e-10, tolR1=3e-10),
       ('A2', 'PA'): dict(R2={'22': '0.21952561492962194', '26': '0.2195256150153055'},
                          R2_stated='0.21952561502', R1_finest='0.21952561498101267', tolR2=2e-10, tolR1=3e-10)}


def mpq(v):
    return mp.mpf(int(v.p)) / int(v.q)


def p1p3(Lmax=240, dps=2000):
    out = {}
    K = E.Field(None)
    for case in P.CASES:
        D1, D2, zB, zC, zA = P.params(case)
        xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
        for integral in ('PC', 'PA'):
            t0 = time.time()
            if integral == 'PA':
                f, pp = P.coefficients(case, Lmax, K, 'PA', pieces_out=True)
            else:
                f, pp = P.coefficients(case, Lmax, K, 'PC'), None
            sec = round(time.time() - t0)
            mp.mp.dps = dps
            keys = set(k for q in f for k in q)
            if pp:
                keys |= set(k for d in pp for part in d.values() for k in part)
            consts = {k: E.value(k, mp) for k in keys}
            tw = P.twopiN(case, integral)
            F = [tw * sum(mpq(v) * consts[k] for k, v in q.items()) for q in f]
            xr = mp.mpf(Fraction(P.REF_X[case]).numerator) / Fraction(P.REF_X[case]).denominator
            Pl = [mp.mpf(1), xr]
            for L in range(1, Lmax):
                Pl.append(((2 * L + 1) * xr * Pl[L] - L * Pl[L - 1]) / (L + 1))
            val = sum(F[L] * Pl[L] for L in range(Lmax + 1))
            ref = REF[(case, integral)]
            diffs = {('R2_n' + k): float(val - mp.mpf(v)) for k, v in ref['R2'].items()}
            diffs['R2_stated'] = float(val - mp.mpf(ref['R2_stated']))
            diffs['R1_finest'] = float(val - mp.mpf(ref['R1_finest']))
            xsm = mp.mpf(xs.numerator) / xs.denominator
            rho = 1 / (xsm + mp.sqrt(xsm * xsm - 1))

            def beta(Fs):
                G = [Fs[L] / rho ** L for L in range(Lmax + 1)]
                return {L: float(mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1)))
                        for L in list(range(20, Lmax, 20)) + [Lmax - 20]}
            rec = dict(value=mp.nstr(val, 30), tail_term=mp.nstr(abs(F[Lmax] * Pl[Lmax]), 5), diffs=diffs,
                       P1=abs(diffs['R2_stated']) <= ref['tolR2'] and abs(diffs['R1_finest']) <= ref['tolR1'],
                       beta_f=beta(F), signs_L_ge_20=sorted(set(int(mp.sign(v)) for v in F[20:])),
                       signs_L_ge_5=sorted(set(int(mp.sign(v)) for v in F[5:])),
                       n_constants=len(set(k for q in f for k in q)), engine_s=sec, dps=dps, Lmax=Lmax,
                       F_first=[mp.nstr(v, 15) for v in F[:6]])
            b220 = rec['beta_f'][Lmax - 20]
            rec['P3'] = (-4.3 <= b220 <= -3.85) and len(rec['signs_L_ge_20']) == 1
            if pp:
                for part_set, name in ((('U2', 'B'), 'dD2sq_U'), (('W0',), 'W0')):
                    Fp = [tw * sum(mpq(v) * consts[k] for part in part_set for k, v in d[part].items()) for d in pp]
                    rec['beta_' + name] = beta(Fp)
            out[f'{case}_{integral}'] = rec
            print(case, integral, json.dumps(rec, indent=1), round(time.time() - t0), 's', flush=True)
            (DATA / f'aapc_p1p3_{dps}.json').write_text(json.dumps(out, indent=1))


def p2(p, Lmax=1000):
    out = {}
    K = E.Field(p)
    for case in P.CASES:
        na_keys = set(k for q in X.pieces(case, 40, K, which=('NA',))['NA'] for k in q)
        C = ABR.candidate_set(case, na_keys)
        ab21 = None
        abf = DATA / f'aabc_p2_{p}.json'
        if abf.exists():
            ab21 = set(json.loads(abf.read_text())[case]['constants'])
        for integral in ('PC', 'PA'):
            t0 = time.time()
            f = P.coefficients(case, Lmax, K, integral)
            sec = round(time.time() - t0)
            sets = [frozenset(q) for q in f]
            L0 = 3 if integral == 'PC' else 2
            S0 = sets[L0]
            changes = [L for L in range(L0 + 1, len(sets)) if sets[L] != S0]
            union = set().union(*sets[L0:])
            rec = dict(seconds=sec, L0=L0, n_constants_L0=len(S0), L_where_set_changes=changes[:20],
                       n_union=len(union), outside_C=sorted(map(repr, union - C)), C_size=len(C),
                       C_unused=sorted(map(repr, C - union)),
                       new_E1_or_log=sorted(repr(k) for k in union - C if k[1] in ('E1', 'ln')),
                       gamma_present=[repr(k) for k in union if k[1] == 'g'],
                       constants=sorted(map(repr, union)),
                       sets_below_L0=[sorted(map(repr, sets[L])) for L in range(L0)])
            if ab21 is not None:
                u = set(map(repr, union))
                rec['equals_AB21'] = u == ab21
                rec['vs_AB21'] = dict(missing=sorted(ab21 - u), extra=sorted(u - ab21))
            rec['P2a'] = not changes
            rec['P2b'] = not (union - C)
            rec['P2c'] = not rec['gamma_present']
            out[f'{case}_{integral}'] = rec
            print(case, integral, p, json.dumps({k: v for k, v in rec.items() if k != 'constants'}, indent=1),
                  flush=True)
            (DATA / f'aapc_p2_{p}.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'p1p3':
        p1p3(*(int(a) for a in sys.argv[2:4]))
    elif cmd == 'p2':
        p2(int(sys.argv[2]))
