"""Record AE2 (docs/predictions/2026-11-15_aa_bc_equal_arms_2.md): equal arms at E3 (zB D = 9) and E4.
Registers E3, E4 in THIS process only; reuses aaeq_run's AE0 pieces.

    python debug/aaeq2_run.py ap0                 : AE2-0 (ii) quadrature, (iii) continuity at 400 digits
    python debug/aaeq2_run.py p1 <case> <name> [Lmax dps] : P1, P3, P4 for one integral at one case
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aaeq_run as R  # noqa: E402  (registers E1, E2 too; harmless)

X, P, E = R.X, R.P, R.E
DATA = HERE / 'data'
X.CASES['E3'] = dict(D1='3', D2='3', zB='3', zC='1', zA='9/10')
X.CASES['E4'] = dict(D1='5/2', D2='5/2', zB='7/10', zC='13/10', zA='11/10')
EQ = ('E3', 'E4')
INTS = R.INTS
PTS = {'E3': ('1/2', '0'), 'E4': ('-1/3', '4/5')}
XKEY = {'1/2': 'x=0.5', '0': 'x=0', '-1/3': 'x=-0.333333', '4/5': 'x=0.8'}
TOL = {'E3': 1e-9, 'E4': 2e-10}
BAND = {('E3', 'S'): (-8.25, -7.75), ('E4', 'S'): (-8.25, -7.75), ('E3', 'PD'): (-10.25, -9.75),
        ('E4', 'PD'): (-10.25, -9.75), ('E4', 'PC'): (-10.25, -9.75), ('E3', 'PC'): (-12.3, -11.7)}


def ap0():
    ok = True
    K = E.Field(None)
    mp.mp.dps = 80
    for case in EQ:
        for name, integral in INTS.items():
            f = P.coefficients(case, 8, K, integral)
            ex = R.numeric(f, case, integral)
            for L in (0, 3, 8):
                q = R.quad_FL(case, name, L)
                with mp.workdps(40):
                    rel = abs(q / ex[L] - 1)
                good = rel < mp.mpf('1e-18')
                ok &= bool(good)
                print(case, name, 'L', L, '(ii) rel', mp.nstr(rel, 3), good, flush=True)
    mp.mp.dps = 400
    eps = Fraction(1, 10 ** 40)
    for case in EQ:
        D1, D2, zB, zC, zA = P.params(case)
        for name, integral in INTS.items():
            f0 = P.coefficients(case, 20, K, integral)
            f1 = P.coefficients(case, 20, K, integral, params_override=(D1 - eps, D2, zB, zC, zA))
            a0, a1 = R.numeric(f0, case, integral), R.numeric(f1, case, integral)
            worst = max(abs(a1[L] / a0[L] - 1) for L in range(21))
            good = worst < mp.mpf('1e-35')
            ok &= bool(good)
            print(case, name, '(iii) continuity worst', mp.nstr(worst, 3), good, flush=True)
    print('AE2-0 self-checks pass:', ok)


def p1(case, name, Lmax=600, dps=2000):
    ctrl = json.loads((DATA / 'aaeq2_controls.json').read_text())
    integral = INTS[name]
    K = E.Field(None)
    t0 = time.time()
    f = P.coefficients(case, Lmax, K, integral)
    sec = round(time.time() - t0)
    mp.mp.dps = dps
    F = R.numeric(f, case, integral)
    vals, diffs = {}, {}
    for xs in PTS[case]:
        x = R.mfq(xs)
        val = sum(F[L] * R.legendre(L, x) for L in range(Lmax + 1))
        ref = ctrl['referee'][f'{case}_{XKEY[xs]}']
        vals[xs] = mp.nstr(val, 30)
        diffs[xs] = dict(R2=float(val - mp.mpf(repr(ref[name]['34']))), R1=float(val - mp.mpf(repr(ref['R1'][name]))))
    tol = TOL[case]
    P1 = all(abs(d['R2']) <= tol for d in diffs.values())
    beta = {L: float(mp.log(abs(F[L + 1] / F[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1))) for L in (220, 400, 580)}
    b = BAND[(case, name)]
    P3 = b[0] <= beta[580] <= b[1]
    if case == 'E3' and name == 'PC':
        P3 = all(b[0] <= beta[L] <= b[1] for L in (220, 400, 580)) and abs(beta[580] - beta[220]) < 0.05
    rec = dict(values=vals, diffs=diffs, P1=P1, beta=beta, P3=P3,
               signs_L3_12=[int(mp.sign(v)) for v in F[3:13]], signs_L575_585=[int(mp.sign(v)) for v in F[575:586]],
               F_580=mp.nstr(F[580], 12), engine_s=sec, dps=dps, Lmax=Lmax)
    if case == 'E3' and name == 'PC':
        pred, lo, hi = ctrl['local_expansion']['E3_PC_pred_F'], 0.9, 1.1
    else:
        pred, lo, hi = ctrl['amplitude'][case][name]['pred_F'], 0.98, 1.02
    rec['P4_ratio'] = {L: float(F[int(L)] / mp.mpf(v)) for L, v in pred.items()}
    rec['P4'] = lo <= rec['P4_ratio']['580'] <= hi
    rec['F_tail'] = {L: mp.nstr(F[L], 25) for L in (220, 400, 579, 580, 581)}   # for PC/PD ratios at E3
    (DATA / f'aaeq2_p1_{dps}_{case}_{name}.json').write_text(json.dumps(rec, indent=1))   # one file per run
    print(case, name, json.dumps(rec, indent=1), flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'ap0':
        ap0()
    elif cmd == 'p1':
        p1(sys.argv[2], sys.argv[3], *(int(v) for v in sys.argv[4:6]))
