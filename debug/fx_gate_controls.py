"""Gate FX0 (debug/fx_design_memo.md section 5) for the depth-2 reducer (debug/fx_reducer.py, debug/fx_exact.py).
Records booleans, agreement magnitudes and timings only; no constant-set composition, no recurrence, no operator.

  V1  exact G_L = F_L/pi^2 (fmpq, L <= 30) evaluated at high precision vs numeric fc_engine F_L (levels A and B),
      cases M1 and M2r; gate <= 1e-25 relative.  Evaluation stability: the same sum at two working precisions.
  V2  a = b, zZ != zW (K4A1 = record AB case A1): G_L zA^3 / 2 == aabc_exact A1 coefficients exactly, L <= 40
      (both sides with ln of rationals expanded over primes).
  V3  mod-p splits on two primes == the fmpq split reduced mod p, L <= 40, cases M1, M2r, K4A1.
  V4  record NA's (R, Q0) = (8, 31) on both primes (na_exact + na_modguess, as aabc_na_posctrl).
  T   per-L cost mod p (case M1, both primes) at L = 100, 300; extrapolated cost to L = 1500.
-> debug/data/fx_gate_controls.json   (written after every stage)
    python debug/fx_gate_controls.py [stage ...]     stages: V1 V2 V3 V4 T (default all)
"""
import json
import math
import random
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import flint
import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import na_exact as NE        # noqa: E402
import na_modguess as G      # noqa: E402
import aabc_exact as AB      # noqa: E402
import fc_engine as FC       # noqa: E402
import fx_reducer as X       # noqa: E402
import fx_exact as FX        # noqa: E402

OUTF = HERE / 'data' / 'fx_gate_controls.json'
PRIMES = [(1 << 61) - 1, 4611686018427387847]
OUT = json.loads(OUTF.read_text()) if OUTF.exists() else {}


def save(key, val):
    OUT[key] = val
    OUTF.write_text(json.dumps(OUT, indent=1))
    print(key, json.dumps(val)[:600], flush=True)


_EXACT = {}


def exact(case, Lmax):
    if case not in _EXACT or len(_EXACT[case][0]) <= Lmax:
        t = time.time()
        _EXACT[case] = (FX.coefficients(case, Lmax, NE.Field(None)), time.time() - t)
    return _EXACT[case]


def v1():
    res = {}
    for case in ('M1', 'M2r'):
        Qs, t_ex = exact(case, 40)
        vals = {}
        for dps in (320, 400):
            mp.mp.dps = dps
            t = time.time()
            cache = {}
            vals[dps] = [X.value(Qs[L], mp, cache) for L in range(31)]
            t_eval = time.time() - t
        mp.mp.dps = 400
        stab = max(float(abs(vals[320][L] / vals[400][L] - 1)) for L in range(31))
        row = dict(t_exact_fmpq_L40=round(t_ex, 1), t_eval_dps400=round(t_eval, 1),
                   eval_stability_320_vs_400=stab)
        for lev in ('A', 'B'):
            t = time.time()
            F, info = FC.coefficients(FX.CASES[case], 30, lev)
            rel = []
            for L in range(31):
                nv = mp.mpf(F[L].mid().str(130, radius=False))
                rel.append(float(abs(vals[400][L] * mp.pi ** 2 / nv - 1)))
            row[lev] = dict(max_rel=max(rel), rel_by_L={L: rel[L] for L in (0, 1, 2, 5, 10, 20, 30)},
                            pass_1e_25=max(rel) <= 1e-25, t_fc=round(time.time() - t, 1))
        res[case] = row
    save('V1_exact_vs_fc_numeric', res)


def _conv_ab(q, K):
    out = {}
    for (e0, kind, arg), v in q.items():
        atoms = () if kind == '1' else ((('g',),) if kind == 'g' else ((kind, Fr(arg)),))
        X.add(out, (Fr(e0), atoms), v)
    return X.canon(out, K)


def v2():
    K = NE.Field(None)
    Qs, t_ex = exact('K4A1', 40)
    t = time.time()
    f = AB.combine(AB.pieces('A1', 40, K), Fr(9, 10), K)
    t_ab = time.time() - t
    s = K(Fr(9, 10) ** 3 / 2)
    eq = [{k: v * s for k, v in Qs[L].items()} == _conv_ab(f[L], K) for L in range(41)]
    save('V2_a_eq_b_vs_aabc_exact_A1', dict(all_equal_L_le_40=all(eq), n_unequal=eq.count(False),
                                             normalisation='G_L * zA^3 / 2 = q_AB', t_fx=round(t_ex, 1),
                                             t_aabc=round(t_ab, 1)))


def v3():
    res = {}
    for case in ('M1', 'M2r', 'K4A1'):
        Qs, _ = exact(case, 40)
        res[case] = {}
        for p in PRIMES:
            t = time.time()
            Qp = FX.coefficients(case, 40, NE.Field(p))
            ok = [Qp[L] == X.reduce_mod(Qs[L], p) for L in range(41)]
            res[case][str(p)] = dict(all_equal_L_le_40=all(ok), n_unequal=ok.count(False),
                                     t_modp_L40=round(time.time() - t, 1))
    save('V3_modp_vs_fmpq', res)


def v4():
    res = {}
    for p in PRIMES:
        t0 = time.time()
        Qs = NE.coefficients('1.8', '2.6', '1.3', '1.0', 1000, NE.Field(p))
        keys = sorted(set(k for q in Qs for k in q), key=str)
        rnd = random.Random(2026)
        w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
        gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]
        B = G.Basis(gen, p)
        r = {R: list(G.min_degree_bisect(B, R, 30, len(gen))) for R in (7, 8)}
        res[str(p)] = dict(order7=r[7], order8=r[8], reproduced=(r[7][0] is None and r[8][0] == 31),
                           seconds=round(time.time() - t0))
    save('V4_NA_8_31', res)


def tcost():
    res = {}
    for case, Ltop, primes in (('M1', 300, PRIMES), ('M2r', 100, PRIMES[:1])):
        for p in primes:
            dt = []
            t0 = time.time()
            last = t0
            for L, _ in FX.Engine(case, NE.Field(p)).run(Ltop):
                now = time.time()
                dt.append(now - last)
                last = now
            row = dict(total_s=round(last - t0, 1))
            for Lc in (100, 300):
                if Lc <= Ltop:
                    win = sorted(dt[Lc - 4:Lc + 1])
                    row['per_L_s_at_%d' % Lc] = round(win[2], 3)
                    row['cum_s_to_%d' % Lc] = round(sum(dt[:Lc + 1]), 1)
            if Ltop >= 300:
                a = math.log(row['per_L_s_at_300'] / row['per_L_s_at_100']) / math.log(3)
                c = row['per_L_s_at_300'] / 300 ** a
                tail = c * (1500 ** (a + 1) - 300 ** (a + 1)) / (a + 1)
                row['local_exponent_100_300'] = round(a, 2)
                row['extrapolated_cum_s_to_1500'] = round(row['cum_s_to_300'] + tail)
                row['extrapolated_cpu_h_to_1500'] = round((row['cum_s_to_300'] + tail) / 3600, 2)
            res.setdefault(case, {})[str(p)] = row
            save('T_cost', res)


STAGES = dict(V1=v1, V2=v2, V3=v3, V4=v4, T=tcost)

if __name__ == '__main__':
    for s in (sys.argv[1:] or list(STAGES)):
        t = time.time()
        STAGES[s]()
        print('stage', s, 'done', round(time.time() - t, 1), 's', flush=True)
