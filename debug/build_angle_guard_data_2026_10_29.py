"""Builds tests/data/paper59_angle_2026_10_29.json for tests/test_paper59_bond_angle_operator.py (a copy of stored
results, prime 2^61 - 1 only).  Contents:
  na: N1 generic sequence L <= 400 (random.Random(2026) weights, as debug/na_points.py) and NA's (8, 31) null vector;
  b2: W1 generic sequence L <= 400 (weights as debug/ao3_run.py) and B2's derived recurrence (debug/data);
  w1: exact rational coefficients of W1 at L in LS with their constant keys, and the K-integral constants at 700 digits.
    python debug/build_angle_guard_data_2026_10_29.py
"""
import json
import random
from fractions import Fraction

import flint
import mpmath as mp

import ao_exact as X
import na_exact as E
import na_modguess as G

p = (1 << 61) - 1
LS = [1, 2, 3, 4, 5, 49, 50, 51, 99, 100, 101]


def gen_of(Qs, keyf):
    keys = sorted(set(k for q in Qs for k in q), key=keyf)
    rnd = random.Random(2026)
    w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
    return [int(sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p))) for q in Qs], keys, w


def fr(x):
    return None if x is None else f'{Fraction(x).numerator}/{Fraction(x).denominator}'


def ser_key(k):
    e0, facs = k
    return [fr(e0), [[f[0]] + [fr(a) if not isinstance(a, str) else a for a in f[1:]] for f in facs]]


# NA, case N1
Qn = E.coefficients('1.8', '2.6', '1.3', '1.0', 400, E.Field(p))
gen_n, _, _ = gen_of(Qn, str)
B = G.Basis([flint.nmod(v, p) for v in gen_n], p)
nul, ks, (N, nrow) = G.nullity(B, 8, 31, 30)
assert nul == 1
na = dict(case='N1', D=['9/5', '13/5'], z=['13/10', '1'], gen=gen_n, keys=[list(k) for k in ks],
          a=[int(N[i, 0]) for i in range(N.nrows())])

# B2, case W1
Qw = X.coefficients('1.7', '2.5', '1.4', '1.1', 400, E.Field(p))
gen_w, keys_w, w = gen_of(Qw, repr)
rec = json.load(open(f'data/b2_recurrence_{p}_4200.json'))
b2 = dict(case='W1', gen=gen_w, polys=rec['polys'], order=rec['order'])

# exact W1 rationals at LS, and the weights of their keys (so the test ties them to the sequence)
Qx = X.coefficients('1.7', '2.5', '1.4', '1.1', max(LS), E.Field(None))
used = sorted(set(k for L in LS for k in Qx[L]), key=repr)
mp.mp.dps = 700
kvals = {}
for k in used:
    for f in k[1]:
        if f[0] not in ('g', 'ln', 'E1') and repr(f) not in kvals:
            kvals[repr(f)] = dict(factor=ser_key((Fraction(0), (f,)))[1][0], value=mp.nstr(X.factor_value(f, mp), 690))
w1 = dict(keys=[ser_key(k) for k in used], weights=[int(w[k]) for k in used],
          coeffs={str(L): {str(used.index(k)): [str(int(v.p)), str(int(v.q))] for k, v in Qx[L].items()} for L in LS},
          K=list(kvals.values()))
json.dump(dict(note='built by debug/build_angle_guard_data_2026_10_29.py; prime 2^61-1', p=p, na=na, b2=b2, w1=w1),
          open('../tests/data/paper59_angle_2026_10_29.json', 'w'))
print('done; K constants', len(kvals), 'keys', len(used))
