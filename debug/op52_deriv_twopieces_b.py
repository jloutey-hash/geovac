"""OP52 memo item 3 (b): is b_L = R(L) a_L for the two non-annihilated pieces (R rational in L, deg <= 6)?
Also: their minimal operators, singular at x*?  Output printed (light)."""
import importlib.util, json
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location('g', ROOT / 'tests' / 'test_paper59_operator_exponents.py')
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
tw = importlib.util.spec_from_file_location('t', ROOT / 'debug' / 'op52_deriv_twopieces.py')
t = importlib.util.module_from_spec(tw); tw.loader.exec_module(t)
BAD = {'A1': ["(Fraction(-162, 25), '1', None)", "(Fraction(-9, 5), '1', None)"],
       'A4': ["(Fraction(-239, 50), '1', None)", "(Fraction(-32, 25), '1', None)"]}
for case, p in g.CP:
    e = g.embedded(case, p)
    a, b = (e['seqs'][k] for k in BAD[case])
    for d in range(0, 7):
        # b_L Q(L) - a_L P(L) = 0 : unknowns Q_0..Q_d, P_0..P_d
        cols = [[b[L] * pow(L, j, p) % p for L in range(300)] for j in range(d + 1)] + \
               [[(-a[L]) * pow(L, j, p) % p for L in range(300)] for j in range(d + 1)]
        rhs = [0] * 300
        # homogeneous: fix Q_d... use nullspace via solving with last column as rhs
        s = t.solve(cols[:-1], [(-c) % p for c in cols[-1]], p)
        if s is not None:
            print(case, p, 'b_L Q(L) = a_L P(L) with degree', d); break
    else:
        print(case, p, 'no rational-in-L ratio up to degree 6')
    for k in BAD[case]:
        c = g.polys(e['piece_ops'][k], p)
        print('  ', k, 'order', len(c) - 1, 'singular at x*:', g.singular_at(c, g.XSTAR[case], p),
              'same op as partner:', e['piece_ops'][BAD[case][0]] == e['piece_ops'][BAD[case][1]])
