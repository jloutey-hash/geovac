"""exch_reverse_lemma.py -- numerical check of Step 1 of the independence lemma for a > b: the log coefficients A_i
of the nine functions at R = 0 (g_i = A_i ln R + entire) with the real branch (LC = ln((a+b)/|b-a|), PV Ei).
Fifth forward difference on e, 2e, .., 6e kills analytic functions to O(e^5); for (anything) x ln R it leaves an
O(1)..O(e^4) remainder.  Reported: |D5(g_i - A_i ln)| (must be O(e^5)) and |D5(g_i)| (size of what is removed).
Run: python debug/exch_reverse_lemma.py -> debug/data/exch_reverse_lemma.json"""
import json
import sys

import mpmath as mp

sys.path.insert(0, "tests")
import _exact_order as EO  # noqa: E402

mp.mp.dps = 50
out = {}
for (a, b) in ((2, 1), (mp.mpf(3) / 2, mp.mpf(4) / 5), (3, 1), (mp.sqrt(5), mp.mpf(1))):
    a, b = mp.mpf(a), mp.mpf(b)
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / abs(b - a))
    g = [
        lambda R: mp.ei(-(a + b) * R),
        lambda R: mp.exp(-(a + b) * R),
        lambda R: mp.exp(-2 * a * R),
        lambda R: mp.exp(-2 * b * R),
        lambda R: mp.exp(2 * a * R) * mp.ei(-2 * (a + b) * R),
        lambda R: mp.exp(2 * b * R) * mp.ei(-2 * (a + b) * R),
        lambda R: mp.exp(-2 * a * R) * (mp.ei(2 * (a - b) * R) - 2 * mp.ei((a - b) * R) - LC),
        lambda R: mp.exp(-2 * b * R) * (mp.ei(2 * (b - a) * R) - 2 * mp.ei((b - a) * R) - LC),
    ]
    fm1 = lambda R: c * (g[5](R) - g[4](R) - g[7](R) + g[6](R))  # noqa: E731
    A, A1, L = EO.log_coefficients(a, b)
    e = mp.mpf("1e-3")
    rows = []
    for i in range(8):
        r1 = EO.fifth_difference(lambda R: g[i](R) - A[i](R) * mp.log(R), e)
        r0 = EO.fifth_difference(g[i], e)
        assert mp.im(r1) == 0
        rows.append({"g": i + 2, "D5(g - A ln)": mp.nstr(abs(r1), 3), "D5(g)": mp.nstr(abs(r0), 3)})
    d0 = EO.f0_fifth_difference(fm1, e)
    d1 = d0 - EO.fifth_difference(lambda R: A1(R) * mp.log(R), e)
    rows.append({"g": 1, "D5(g - A ln)": mp.nstr(abs(d1), 3), "D5(g)": mp.nstr(abs(d0), 3)})
    tag = "(%s, %s)" % (mp.nstr(a, 8), mp.nstr(b, 8))
    out[tag] = rows
    ok = all(mp.mpf(r["D5(g - A ln)"]) < mp.mpf("1e-10") for r in rows)  # analytic floor ~ (2 max(a,b) e)^5
    print(tag, "all residuals at the analytic O(e^5) floor (< 1e-10; removed ln parts 1e-7..0.09):", ok, flush=True)
    for r in rows:
        print("   ", r, flush=True)
json.dump(out, open("debug/data/exch_reverse_lemma.json", "w"), indent=1)
