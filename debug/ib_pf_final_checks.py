"""IB4 final checks on the lifted PF operator: agreement with every prime's operator, and a 1/200 change in
one coefficient (the largest-degree coefficient of the order-5 polynomial) must fail on the series."""
import json, glob
from fractions import Fraction
import mpmath as mp
from ib_crt import series_residual
d = json.load(open("data/ib_pf_op_Q.json"))
ip = [[int(v) for v in p] for p in d["polys"]]
agree = n = 0
for f in sorted(glob.glob("data/ib_pf_p*_op.json")):
    o = json.load(open(f)); p = o["P"]; n += 1
    il = pow(ip[-1][-1] % p, -1, p)
    agree += [[v * il % p for v in q] for q in ip] == [[v % p for v in q] for q in o["polys"]]
print("agrees with the operator at", agree, "of", n, "primes")
ser = json.load(open("data/pf_series_direct_explore_N300_d280_460x460_sub.json"))["phi_n"]
op = [[Fraction(v) for v in p] for p in ip]
op[5][-1] = op[5][-1] * (1 + Fraction(1, 200))
print("1/200 change, residual:", mp.nstr(series_residual(op, ser, dps=300), 3))
