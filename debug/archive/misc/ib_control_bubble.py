"""IB record control: the D = 3 one-loop bubble with masses (m1, m2) through ib_reduce.
Reduce dB/dz to masters at rational z, rebuild the rational coefficients, and check the DE on the closed
form B = (2 pi^2/sqrt z) arctan(sqrt z/(m1+m2)), tadpoles T(m) = -2 pi^2 m (dimensional regularization)."""
from fractions import Fraction as F
import mpmath as mp
from ib_reduce import Family, Reducer, shifts_box, dz_terms, ratrec, fr

mp.mp.dps = 40
m1, m2 = F(1), F(3, 2)
fam = Family(["l"], ["p"], [({"l": 1}, m1 * m1, F(1)), ({"l": 1, "p": 1}, m2 * m2, F(1))],
             zero_sector=lambda pos: not any(pos))
worst = 0
for z in (F(3), F(7, 2), F(11)):
    fam.setup({("p", "p"): z}, 3)
    R = Reducer(fam)
    R.seed(shifts_box(2, 3, 2, fam.integer))
    red = R.reduce(dz_terms(fam, (0, 0), "p", fr(z)))
    masters = sorted(R.reduce({(0, 0): 1}).keys())
    coef = {n: ratrec(c) for n, c in red.items()}
    sz = mp.sqrt(mp.mpf(z.numerator) / z.denominator)
    M = mp.mpf(m1.numerator) / m1.denominator + mp.mpf(m2.numerator) / m2.denominator
    B = lambda zz: 2 * mp.pi**2 / mp.sqrt(zz) * mp.atan(mp.sqrt(zz) / M)
    zf = sz**2
    val = {(0, 0): B(zf), (1, -1): None}
    T = {(0, -1): -2 * mp.pi**2 * mp.mpf(m1.numerator) / m1.denominator,
         (-1, 0): -2 * mp.pi**2 * mp.mpf(m2.numerator) / m2.denominator}
    lhs = mp.diff(B, zf)
    rhs = 0
    for n, c in coef.items():
        v = B(zf) if n == (0, 0) else T[n]
        rhs += mp.mpf(c.numerator) / c.denominator * v
    rel = abs(lhs - rhs) / abs(lhs)
    worst = max(worst, rel)
    print("z", z, "masters", masters, "coeffs", {n: str(c) for n, c in coef.items()}, "rel", mp.nstr(rel, 3))
print("CONTROL bubble worst rel", mp.nstr(worst, 3), "PASS" if worst < 1e-30 else "FAIL")
