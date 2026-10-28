"""WG derivation part 3b (phase 1; no scored F_l): q_l to order l^-5 on the symmetric line T0 = T1 = T.
R_l expanded directly in eps = 1/l term by term (series of each factor), then q = R(l+1)/R(l)."""
import sympy as sp
e, T = sp.symbols('epsilon T', positive=True)
N = 7
def ser(x):
    return sp.series(x, e, 0, N).removeO()
L = 1 / e
R = 0
for k in range(N):
    b = L + 2 * k + sp.Rational(3, 2)
    M = 0
    for j in range(N):
        M += T ** j * ser(1 / sp.rf(b, j)) if j else 1
    dbl = 1
    for i in range(1, k + 1):
        dbl = ser(dbl * ser(1 / (2 * L + 2 * i + 1)))
    pre = ser(dbl * ser((2 * L + 1) / (2 * L + 4 * k + 1)))
    R += sp.expand(ser(T ** (2 * k) / (2 ** k * sp.factorial(k)) * pre * M))
R = sp.expand(R)
Rc = [sp.expand(R.coeff(e, j)) for j in range(N)]
# R(l+1): eps -> eps/(1+eps)
Rp = sp.expand(ser(sum(c * (e / (1 + e)) ** j for j, c in enumerate(Rc))))
q = sp.expand(ser(Rp * ser(1 / sum(c * e ** j for j, c in enumerate(Rc)))))
qc = [sp.factor(q.coeff(e, j)) for j in range(6)]
for j, c in enumerate(qc):
    print(j, c, [float(c.subs(T, v)) for v in (1.6362405, 5.69127130435)])

# toy check, symmetric T = 1.2 (not a scored case), against the exact series form
import mpmath as mp
import wg_derive as W
mp.mp.dps = 50
Tt = mp.mpf('1.2')
for l in (20, 40, 80):
    qa = W.F_closed(Tt, Tt, l + 1) * (2 * l + 3) / (W.F_closed(Tt, Tt, l) * Tt)
    qs4 = sum(mp.mpf(float(c.subs(T, 1.2))) / mp.mpf(l) ** j for j, c in enumerate(qc[:5]))
    qs5 = qs4 + mp.mpf(float(qc[5].subs(T, 1.2))) / mp.mpf(l) ** 5
    print(f"toy T=1.2 l={l}: resid after l^-4 {mp.nstr(qa - qs4, 5)}, after l^-5 {mp.nstr(qa - qs5, 5)}")
