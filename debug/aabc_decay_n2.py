import sys
sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
import mpmath as mp
import na_exact as E
D1, D2, zX, zH = '1.5', '2.4', '1.0', '1.2'
Lmax = 240
Q = E.coefficients(D1, D2, zX, zH, Lmax, E.Field(None))
for dps in (2000, 2600):
    mp.mp.dps = dps
    rho = mp.mpf(5) / 8
    F = [sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items()) for q in Q]
    G = [F[L] / rho ** L for L in range(Lmax + 1)]
    beta = {L: round(float(mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1))), 4)
            for L in range(40, Lmax, 60)}
    print(dps, beta, sorted(set(int(mp.sign(f)) for f in F[5:])), flush=True)
