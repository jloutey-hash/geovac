"""Decisive saturation test: exact-closure (flat vs degree) vs approximation (monotone).
Hyp A: S = {L,L',L'',L'''} x poly(D)                    -- expect SATURATE at deg2 (in-module)
Hyp C2: S = {L,L'} x poly + {K0(D),K1(D),J0(wD),Y0(wD)} x poly  (both Bessel sectors)
        -- if it also saturates at <=12 params, an elementary inhomogeneity EXISTS.
"""
import mpmath as mp
mp.mp.dps = 45
RHO = mp.mpf('0.37')
W = mp.sqrt((1 - RHO) / RHO)


def _q(gu):
    return mp.quad(gu, [0, 1, 3, 8, 20, mp.inf])


def data_at(D):
    D = mp.mpf(D); rho = RHO
    def mom(n):
        return _q(lambda u: (1 + u * u) ** n * 2 * mp.e ** (-D * (1 + u * u)) / mp.sqrt((2 + u * u) * (rho * (1 + u * u) ** 2 + 1 - rho)))
    m0, m1, m2, m3 = mom(0), mom(1), mom(2), mom(3)
    Lk = [m0, -m1, m2, -m3]
    dLdr = _q(lambda u: -u * u * mp.e ** (-D * (1 + u * u)) * mp.sqrt(2 + u * u) / (rho * (1 + u * u) ** 2 + 1 - rho) ** mp.mpf('1.5'))
    d2 = _q(lambda u: (mp.mpf(3) / 2) * u ** 4 * mp.e ** (-D * (1 + u * u)) * (2 + u * u) ** mp.mpf('1.5') / (rho * (1 + u * u) ** 2 + 1 - rho) ** mp.mpf('2.5'))
    S = rho * (1 - rho) * d2 - (2 * rho - 1) * dLdr - mp.mpf(1) / 4 * Lk[0]
    bess = [mp.besselk(0, D), mp.besselk(1, D), mp.besselj(0, W * D), mp.bessely(0, W * D)]
    return Lk, bess, S


def fit(kind, deg, verify=8):
    if kind == 'A':
        gens = lambda Lk, bs: Lk                 # 4 generators
    else:
        gens = lambda Lk, bs: [Lk[0], Lk[1]] + bs  # 2 + 4 = 6 generators
    # probe generator count
    ng = len(gens([0, 0, 0, 0], [0, 0, 0, 0]))
    ncols = ng * (deg + 1)
    Dall = [mp.mpf('0.4') + mp.mpf('0.31') * i for i in range(ncols + verify)]
    dat = [data_at(D) for D in Dall]

    def row(D, Lk, bs):
        r = []
        for ggen in gens(Lk, bs):
            for k in range(deg + 1):
                r.append(D ** k * ggen)
        return r

    A = mp.matrix([row(Dall[i], dat[i][0], dat[i][1]) for i in range(ncols)])
    b = mp.matrix([dat[i][2] for i in range(ncols)])
    coef = mp.lu_solve(A, b)
    mr = mp.mpf(0)
    for i in range(ncols, ncols + verify):
        pred = sum(coef[j] * row(Dall[i], dat[i][0], dat[i][1])[j] for j in range(ncols))
        mr = max(mr, abs(pred - dat[i][2]))
    return ncols, mr


print(f"rho={RHO}, dps=45, floor ~1e-40\n")
print("Hyp A  {L,L',L'',L'''} x poly(D):")
for deg in [1, 2, 3, 4]:
    nc, r = fit('A', deg)
    print(f"   deg {deg} ({nc:2d} params): {mp.nstr(r,4)}")
print("\nHyp C2 {L,L'} + {K0(D),K1(D),J0(wD),Y0(wD)}, x poly(D):")
for deg in [1, 2, 3, 4, 5]:
    nc, r = fit('C', deg)
    print(f"   deg {deg} ({nc:2d} params): {mp.nstr(r,4)}")
print("\nRead: FLAT across degree = exact finite closure; MONOTONE-decreasing = approximation.")
