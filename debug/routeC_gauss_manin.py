"""Route C -- the Gauss-Manin connection for the one-mass slice (ABW sprint step 1).

L(D,rho) = int_1^inf e^{-Dx} dx / sqrt(Q),   Q=(x^2-1)(rho x^2+1-rho).
Verified elsewhere: 4th-order PF ODE in D; 2nd-order modulus PF
  M_rho = rho(1-rho) d^2/drho^2 + (1-2rho) d/drho - 1/4,  periods {K(rho),K(1-rho)};
  M_rho[L] = S(D,rho) != 0.

Because L is holonomic in (D,rho) and satisfies a 4th-order D-ODE, the FLAT
Gauss-Manin connection expresses the rho-derivatives via the D-derivative basis
{L, L', L'', L'''}.  So the source should close as
    S(D,rho) = sum_{k=0}^{3} P_k(D,rho) L^{(k)}(D)      (P_k polynomial in D)
with NO separate boundary term (the L^{(k)} already carry the x=1 branch data).
We test that hypothesis (clean D-module), and a variant with an e^{-D} boundary,
by exact square solve + held-out verification.

L^{(k)}(D) = (-1)^k int x^k e^{-Dx} Q^{-1/2} dx = (-1)^k moment(k).
"""
from __future__ import annotations
import mpmath as mp
mp.mp.dps = 35


def _quad_u(gu):
    return mp.quad(gu, [0, 1, 3, 8, 20, mp.inf])


def moments_and_source(D, rho):
    """Return (L, L', L'', L''', S) at (D,rho); all from x=1+u^2 quadratures."""
    D = mp.mpf(D); rho = mp.mpf(rho)
    def base(u):
        x = 1 + u * u; P = rho * x * x + 1 - rho
        return x, P, 2 * mp.e ** (-D * x) / mp.sqrt((x + 1) * P)   # 2 e^{-Dx}/sqrt((x+1)P)
    def mom(n):
        return _quad_u(lambda u: (1 + u * u) ** n * base(u)[2])
    m0, m1, m2, m3 = mom(0), mom(1), mom(2), mom(3)
    Lk = [m0, -m1, m2, -m3]                       # L, L', L'', L'''
    # source S = M_rho[L]:  d_rho L = -1/2 int (x^2-1)^{1/2} P^{-3/2} e ; d2 = 3/2 ...
    dLdr = _quad_u(lambda u: -u * u * mp.e ** (-D * (1 + u * u)) * mp.sqrt(2 + u * u)
                   / (rho * (1 + u * u) ** 2 + 1 - rho) ** mp.mpf('1.5'))
    d2Ldr2 = _quad_u(lambda u: (mp.mpf(3) / 2) * u ** 4 * mp.e ** (-D * (1 + u * u))
                     * (2 + u * u) ** mp.mpf('1.5')
                     / (rho * (1 + u * u) ** 2 + 1 - rho) ** mp.mpf('2.5'))
    S = rho * (1 - rho) * d2Ldr2 - (2 * rho - 1) * dLdr - mp.mpf(1) / 4 * Lk[0]
    return Lk, S


def fit(rho, deg=3, basis='Dmodule', verify=6):
    """basis:
       'Dmodule' -> {L,L',L'',L'''} x poly(D)             (Hyp A: pure D-module)
       'ABW'     -> {L,L'} x poly(D) + {K0(D),K1(D)} x poly(D)  (Hyp C: period + Bessel subtopology)
    """
    rho = mp.mpf(rho)
    if basis == 'Dmodule':
        ncols = 4 * (deg + 1)
    elif basis == 'ABW':
        ncols = 2 * (deg + 1) + 2 * (deg + 1)
    Dall = [mp.mpf('0.4') + mp.mpf('0.33') * i for i in range(ncols + verify)]
    data = [moments_and_source(D, rho) for D in Dall]

    def row(D, Lk):
        r = []
        if basis == 'Dmodule':
            for Lval in Lk:
                for k in range(deg + 1):
                    r.append(D ** k * Lval)
        elif basis == 'ABW':
            for Lval in (Lk[0], Lk[1]):                     # L, L' (period level)
                for k in range(deg + 1):
                    r.append(D ** k * Lval)
            for Bs in (mp.besselk(0, D), mp.besselk(1, D)):  # Bessel subtopology
                for k in range(deg + 1):
                    r.append(D ** k * Bs)
        return r

    A = mp.matrix([row(Dall[i], data[i][0]) for i in range(ncols)])
    b = mp.matrix([data[i][1] for i in range(ncols)])
    coef = mp.lu_solve(A, b)
    maxres = mp.mpf(0)
    for i in range(ncols, ncols + verify):
        pred = sum(coef[j] * row(Dall[i], data[i][0])[j] for j in range(ncols))
        maxres = max(maxres, abs(pred - data[i][1]))
    return coef, maxres


def _tag(res):
    return "CLOSES" if res < mp.mpf(10) ** (-18) else ("partial" if res < mp.mpf(10) ** (-6) else "NO")


def main():
    print("Route C -- Gauss-Manin connection / ABW-structure discriminator (one-mass slice)\n")
    print("Hyp A (in-module):  S = sum_{k<=3} P_k(D) L^{(k)}       [circular for VoP if it closes]")
    for deg in [1, 2, 3]:
        _, res = fit('0.37', deg=deg, basis='Dmodule')
        print(f"   deg {deg}: max held-out residual = {mp.nstr(res,4):>12}   [{_tag(res)}]")

    print("\nHyp C (ABW inhomogeneity):  S = {L,L'} poly(D)  +  {K0(D),K1(D)} poly(D)")
    print("   closes  => genuine period + elementary-Bessel subtopology => ELi reachable")
    print("   fails   => no elementary inhomogeneity => object is a new elliptic period")
    for rho in ['0.37', '0.5', '0.23']:
        line = []
        for deg in [2, 3, 4]:
            _, res = fit(rho, deg=deg, basis='ABW')
            line.append(f"deg{deg}:{mp.nstr(res,3)}[{_tag(res)}]")
        print(f"   rho={rho:>5}:  " + "   ".join(line))


if __name__ == '__main__':
    main()
