"""Route C -- ABW-route diagnostic on the 4th-order Picard-Fuchs ODE.

Object (Paper 59 eq:laplace, memo achieved-item 2/4):
    L(D) = int_1^inf e^{-Dx} dx / sqrt(Q(x)),   Q(x) = (x^2-1)(rho x^2 + 1 - rho)
= sqrt(c1) * N(D), the Laplace transform of the holomorphic elliptic differential
omega = dx/sqrt(Q) over the cycle [1, inf).  rho = c2/c1 in (0,1).

Memo's 4th-order PF ODE in D (L = sqrt(c1) N):
    D*rho*L'''' + 2*rho*L''' + D*(1-2rho)*L'' + (1-2rho)*L' - D*(1-rho)*L = 0
char poly rho(s^2-1)(s^2+(1-rho)/rho): rates s = +-1 and +- i*omega, omega=sqrt((1-rho)/rho).

KEY: derivatives are EXACT moment quadratures (differentiate under the integral):
    L^{(n)}(D) = (-1)^n int_1^inf x^n e^{-Dx}/sqrt(Q) dx.
No numerical differentiation of a quadrature-defined function.

Legs:
  1. ODE residual (state verification -- guard vs stale memo).
  2. Single-Bessel test: does {K0(D), I0(D), J0(wD), Y0(wD)} solve the D-operator?
     (memo says the operator does NOT factor into elementary Bessel ops; verify.)
  3. Are the ELLIPTIC PERIODS D-solutions?  A period is a constant in D; O[const] =
     -D(1-rho)*const != 0.  So periods live in the MODULUS (rho) variable, not D.
  4. Limit/period checks: L(0) = complete elliptic K;  rho->0 gives L(D)->K0(D).
"""
import mpmath as mp
mp.mp.dps = 40

RHO = mp.mpf('0.37')      # generic modulus
OMEGA = mp.sqrt((1 - RHO) / RHO)


def Qx(x, rho=RHO):
    return (x * x - 1) * (rho * x * x + 1 - rho)


def moment(n, D, rho=RHO):
    """(-1)^n L^{(n)}(D) = int_1^inf x^n e^{-Dx}/sqrt(Q) dx, via x = 1+u^2 (kills the
    1/sqrt(x-1) endpoint)."""
    D = mp.mpf(D)
    def f(u):
        x = 1 + u * u
        # dx = 2u du, sqrt(x-1)=u cancels -> 2 e^{-Dx} x^n / sqrt((x+1)(rho x^2+1-rho))
        return 2 * mp.e ** (-D * x) * x ** n / mp.sqrt((x + 1) * (rho * x * x + 1 - rho))
    return mp.quad(f, [0, 1, 3, 8, 20, mp.inf])


def L(D, rho=RHO):
    return moment(0, D, rho)


def Ln(n, D, rho=RHO):
    return (-1) ** n * moment(n, D, rho)


def ode_residual(D, rho=RHO):
    L0 = Ln(0, D, rho); L1 = Ln(1, D, rho); L2 = Ln(2, D, rho)
    L3 = Ln(3, D, rho); L4 = Ln(4, D, rho)
    return (D * rho * L4 + 2 * rho * L3 + D * (1 - 2 * rho) * L2
            + (1 - 2 * rho) * L1 - D * (1 - rho) * L0)


def apply_op(y, D, rho=RHO):
    """Apply the D-operator to an arbitrary smooth y(D) via mp.diff (closed-form y)."""
    y0 = y(D); y1 = mp.diff(y, D, 1); y2 = mp.diff(y, D, 2)
    y3 = mp.diff(y, D, 3); y4 = mp.diff(y, D, 4)
    return (D * rho * y4 + 2 * rho * y3 + D * (1 - 2 * rho) * y2
            + (1 - 2 * rho) * y1 - D * (1 - rho) * y0)


def main():
    print(f"rho = {RHO}, omega = sqrt((1-rho)/rho) = {mp.nstr(OMEGA, 12)}\n")

    print("LEG 1 -- 4th-order PF ODE residual (exact moment quadratures):")
    for D in ['0.5', '1.0', '2.0', '3.5']:
        r = ode_residual(mp.mpf(D))
        print(f"   D={D:>4}:  residual = {mp.nstr(r, 6)}")
    print("   (expect ~0 to full precision -> confirms current state of the ODE)\n")

    print("LEG 2 -- single-Bessel test: apply the D-operator to each single Bessel")
    w = OMEGA
    bessels = {
        'K0(D)':      lambda D: mp.besselk(0, D),
        'I0(D)':      lambda D: mp.besseli(0, D),
        'J0(w D)':    lambda D: mp.besselj(0, w * D),
        'Y0(w D)':    lambda D: mp.bessely(0, w * D),
    }
    for name, y in bessels.items():
        res = apply_op(y, mp.mpf('1.7'))
        val = y(mp.mpf('1.7'))
        print(f"   O[{name:8}](1.7) = {mp.nstr(res, 6):>16}   (fn value {mp.nstr(val,6)})")
    print("   nonzero => single Bessels are NOT solutions => the K0/I0 and J0/Y0")
    print("   sectors are COUPLED (genuine ellipticity), matching the memo.\n")

    print("LEG 3 -- are the elliptic PERIODS solutions of the D-operator?")
    const = lambda D: mp.mpf(1)
    print(f"   O[const=1](1.7) = {mp.nstr(apply_op(const, mp.mpf('1.7')), 6)}")
    print("   = -D(1-rho) != 0.  A period is D-constant, so periods are NOT D-solutions;")
    print("   they solve the MODULUS (rho) Picard-Fuchs, not this D-ODE.\n")

    print("LEG 4 -- period / limit sanity:")
    # L(0) = complete elliptic integral.  Q=(x^2-1)(rho x^2+1-rho); with c1=1 the memo
    # period formula: int_0^inf dk/sqrt((c1 k^2+1)(c2 k^2+1)); here in x-vars L(0)=int_1^inf dx/sqrt(Q).
    L0 = L(mp.mpf('0'))
    # closed period: int_1^inf dx/sqrt((x^2-1)(rho x^2+1-rho)).  Standard reduction:
    # = (1/sqrt(1)) K(m) form.  Compare against direct high-precision quad + ellipk.
    m_param = 1 - RHO           # per memo modulus m = 1 - c_min/c_max with c1=1,c2=rho (rho<1)
    Kell = mp.ellipk(m_param)
    print(f"   L(0)              = {mp.nstr(L0, 20)}")
    print(f"   ellipk(1-rho)     = {mp.nstr(Kell, 20)}   ratio = {mp.nstr(L0/Kell, 12)}")
    # rho -> 0 : L(D) -> K0(D)
    for rr in ['0.1', '0.01', '0.001']:
        Ld = L(mp.mpf('1.3'), rho=mp.mpf(rr))
        print(f"   rho={rr:>6}: L(1.3)={mp.nstr(Ld,12)}   K0(1.3)={mp.nstr(mp.besselk(0,mp.mpf('1.3')),12)}")
    print("   (L(D) -> K0(D) as rho->0: coincident Fock scales, genus-0 limit)")


if __name__ == '__main__':
    main()
