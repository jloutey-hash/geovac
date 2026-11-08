"""Leg 5 -- the CORRECTED ABW target: the 2nd-order elliptic PF in the MODULUS rho.

L(D,rho) = int_1^inf e^{-Dx} dx / sqrt(Q),   Q=(x^2-1)(rho x^2 + 1-rho),  m = 1-rho.
Period (D=0): L(0,rho) = K(1-rho) = K(m)  (Leg 4, 19 digits).

Legendre PF for K(m):  m(1-m)K'' + (1-2m)K' - 1/4 K = 0.
In rho (m=1-rho, d/dm = -d/drho) the modulus operator is
    M_rho[f] = rho(1-rho) f'' - (2rho-1) f' - 1/4 f.
Homogeneous solutions: the two periods  K(1-rho)  and  K(rho).

Claims to verify numerically (all derivatives EXACT under the integral;
 d_rho Q = (x^2-1)^2, d_rho^2 Q = 0, so
   d_rho   L = -1/2 int e^{-Dx}(x^2-1)^2 Q^{-3/2} dx
   d_rho^2 L = +3/4 int e^{-Dx}(x^2-1)^4 Q^{-5/2} dx ):

 (A) M_rho annihilates the period  L(0,rho)=K(1-rho)     -> ~0.
 (B) M_rho annihilates the SECOND period K(rho)          -> ~0  (independent hom. soln).
 (C) applied to L(D,rho), D!=0, M_rho gives an explicit nonzero SOURCE S(D,rho)
     -- the ABW inhomogeneity.  This is the object variation-of-parameters against
     {K(1-rho), K(rho)} integrates into an elliptic dilogarithm.
 (D) exhibit the source's small-D behaviour (the D ln D fingerprint should live here).
"""
import mpmath as mp
mp.mp.dps = 40


def _quad_u(gu):
    """int_0^inf gu(u) du, gu already the u-space integrand (x=1+u^2, Jacobian folded
    in and the (x-1)^{-1/2}-type endpoint cancelled analytically)."""
    return mp.quad(gu, [0, 1, 3, 8, 20, mp.inf])


def L(D, rho):
    # 2u * e^{-Dx}/(u sqrt(x+1) sqrt(P)) = 2 e^{-Dx}/sqrt((x+1) P),  P=rho x^2+1-rho
    D = mp.mpf(D); rho = mp.mpf(rho)
    def gu(u):
        x = 1 + u * u; P = rho * x * x + 1 - rho
        return 2 * mp.e ** (-D * x) / mp.sqrt((x + 1) * P)
    return _quad_u(gu)


def dL_drho(D, rho):
    # -1/2 (x^2-1)^2 Q^{-3/2} = -1/2 (x^2-1)^{1/2} P^{-3/2};  *2u Jac -> -u^2 sqrt(x+1) P^{-3/2}
    D = mp.mpf(D); rho = mp.mpf(rho)
    def gu(u):
        x = 1 + u * u; P = rho * x * x + 1 - rho
        return -u * u * mp.e ** (-D * x) * mp.sqrt(x + 1) / P ** mp.mpf('1.5')
    return _quad_u(gu)


def d2L_drho2(D, rho):
    # 3/4 (x^2-1)^4 Q^{-5/2} = 3/4 (x^2-1)^{3/2} P^{-5/2}; *2u Jac -> 3/2 u^4 (x+1)^{3/2} P^{-5/2}
    D = mp.mpf(D); rho = mp.mpf(rho)
    def gu(u):
        x = 1 + u * u; P = rho * x * x + 1 - rho
        return (mp.mpf(3) / 2) * u ** 4 * mp.e ** (-D * x) * (x + 1) ** mp.mpf('1.5') / P ** mp.mpf('2.5')
    return _quad_u(gu)


def M_rho_on_L(D, rho):
    rho = mp.mpf(rho)
    f, f1, f2 = L(D, rho), dL_drho(D, rho), d2L_drho2(D, rho)
    return rho * (1 - rho) * f2 - (2 * rho - 1) * f1 - mp.mpf(1) / 4 * f


def M_rho_on(fn, rho, h=mp.mpf('1e-12')):
    """Apply M_rho to a closed-form f(rho) via central finite differences."""
    rho = mp.mpf(rho)
    f = fn(rho)
    f1 = (fn(rho + h) - fn(rho - h)) / (2 * h)
    f2 = (fn(rho + h) - 2 * f + fn(rho - h)) / h ** 2
    return rho * (1 - rho) * f2 - (2 * rho - 1) * f1 - mp.mpf(1) / 4 * f


def main():
    RHO = mp.mpf('0.37')
    print(f"rho = {RHO}\n")

    print("(A) M_rho on the period L(0,rho) = K(1-rho):")
    print(f"    M_rho[L(0,.)]  = {mp.nstr(M_rho_on_L(mp.mpf('0'), RHO), 6)}")
    print(f"    M_rho[K(1-.)]  = {mp.nstr(M_rho_on(lambda r: mp.ellipk(1 - r), RHO), 6)}  (closed-form period)")

    print("\n(B) M_rho on the SECOND period K(rho):")
    print(f"    M_rho[K(.)]    = {mp.nstr(M_rho_on(lambda r: mp.ellipk(r), RHO), 6)}")
    print("    both ~0 => {K(1-rho), K(rho)} are the two homogeneous solutions.\n")

    print("(C) M_rho on L(D,rho), D!=0 -> the ABW SOURCE S(D,rho):")
    for D in ['0.5', '1.0', '2.0']:
        S = M_rho_on_L(mp.mpf(D), RHO)
        print(f"    D={D:>4}:  S = {mp.nstr(S, 8)}")
    print("    nonzero => L(D,.) is the INHOMOGENEOUS solution of the elliptic PF.\n")

    print("(D) small-D structure of the source (D ln D fingerprint):")
    print("     D          S(D)          S(D)/D        S(D)/(D ln D)")
    for D in ['0.2', '0.1', '0.05', '0.025', '0.0125']:
        Dm = mp.mpf(D)
        S = M_rho_on_L(Dm, RHO)
        print(f"   {D:>7}   {mp.nstr(S,8):>12}   {mp.nstr(S/Dm,8):>12}   {mp.nstr(S/(Dm*mp.log(Dm)),8):>12}")
    print("    if S/(D ln D) -> const, the source carries the D ln D non-analyticity")
    print("    (the inhomogeneous 'regulator' term that becomes the elliptic dilog).")


if __name__ == '__main__':
    main()
