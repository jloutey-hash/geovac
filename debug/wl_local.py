"""WL derivation: the leading bond-angle singular term of (XY|XZ) at unequal bond lengths, from the Landau pinch.

The cusp-cusp singular part is S ~ Phi |w|^7, w = Y - Z, near |w|^2 = 0, i.e. x* = (D1^2 + D2^2)/(2 D1 D2).
At a null w every Taylor degree of the cusp functions G1, G2 contributes at the same order; the degree-(d1, d2)
coefficient is (pi^2/630) (-1)^d1 (d1+1)! (d2+1)! 24/(d1+d2+4)! (checked against debug/wa6_local.py for five
monomials), which are the moments of 24 s u on the simplex s, u >= 0, s + u <= 1. So
    Phi = (pi^2/630) 24 int int s u G1(-s w) G2(u w) ds du,
the cusp functions averaged over the segment between the two H atoms (the collinear pinch). On a null w,
|Y - s w|^2 = (1-s) D1^2 + s D2^2, and sinh(zH|a|)/|a| = zH there, so
    Phi = (pi^2/630) 24 zH^2 NH^2 NX^2 int int s u exp(-zX (R1(s) + R2(u))),
    R1(s) = sqrt((1-s) D1^2 + s D2^2), R2(u) = sqrt((1-u) D2^2 + u D1^2).
At D1 = D2 this is WA6's C. Leading term: L_l = Legendre coefficients of Phi (2 D1 D2)^(7/2) (x* - x)^(7/2).
Usage: python debug/wl_local.py
"""
import mpmath as mp

mp.mp.dps = 40


def phi(zX, zH, D1, D2):
    zX, zH, D1, D2 = map(mp.mpf, (zX, zH, D1, D2))
    N2 = (zH ** 3 / mp.pi) * (zX ** 3 / mp.pi)
    R1 = lambda s: mp.sqrt((1 - s) * D1 ** 2 + s * D2 ** 2)
    R2 = lambda u: mp.sqrt((1 - u) * D2 ** 2 + u * D1 ** 2)
    I = mp.quad(lambda s: s * mp.exp(-zX * R1(s)) * mp.quad(lambda u: u * mp.exp(-zX * R2(u)), [0, 1 - s]), [0, 1])
    return mp.pi ** 2 / 630 * 24 * zH ** 2 * N2 * I


def naive(zX, zH, D1, D2):
    zX, zH, D1, D2 = map(mp.mpf, (zX, zH, D1, D2))
    return zH ** 5 * zX ** 3 * mp.exp(-zX * (D1 + D2)) / 630


def leg_coeffs(A, D1, D2, ls):
    D1, D2 = mp.mpf(D1), mp.mpf(D2)
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    f = lambda x: A * (2 * D1 * D2 * (xs - x)) ** mp.mpf(3.5)
    return {l: (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: f(x) * mp.legendre(l, x), mp.linspace(-1, 1, 2 * l // 3 + 3))
            for l in ls}


CASES = {'E1 2.2': (1, 1, 1.809, 2.2), 'E1 1.4': (1, 1, 1.809, 1.4),
         'E2 2.2': (2.25, 1.24, 1.809, 2.2), 'E2 1.4': (2.25, 1.24, 1.809, 1.4)}

if __name__ == "__main__":
    print("symmetric check E1:", mp.nstr(phi(1, 1, 1.809, 1.809) / naive(1, 1, 1.809, 1.809), 20))
    for k, a in CASES.items():
        p, n = phi(*a), naive(*a)
        print(k, "Phi =", mp.nstr(p, 16), " Phi/naive =", mp.nstr(p / n, 10))
