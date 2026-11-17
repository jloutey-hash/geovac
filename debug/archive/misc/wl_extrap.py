"""WL extrapolation rule and its synthetic control (record docs/predictions/2026-10-28_water_landau_unequal.md).
r_l = F_l / L_l = 1 + c1/l + c2/l^2 + ...; r_inf = the degree-3 polynomial in 1/l through r at l = 16, 24, 32, 40,
evaluated at 1/l = 0 (Neville). Control: f = (x*-x)^(7/2) (1 + beta (x*-x) + gamma (x*-x)^2), amplitude 1 known."""
import mpmath as mp

mp.mp.dps = 40
LS = (16, 24, 32, 40)


def r_inf(r):
    """Neville extrapolation to h = 1/l = 0 through the points (1/l, r_l), l in LS."""
    h = [mp.mpf(1) / l for l in LS]
    P = [mp.mpf(r[l]) for l in LS]
    n = len(h)
    for m in range(1, n):
        for i in range(n - m):
            P[i] = (h[i] * P[i + 1] - h[i + m] * P[i]) / (h[i] - h[i + m])
    return P[0]


def leg(f, l):
    return (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: f(x) * mp.legendre(l, x), mp.linspace(-1, 1, 2 * l // 3 + 3))


def control(xs, beta, gamma):
    xs = mp.mpf(xs)
    base = lambda x: (xs - x) ** mp.mpf(3.5)
    full = lambda x: base(x) * (1 + beta * (xs - x) + gamma * (xs - x) ** 2)
    r = {l: leg(full, l) / leg(base, l) for l in LS}
    return r, r_inf(r)


if __name__ == "__main__":
    for D2 in (1.4, 2.2, 2.6):
        xs = (1.809 ** 2 + D2 ** 2) / (2 * 1.809 * D2)
        for beta, gamma in ((-3, 0), (3, 0), (-3, 4), (5, -6)):
            r, ri = control(xs, beta, gamma)
            print(f"D2={D2} x*={xs:.4f} beta={beta} gamma={gamma}: r16={float(r[16]):.4f} r40={float(r[40]):.4f}"
                  f"  r_inf={float(ri):.6f}  true={1 + beta * 0 + gamma * 0}")
