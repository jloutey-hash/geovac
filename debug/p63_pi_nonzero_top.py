"""Top-term constant for non-homogeneous trees by independent exterior-tail quadrature (own code; no model).
Hat normalisation: channel-L exterior tail E_L = (1/(2L+1))[-r^-L-1 A + r^L B], A = int_r^oo s^(L+2) W, B = int_r^oo s^(1-L) W.
Everything scaled by e^{kx}.  Predicted C = (-1)^(n-1)(pi/2)^((n-1)/2) Pi / prod k_u^2, exponent (5-n)/2 = 1."""
import math

import numpy as np
from scipy.special import kve, roots_laguerre

VN, VW = roots_laguerre(90)


def lag(f, x, k):
    v = VN / k
    return float(np.sum(VW / k * f(x + v)))


def tail(Wsc, L, x, k):
    """(E, E') scaled by e^{kx}, for a source W = e^{-ks} Wsc(s)."""
    A = lag(lambda s: s ** (L + 2) * Wsc(s), x, k)
    B = lag(lambda s: s ** (1 - L) * Wsc(s), x, k)
    E = (-x ** (-L - 1) * A + x ** L * B) / (2 * L + 1)
    dE = ((L + 1) * x ** (-L - 2) * A + L * x ** (L - 1) * B) / (2 * L + 1)   # d/dx of the unscaled tail, rescaled
    return E, dE


def vec(fn):
    return lambda s: np.array([fn(si) for si in np.atleast_1d(s)])


def Es(x):   # s leaf: channel 0, W = K0
    return tail(lambda s: kve(0, s), 0, x, 1.0)


def E1(x):   # p_z leaf: d_z K0 = -mu K1, channel 1, W = -K1
    return tail(lambda s: -kve(1, s), 1, x, 1.0)


def G_N2(x):
    """root 1 + d_z (P_0(-d) = 1 - d_z), leaves 1 + d_z and s.  G = x^2 [Phi0 Phis - (1/3)(g' + 2g/x)], g = Phi1 Phis."""
    e0, d0 = Es(x)
    e1, d1 = E1(x)
    g, dg = e1 * e0, d1 * e0 + e1 * d0
    return x ** 2 * (e0 * e0 - (dg + 2 * g / x) / 3)


def G_N3(x):
    """root 1 + d_z^2 (P_0(-d) = 1 + d_z^2), two s leaves.  G = x^2 [Psi + (1/3)(Psi'' + 2 Psi'/x)], Psi = Phis^2;
    Phis'' = -W - 2 Phis'/x (channel-0 ODE), W = K0."""
    e, de = Es(x)
    w = kve(0, x)
    dde = -w - 2 * de / x
    Psi, dPsi, ddPsi = e * e, 2 * e * de, 2 * de * de + 2 * e * dde
    return x ** 2 * (Psi + (ddPsi + 2 * dPsi / x) / 3)


def G_N6(x):
    """root 1 + d_z + d_z^2 (P_0(-d) = 1 - d_z + d_z^2), middle d_z + lap/2, leaf p_z (path, root at end).
    middle sources: W0 = -K1 Phi1/3, W1 = (1/2)(K0 - K1/r) Phi1, W2 = -(2/3) K1 Phi1.
    G = x^2 [M0 - (1/3)(M1' + 2M1/x) + (1/3)(M0'' + 2M0'/x) + (2/15)(g' + 2g/x)], g = M2' + 3 M2/x."""
    W0 = vec(lambda s: -kve(1, s) * E1(s)[0] / 3)
    W1 = vec(lambda s: 0.5 * (kve(0, s) - kve(1, s) / s) * E1(s)[0])
    W2 = vec(lambda s: -(2.0 / 3.0) * kve(1, s) * E1(s)[0])
    M0, dM0 = tail(W0, 0, x, 2.0)
    M1, dM1 = tail(W1, 1, x, 2.0)
    M2, dM2 = tail(W2, 2, x, 2.0)
    w0, w2 = W0(x)[0], W2(x)[0]
    ddM0 = -w0 - 2 * dM0 / x
    ddM2 = -w2 - 2 * dM2 / x + 6 * M2 / x ** 2
    g = dM2 + 3 * M2 / x
    dg = ddM2 + 3 * dM2 / x - 3 * M2 / x ** 2
    return x ** 2 * (M0 - (dM1 + 2 * M1 / x) / 3 + (ddM0 + 2 * dM0 / x) / 3 + (2.0 / 15.0) * (dg + 2 * g / x))


def richardson(G, b=1.0):
    xs = [50.0, 100.0, 200.0, 400.0]
    vals = [G(x) / x ** b for x in xs]
    r1 = [2 * vals[i + 1] - vals[i] for i in range(3)]
    r2 = [(4 * r1[i + 1] - r1[i]) / 3 for i in range(2)]
    C = (8 * r2[1] - r2[0]) / 7
    slope = math.log(G(800.0) / G(400.0)) / math.log(2.0)
    return C, slope


if __name__ == "__main__":
    # (name, G, Pi, prod k^2, principal average x (-(n-1))^{max degree of root} -- the 'homogeneous' formula)
    cases = [("N2", G_N2, 1 / 3, 1, (-2) * 1 / 3), ("N3", G_N3, 7 / 3, 1, 4 * 1 / 3), ("N6", G_N6, 4 / 5, 4, 0.0)]
    for nm, G, Pi, k2, wrong in cases:
        C, slope = richardson(G)
        pred = (math.pi / 2) * Pi / k2
        alt = (math.pi / 2) * wrong / k2
        print(f"{nm}: C numeric {C:.12f}  predicted C*Pi {pred:.12f}  rel {abs(C / pred - 1):.1e}  "
              f"principal-symbol formula {alt:.6f}  local exponent {slope:.4f}")
