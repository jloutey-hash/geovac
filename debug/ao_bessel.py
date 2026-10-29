"""Modified spherical Bessel functions at half-integer order, elementary (record AO2).
i_L(x) = x^L sum_m (x^2/2)^m / (m! (2L+2m+1)!!)  (series, all terms positive);
k_L(x): k_0 = e^{-x}/x, k_1 = e^{-x}(x+1)/x^2, k_{L+1} = k_{L-1} + (2L+1)/x k_L (upward, stable).  Conventions as
debug/ao_radial.py: i_0 = sinh x / x, k_0 = e^{-x}/x."""
import flint
from flint import arb


def il_pair(L, x):
    """(i_L(x), i_{L+1}(x)) by the series."""
    out = []
    for M in (L, L + 1):
        df = arb(1)
        for k in range(1, 2 * M + 2, 2):
            df *= k
        term = x ** M / df
        s = term
        h = x * x / 2
        m = 0
        tol = 2.0 ** (-flint.ctx.prec - 10)
        while True:
            m += 1
            term = term * h / (m * (2 * M + 2 * m + 1))
            s += term
            if float((term / s).mid()) <= tol:      # ratio, not values: i_L underflows a float at large L
                break
        out.append(s)
    return out


def kl_pair(L, x):
    """(k_L(x), k_{L+1}(x)) by upward recurrence."""
    e = (-x).exp()
    k0, k1 = e / x, e * (x + 1) / (x * x)
    for l in range(1, L + 1):
        k0, k1 = k1, k0 + (2 * l + 1) / x * k1
    return k0, k1
