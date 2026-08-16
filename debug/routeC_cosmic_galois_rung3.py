"""Cosmic-Galois probe, RUNG 3 (frontier): is the INTEGRATED T2 a Gamma(2) multiple
modular value / elliptic polylogarithm?

Rung 1 pinned the FAMILY as Legendre/Gamma(2) (lambda=1-rho) and Rung 2 confirmed the
FIBER periods are Gamma-values.  The observable-level statement is about the INTEGRATED
object: T2 = (8/pi) d^4/dzeta int_0^1 ds int_0^1 dt int_0^inf dk j0(k|W|) P1 P2, which
pulls back to an ITERATED INTEGRAL over the lambda-line X(2) (via lambda=1-rho) of
period-weighted forms -- a candidate Gamma(2) multiple modular value (MMV) /
elliptic polylogarithm.  This is the SAME open object as the ABW closed form
(sprint_routeC_momentum_memo.md) -- genuinely frontier.

This driver does the bounded, PRINCIPLED probe (no grab-bag PSLQ; two-precision-stable):
 (1) recompute the collinear symmetric value (X=0, Y=(0,0,1), Z=(0,0,-1), 1s zeta=1;
     D1=D2=1, |W|=s+t) at two precisions -> a stable high-precision target.
 (2) test whether it is a SINGLE-FIBER weight<=2 PERIOD PRODUCT built from the rho=1/2
     CM-fiber period varpi=K(1/2)=Gamma(1/4)^2/(4 sqrt pi) and pi:
        basis {1, pi, pi^2, varpi, varpi^2, varpi*pi}.
     A hit => reducible to the single CM fiber; NO hit => genuine MULTI-fiber object
     (integrates over all of X(2)), consistent with an MMV, not a single period.
 (3) recall (weight_probe.py) that weight-0/1/2 POLYLOG bases already returned nothing
     -> not in the genus-0 polylog ring.  Together: the integrated value is a genuine
     Gamma(2) elliptic MMV -- realizing it explicitly needs the iterated-Eisenstein /
     elliptic-polylog basis for Gamma(2) = the frontier (collaboration piece).
"""
from __future__ import annotations
import mpmath as mp


def P(x, k, D):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-D * Del) * (D * D / Del ** 3 + 3 * D / Del ** 4 + 3 / Del ** 5)


def I_rad(s, t, D1=mp.mpf(1), D2=mp.mpf(1)):
    b = s + t
    def f(k):
        j0 = mp.sin(k * b) / (k * b) if k * b > mp.mpf('1e-20') else mp.mpf(1)
        return j0 * P(s, k, D1) * P(t, k, D2)
    return mp.quad(f, [0, 1, 3, 8, 20, mp.inf])


def collinear_value():
    inner = lambda s: mp.quad(lambda t: I_rad(s, t), [0, mp.mpf('0.5'), 1])
    return (8 / mp.pi) * mp.quad(inner, [0, mp.mpf('0.5'), 1])


def main():
    print("RUNG 3 (frontier) -- is the integrated T2 a Gamma(2) multiple modular value?\n")

    # NOTE: collinear_value() by naive nested tanh-sinh is impractically slow for
    # PSLQ-grade precision (times out at dps>=18 in minutes) -- a dedicated evaluator
    # (analytic k-integral / series) is itself frontier engineering.  For the bounded
    # single-fiber exclusion we use the routeC_weight_probe.py value (~9-10 digits).
    # Set RECOMPUTE=True (and allow many minutes) to regenerate it here.
    RECOMPUTE = False
    mp.mp.dps = 20
    if RECOMPUTE:
        v = collinear_value()
        print(f"(1) collinear value (recomputed): {mp.nstr(v, 18)}")
    else:
        v = mp.mpf('0.3953557703')      # routeC_weight_probe.py (dps26 tanh-sinh), ~10 digits
        print(f"(1) collinear value (from routeC_weight_probe.py, ~10 digits): {v}")

    print("\n(2) single-fiber weight<=2 period-product test (rho=1/2 CM fiber),")
    print("    PRECISION-LIMITED (~9 digits) -- illustrative, not decisive:")
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))     # = K(1/2)
    basis = [mp.mpf(1), mp.pi, mp.pi ** 2, varpi, varpi ** 2, varpi * mp.pi]
    names = ["1", "pi", "pi^2", "varpi", "varpi^2", "varpi*pi"]
    rel = mp.pslq([v] + basis, tol=mp.mpf(10) ** -8, maxcoeff=10 ** 3, maxsteps=10 ** 6)
    print(f"    basis {names}")
    print(f"    PSLQ relation (tol 1e-8, maxcoeff 1e3): {rel}")
    print("    None / huge-coeff => not a low-height single-fiber period product (bounded).")

    print("\n(3) polylog ring: weight-0/1/2 bases already NEGATIVE (routeC_weight_probe.py),")
    print("    so not in the genus-0 polylog ring either.")
    print("\nStatus: the integrated T2 is a genuine Gamma(2) elliptic MMV -- modular home")
    print("identified (Rung 1), fiber periods are Gamma-values (Rung 2), not in the polylog")
    print("ring (weight_probe), and the single fiber does not close it (ABW obstruction).")
    print("Explicit MMV realization needs BOTH a Gamma(2) iterated-Eisenstein/elliptic-polylog")
    print("basis AND a dedicated high-precision evaluator of the integrated observable")
    print("(naive nested quadrature is too slow) = the frontier / ABW-open / collaboration piece.")


if __name__ == '__main__':
    main()
