"""Fast K0, K1 and L(X) = int_0^X K0 at high precision for the P3e-2W grids.

Same functions as debug/three_electron_p3e_check.py (K0, K1 via mpmath.besselk; L via the Struve
form), evaluated by a faster route because mpmath's integer-order besselk costs 1.5-3.5 s per call
at 430 digits for 20 < X < 600 (measured 2026-10-19), which would put a triple at ~25 h.
  small/moderate X: fixed-point (gmpy2 integer) power series, shared terms t_k = (X/2)^{2k}/k!^2,
      K0 = -(ln(X/2)+gamma) I0 + sum t_k H_k,   K1 = (1/X - I1 K0)/I0 (Wronskian),
      L  = X sum t_k/(2k+1) [-(ln(X/2)+gamma) + 1/(2k+1) + H_k]  (termwise integral of K0),
      guard bits for the e^{2X} cancellation;
  large X (2X/ln10 > p + 10): Hankel asymptotic series for K0, K1 and for Ki(X) = int_X^inf K0,
      L = pi/2 - Ki, coefficients b_k = a_k(0) - (k - 1/2) b_{k-1}.
p = relative digits required.  Validated against mpmath in selftest().
"""
from __future__ import annotations
import sys
import gmpy2
import mpmath as mp

mpz = gmpy2.mpz


def _series(X, p: int):
    """(K0, K1, L) at X > 0 to ~p relative digits by the fixed-point power series."""
    Xf = float(X)
    extra = int(2 * Xf / 2.302585) + 30
    bits = int((p + extra) * 3.3219281) + 64
    with mp.workprec(bits + 32):
        Xm = mp.mpf(X)
        one = mpz(1) << bits
        Y = mpz(mp.floor(Xm * Xm / 4 * mp.mpf(2) ** bits))
        lg = mp.log(Xm / 2) + mp.euler
        t = one
        H = mpz(0)
        I0 = mpz(0); SH = mpz(0); I1 = mpz(0); A1 = mpz(0); A2 = mpz(0); A3 = mpz(0)
        k = 0
        while True:
            I0 += t
            I1 += t // (k + 1)
            tH = (t * H) >> bits
            SH += tH
            o = 2 * k + 1
            A1 += t // o
            A2 += t // (o * o)
            A3 += tH // o
            k += 1
            t = ((t * Y) >> bits) // (k * k)
            if t == 0 and k > Xf:
                break
            H += one // k
        sc = mp.mpf(2) ** (-bits)
        I0m, SHm = mp.mpf(I0) * sc, mp.mpf(SH) * sc
        I1m = mp.mpf(I1) * sc * Xm / 2
        K0 = -lg * I0m + SHm
        K1 = (1 / Xm - I1m * K0) / I0m
        L = Xm * (-lg * mp.mpf(A1) * sc + mp.mpf(A2) * sc + mp.mpf(A3) * sc)
    return +K0, +K1, +L


def _asym(X, p: int):
    """(K0, K1, L) by the large-X asymptotic series; requires 2X/ln10 > p + 10.
    K0, K1 to p relative digits; Ki = int_X^inf K0 to absolute 10^-(p+8) (L = pi/2 - Ki ~ 1)."""
    with mp.workdps(p + 15):
        Xm = mp.mpf(X)
        tol = mp.mpf(10) ** (-(p + 8))
        pre = mp.sqrt(mp.pi / (2 * Xm)) * mp.exp(-Xm)

        def hankel(nu2):
            s, a, xp, k, last = mp.mpf(0), mp.mpf(1), mp.mpf(1), 0, mp.inf
            while True:
                term = abs(a) * xp
                if term < tol:
                    return s
                if term > last and k > 4:
                    raise ValueError("Hankel series not converged")
                last = term
                s += a * xp
                k += 1
                a = a * (nu2 - (2 * k - 1) ** 2) / (8 * k)
                xp = xp / Xm
        s0, s1 = hankel(0), hankel(4)
        # Ki: b_0 = 1, b_k = a_k(0) - (k - 1/2) b_{k-1}; absolute tolerance on pre * sum
        sk, a0, b, xp, k, last = mp.mpf(0), mp.mpf(1), mp.mpf(1), mp.mpf(1), 0, mp.inf
        atol = tol / pre
        while True:
            term = abs(b) * xp
            if term < atol:
                break
            if term > last and k > 4:
                raise ValueError("Ki series not converged")
            last = term
            sk += b * xp
            k += 1
            a0 = a0 * (0 - (2 * k - 1) ** 2) / (8 * k)
            b = a0 - (k - mp.mpf(1) / 2) * b
            xp = xp / Xm
        K0, K1, Ki = pre * s0, pre * s1, pre * sk
        L = mp.pi / 2 - Ki
    return +K0, +K1, +L


def KKL(X, p: int):
    """K0(X), K1(X), L(X) to p relative digits (results at the current mp precision)."""
    if 2 * float(X) / 2.302585 > p + 10:
        return _asym(X, p)
    return _series(X, p)


def selftest(dps: int = 430) -> float:
    sys.path.insert(0, "debug")
    from three_electron_p3e_check import L as Lref
    mp.mp.dps = dps
    worst = mp.mpf(0)
    for xs in ["1e-30", "1e-8", "0.001", "0.37", "1", "2.5", "7.3", "19", "44.4", "130", "333", "480",
               "520", "560", "610", "900"]:
        X = mp.mpf(xs)
        K0, K1, L = KKL(X, dps)
        r0, r1 = mp.besselk(0, X), mp.besselk(1, X)
        rL = Lref(X)
        e = max(abs(K0 / r0 - 1), abs(K1 / r1 - 1), abs(L / rL - 1))
        worst = max(worst, e)
        print(f"X={xs:6s} {'asym' if 2*float(X)/2.302585 > dps+10 else 'series'}  rel err {mp.nstr(e, 3)}", flush=True)
    print("worst", mp.nstr(worst, 3))
    return worst


if __name__ == "__main__":
    selftest(int(sys.argv[1]) if len(sys.argv) > 1 else 430)
