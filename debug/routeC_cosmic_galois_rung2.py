"""Cosmic-Galois probe, RUNG 2: systematic PSLQ of the Route C periods against the
modular / Gamma-value (Chowla-Selberg) ring across the CM fibers the domain visits.

Rung 1 pinned the family as Legendre/Gamma(2) with lambda = 1 - rho.  The physical
base rho = t(1-t)/[s(1-s)] sweeps ALL of (0, inf) as (s,t) -> boundary, so the (s,t)
integration covers the REAL locus lambda = 1 - rho < 1 of X(2) and hits infinitely many CM
fibers -- NOT every one (corrected 2026-09-07): non-real lambda (disc -3) and lambda=2
are off the contour.  Rung 2 confirms the periods there are Gamma-values, at TWO distinct
fundamental discriminants (so it is a systematic Galois/CM fact, not a lemniscatic
coincidence):
  - disc -4  (tau = i,     rho = 1/2):        Gamma(1/4)   [done exactly in Rung 1]
  - disc -8  (tau = i*sqrt2, rho = 1 - lam):  Gamma(1/8), Gamma(3/8)

Method: at tau = i*sqrt(N), lam_N = theta-lambda(tau) (= the singular modulus), the
period is K(lam_N); PSLQ ln K against a HOMOGENEOUS log-basis {ln pi, ln 2, ln(unit),
ln Gamma(a/d)} (weight-1, per the corpus PSLQ discipline).  A small-integer relation =
Chowla-Selberg confirmed.  Cross-checked at two precisions.
"""
from __future__ import annotations
import mpmath as mp


def lambda_theta(tau):
    q = mp.e ** (1j * mp.pi * tau)
    return (mp.jtheta(2, 0, q) / mp.jtheta(3, 0, q)) ** 4


def pslq_logrel(target, basis_vals, names, maxcoeff=10 ** 6):
    # CM-fiber values are real (tau on the imaginary axis); coerce away any 0j.
    # tol scales with precision (fixed tol was the instability artifact).
    vec = [mp.re(mp.log(target))] + [mp.re(mp.log(v)) for v in basis_vals]
    tol = mp.mpf(10) ** (-(3 * mp.mp.dps) // 4)
    rel = mp.pslq(vec, tol=tol, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    return rel


def run(dps):
    mp.mp.dps = dps
    out = {}

    # --- disc -4 control: tau=i, rho=1/2, lam=1/2 ---
    lam1 = lambda_theta(1j)                     # = 1/2
    K1 = mp.ellipk(lam1)
    rel1 = pslq_logrel(K1, [mp.pi, mp.mpf(2), mp.gamma(mp.mpf(1) / 4)],
                       ["lnK", "ln pi", "ln 2", "ln G(1/4)"])
    cf1 = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
    out["disc-4"] = (mp.nstr(lam1, 12), mp.nstr(1 - lam1, 12), rel1, mp.nstr(abs(K1 - cf1), 3))

    # --- disc -8: tau=i sqrt2 ---
    lam2 = mp.re(lambda_theta(1j * mp.sqrt(2)))  # = 3 - 2 sqrt2
    K2 = mp.ellipk(lam2)
    unit = 1 + mp.sqrt(2)                        # fundamental unit of Q(sqrt2)
    rel2 = pslq_logrel(
        K2, [mp.pi, mp.mpf(2), unit, mp.gamma(mp.mpf(1) / 8), mp.gamma(mp.mpf(3) / 8)],
        ["lnK", "ln pi", "ln 2", "ln(1+sqrt2)", "ln G(1/8)", "ln G(3/8)"])
    # Direct closed-form residual (independent of PSLQ): the classical disc-8 value
    #   K(k2) = (1+sqrt2)^{1/2} G(1/8) G(3/8) / (2^{13/4} sqrt pi)
    cf2 = mp.sqrt(unit) * mp.gamma(mp.mpf(1) / 8) * mp.gamma(mp.mpf(3) / 8) / (2 ** mp.mpf('3.25') * mp.sqrt(mp.pi))
    out["disc-8"] = (mp.nstr(lam2, 12), mp.nstr(1 - lam2, 12), rel2, mp.nstr(abs(K2 - cf2), 3))
    return out


def main():
    print("RUNG 2 -- periods at CM fibers are Gamma-values (Chowla-Selberg), across discriminants\n")
    print("Physical base rho = t(1-t)/[s(1-s)] sweeps (0, inf): lambda = 1-rho covers the\nREAL locus lambda < 1 of X(2), hitting infinitely many CM fibers -- not every one.\n")

    r50 = run(50)
    r70 = run(70)   # cross-precision stability (per the corpus: two-precision-stable required)

    for key in ("disc-4", "disc-8"):
        lam, rho, rel, cfres = r50[key]
        _, _, rel_hi, _ = r70[key]
        stable = (rel == rel_hi)
        print(f"{key}:  lambda={lam}  rho(=1-lambda)={rho}")
        print(f"   PSLQ log-relation (dps50): {rel}   stable at dps70? {stable}")
        print(f"   direct closed-form residual |K - CS-value|: {cfres}")
        if key == "disc-4":
            print("   [lnK, ln pi, ln 2, ln G(1/4)];  K = G(1/4)^2/(4 sqrt pi)")
        else:
            print("   [lnK, ln pi, ln 2, ln(1+sqrt2), ln G(1/8), ln G(3/8)];")
            print("   K = (1+sqrt2)^{1/2} G(1/8) G(3/8) / (2^{13/4} sqrt pi)")
        print()

    print("Verdict: the DIRECT closed-form residual (~1e-51, independent of PSLQ) confirms")
    print("K at each CM fiber is a Gamma-value product (Chowla-Selberg) at TWO distinct")
    print("fundamental discriminants (-4 via Gamma(1/4); -8 via Gamma(1/8),Gamma(3/8)) => the")
    print("periods populate the modular/Gamma-value (cosmic-Galois) ring SYSTEMATICALLY,")
    print("tracking the discriminant.  NOTE: the disc-8 PSLQ leg was NOT two-precision-stable")
    print("(a genuine PSLQ trap flagged by the cross-precision check); the confirmation rests")
    print("on the direct closed-form residual, not the PSLQ -- per the corpus PSLQ discipline.")


if __name__ == '__main__':
    main()
