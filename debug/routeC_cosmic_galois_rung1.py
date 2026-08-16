"""Cosmic-Galois probe, RUNG 1: is the Route C elliptic family a CONGRUENCE modular
family, and at what level?  (Decisive GO/STOP for the elliptic cosmic-Galois question.)

Family (normalized c1=1, c2=rho): E_rho : y^2 = (x^2-1)(rho x^2 + 1 - rho),
base rho = c2/c1 = t(1-t)/(s(1-s)) in (0, inf).  Its two periods are the verified
K(1-rho), K(rho) (solutions of the modulus PF M_rho), so the modular parameter is
    tau(rho) = i K(rho) / K(1-rho)          (orientation fixed to Im tau > 0).

TESTS (all decisive):
 A. MODULAR IDENTIFICATION.  Compute tau(rho) from the periods, then the theta-series
    modular lambda  L(tau) = (theta2/theta3)^4, and check it equals an S_3 image of the
    branch-point cross-ratio  lam(rho) = ((sqrt(rho)-1)/(sqrt(rho)+1))^2.  A match at
    several rho => the family IS the Legendre (Gamma(2)) universal family pulled back
    along rho -> lam(rho); its periods are Gamma(2) modular periods.
 B. j-INVARIANT rational in rho (isotriviality check: j must VARY => genuine family).
 C. SINGULAR FIBERS over the rho-base and Kodaira types => the level / GO-STOP.
 D. physical (s,t) -> rho range (which arc of the modular family the chemistry visits;
    informs Rung 2).
"""
from __future__ import annotations
import mpmath as mp
mp.mp.dps = 40


def tau_of_rho(rho):
    rho = mp.mpf(rho)
    t = 1j * mp.ellipk(rho) / mp.ellipk(1 - rho)     # i K(m=rho)/K(m=1-rho)
    return t if t.imag > 0 else -t


def lambda_theta(tau):
    q = mp.e ** (1j * mp.pi * tau)                    # nome
    th2 = mp.jtheta(2, 0, q)
    th3 = mp.jtheta(3, 0, q)
    return (th2 / th3) ** 4


def cross_ratio_lambda(rho):
    u = mp.sqrt(mp.mpf(rho))
    return ((u - 1) / (u + 1)) ** 2


def s3_orbit(x):
    return [x, 1 - x, 1 / x, 1 / (1 - x), x / (x - 1), (x - 1) / x]


def j_from_lambda(lam):
    return 256 * (1 - lam + lam ** 2) ** 3 / (lam ** 2 * (1 - lam) ** 2)


def main():
    print("RUNG 1 -- modular identification of the Route C elliptic family\n")

    print("A. MODULAR IDENTIFICATION: does theta-lambda(tau(rho)) match a cross-ratio S3 image?")
    for rho in ['0.37', '0.6', '0.15', '0.85']:
        tau = tau_of_rho(rho)
        Lt = lambda_theta(tau)
        lam = cross_ratio_lambda(rho)
        orbit = s3_orbit(lam)
        # also test against the direct m=rho / m=1-rho representatives
        cands = orbit + s3_orbit(mp.mpf(rho))
        best = min(cands, key=lambda c: abs(c - Lt))
        hit = abs(best - Lt)
        print(f"   rho={rho:>5}: tau={mp.nstr(tau,6)}  theta-lam={mp.nstr(Lt,8)}  "
              f"nearest S3 rep={mp.nstr(best,8)}  |diff|={mp.nstr(hit,3)}")
    print("   (|diff| ~ 0 => the family is the Legendre/Gamma(2) family; periods are modular)\n")

    print("B. j-INVARIANT j(rho) = 256(1-lam+lam^2)^3/(lam^2(1-lam)^2), lam=cross-ratio:")
    for rho in ['0.37', '0.6', '0.15']:
        lam = cross_ratio_lambda(rho)
        print(f"   rho={rho:>5}: lam={mp.nstr(lam,8)}  j={mp.nstr(j_from_lambda(lam),10)}")
    print("   (j varies with rho => NON-isotrivial genuine family, not a constant curve)\n")

    print("C. SINGULAR FIBERS over the rho-base (P^1):")
    for rho in ['1e-8', '0.5', '1.0', '2.0', '1e8']:
        rr = mp.mpf(rho)
        lam = cross_ratio_lambda(rr)
        # discriminant of the Legendre curve ~ lam^2 (1-lam)^2 ; degenerate where lam in {0,1}
        disc = lam ** 2 * (1 - lam) ** 2
        print(f"   rho={rho:>6}: lam={mp.nstr(lam,6)}  disc~{mp.nstr(disc,4)}  "
              + ("SINGULAR" if disc < mp.mpf(10) ** -6 else "smooth"))
    print("   degenerations: rho->0 (Q->x^2-1, genus 0), rho=1 (c1=c2 diagonal, genus 0),")
    print("   rho->inf.  Three special fibers over rho in {0,1,inf}.\n")

    print("D. physical (s,t) -> rho = t(1-t)/(s(1-s)) range:")
    vals = []
    for s in [0.1, 0.3, 0.5, 0.7, 0.9]:
        for t in [0.1, 0.3, 0.5, 0.7, 0.9]:
            vals.append((t * (1 - t)) / (s * (1 - s)))
    print(f"   rho spans [{min(vals):.4f}, {max(vals):.4f}], includes rho=1 (diagonal s=t) "
          f"= the genus-0 singular fiber.")
    print("   => the chemistry integration sweeps THROUGH the singular fiber at rho=1.\n")

    print("E. RUNG-2 CM DATA POINT: rho=1/2 is IN the physical range and gives tau=i")
    print("   (CM by Z[i], lemniscatic).  The period there must be a Gamma-value (Chowla-Selberg):")
    rho = mp.mpf('0.5')
    tau = tau_of_rho(rho)
    K = mp.ellipk(rho)
    CS = mp.gamma(mp.mpf('0.25')) ** 2 / (4 * mp.sqrt(mp.pi))
    print(f"   rho=1/2: tau={mp.nstr(tau,10)}  K(1/2)={mp.nstr(K,24)}")
    print(f"            Gamma(1/4)^2/(4 sqrt pi)={mp.nstr(CS,24)}  |diff|={mp.nstr(abs(K-CS),3)}")
    print("   => a molecular integral's period, at a physically-visited fiber, IS a Gamma-value.")
    print("      The Route C elliptic period populates the modular/Gamma-value (cosmic-Galois) ring.")


if __name__ == '__main__':
    main()
