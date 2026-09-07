"""Paper 61 (companion of Paper 59) sec:modular -- the conductor-4 arithmetic of the Legendre fibre period.

Backs the corrected justification for admitting G = beta(2) = L(2, chi_-4) into the
PSLQ ring (v4.102 audit, debug/sprint_qi_seam_audit_memo.md):

  (i)   theta_3(0,q)^2 is the theta series of Z[i]: r_2(n) = 4 sum_{d|n} chi_-4(d)
        (Jacobi two-squares), so theta_3^2 is the weight-1 level-4 Eisenstein series
        with character chi_-4 and L(theta_3^2, s) = 4 zeta(s) beta(s).
  (ii)  The fibre period IS that form: K(lambda) = (pi/2) theta_3(0,q)^2.
  (iii) The physical Bessel twist promotes the level 2 -> 4:
        sqrt(rho) = sqrt(1-lambda) = theta_4^2/theta_3^2 (a Gamma(4)-level function),
        via the decay rate A = sqrt(c1)+sqrt(c2) = sqrt(u)(1+sqrt(rho)).
  (iv)  Convention pin: the paper's K(1/2) = Gamma(1/4)^2/(4 sqrt pi) is the PARAMETER
        convention, in which int_0^1 K dm = 2 and int_0^1 E dm = 4/3 (rational);
        the Catalan identities int_0^1 K(k) dk = 2G, int_0^1 E(k) dk = G + 1/2 hold in
        the MODULUS convention only.  (The physical co-area measure d(rho) is the
        parameter one.)
"""
from mpmath import (mp, mpf, ellipk, ellipe, gamma, sqrt, pi, quad, catalan,
                    jtheta, exp)


def _chi4(d: int) -> int:
    return 0 if d % 2 == 0 else (1 if d % 4 == 1 else -1)


def test_theta3_sq_is_Zi_theta_series():
    """theta_3(0,q)^2 = 1 + sum_n r_2(n) q^n with r_2(n) = 4 sum_{d|n} chi_-4(d)."""
    mp.dps = 30
    q = mpf('0.13')
    r2 = lambda n: 4 * sum(_chi4(d) for d in range(1, n + 1) if n % d == 0)
    series = 1 + sum(r2(n) * q**n for n in range(1, 220))
    assert abs(jtheta(3, 0, q)**2 - series) < mpf(10)**(-25)


def test_fibre_period_is_the_chi4_form():
    """K(lambda(tau)) = (pi/2) theta_3(0,q)^2 at a generic (non-CM) tau."""
    mp.dps = 40
    tau = mpf('1.3') * 1j
    q = exp(1j * pi * tau)
    th2, th3 = jtheta(2, 0, q), jtheta(3, 0, q)
    lam = (th2 / th3)**4
    assert abs(ellipk(lam) - pi / 2 * th3**2) < mpf(10)**(-35)


def test_bessel_twist_promotes_level_2_to_4():
    """sqrt(1-lambda) = theta_4^2/theta_3^2 (Gamma(4)); and the twist's decay rate
    A = sqrt(c1)+sqrt(c2) = sqrt(u)(1+sqrt(rho)) makes sqrt(rho) physical."""
    mp.dps = 40
    tau = mpf('1.3') * 1j
    q = exp(1j * pi * tau)
    th2, th3, th4 = jtheta(2, 0, q), jtheta(3, 0, q), jtheta(4, 0, q)
    lam = (th2 / th3)**4
    assert abs(sqrt(1 - lam) - th4**2 / th3**2) < mpf(10)**(-35)
    # decay rate: c1 = u, c2 = u*rho  =>  sqrt(c1)+sqrt(c2) = sqrt(u)(1+sqrt(rho))
    u, rho = mpf('0.37'), mpf('0.22')
    assert abs((sqrt(u) + sqrt(u * rho)) - sqrt(u) * (1 + sqrt(rho))) < mpf(10)**(-38)


def test_convention_pin_parameter_vs_modulus():
    """The paper's K(1/2) is parameter-convention; the G-identities are modulus-only."""
    mp.dps = 40
    lemniscatic = gamma(mpf(1) / 4)**2 / (4 * sqrt(pi))
    assert abs(ellipk(mpf(1) / 2) - lemniscatic) < mpf(10)**(-38)      # parameter conv
    assert abs(quad(lambda m: ellipk(m), [0, 1]) - 2) < mpf(10)**(-30)             # rational
    assert abs(quad(lambda m: ellipe(m), [0, 1]) - mpf(4) / 3) < mpf(10)**(-30)    # rational
    assert abs(quad(lambda k: ellipk(k**2), [0, 1]) - 2 * catalan) < mpf(10)**(-30)         # 2G
    assert abs(quad(lambda k: ellipe(k**2), [0, 1]) - (catalan + mpf(1) / 2)) < mpf(10)**(-30)
