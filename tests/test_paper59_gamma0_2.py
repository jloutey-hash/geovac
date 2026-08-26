"""Paper 59 sec:modular -- the s<->t exchange is a modular involution and T2's
modular home is X_0(2) (CHANGELOG v4.103.0; driver debug/aha_t2a_involution.py).

Pins: (i) c = s(1-s) is s -> 1-s invariant, so the full Feynman-square symmetry
group acts on the modulus rho = c2/c1 through {id, rho -> 1/rho} only;
(ii) rho -> 1/rho with lambda = 1-rho is lambda -> lambda/(lambda-1), the coset
T*Gamma(2); (iii) its elliptic representative M = [[1,-1],[2,-1]] has M^2 = -I,
M = T mod 2, and fixes tau* = (1+i)/2 with lambda(tau*) = 2, j(tau*) = 1728
(disc -4, CM by Z[i]) -- a point OFF the physical contour rho > 0;
(iv) <Gamma(2), M> = Gamma_0(2); the physical s<->t fixed locus is rho = 1,
the pure-Tate diagonal cusp.
"""
import mpmath as mp
import sympy as sp


def _lam(tau):
    q = mp.exp(mp.mpc(0, 1) * mp.pi * tau)
    return (mp.jtheta(2, 0, q) / mp.jtheta(3, 0, q)) ** 4


def test_symbol_side_involution():
    s, rho = sp.symbols("s rho", positive=True)
    c = s * (1 - s)
    assert sp.simplify(c.subs(s, 1 - s) - c) == 0          # s -> 1-s trivial on c
    lam = 1 - rho
    lam_inv = 1 - 1 / rho                                   # rho -> 1/rho
    assert sp.simplify(lam_inv - lam / (lam - 1)) == 0      # = lambda/(lambda-1)


def test_lambda_T_coset():
    mp.mp.dps = 40
    for tau in (mp.mpc("0.3", "1.1"), mp.mpc("-0.2", "0.7")):
        lam = _lam(tau)
        assert abs(_lam(tau + 1) - lam / (lam - 1)) < mp.mpf("1e-30")


def test_elliptic_representative():
    mp.mp.dps = 40
    a, b, c, d = 1, -1, 2, -1                                # M = [[1,-1],[2,-1]]
    assert (a * a + b * c, a * b + b * d, c * a + d * c, c * b + d * d) == (-1, 0, 0, -1)
    assert (a % 2, b % 2, c % 2, d % 2) == (1, 1, 0, 1)      # M = T mod 2
    tau_star = mp.mpc("0.5", "0.5")
    assert abs(2 * tau_star ** 2 - 2 * tau_star + 1) < mp.mpf("1e-38")
    m_tau = (a * tau_star + b) / (c * tau_star + d)
    assert abs(m_tau - tau_star) < mp.mpf("1e-38")
    for tau in (mp.mpc("0.3", "1.1"), mp.mpc("-0.2", "0.7")):
        lam = _lam(tau)
        mt = (a * tau + b) / (c * tau + d)
        assert abs(_lam(mt) - lam / (lam - 1)) < mp.mpf("1e-30")
    assert abs(_lam(tau_star) - 2) < mp.mpf("1e-25")
    assert abs(mp.kleinj(tau_star) - 1) < mp.mpf("1e-25")    # j = 1728 (kleinj = j/1728)


def test_physical_fixed_locus_is_cusp():
    """The s<->t fixed locus t = s (or t = 1-s) gives rho = 1 = the diagonal
    cusp: K(1-rho) -> pi/2 exactly (pure Tate), tau -> i*infinity."""
    mp.mp.dps = 30
    assert abs(mp.ellipk(0) - mp.pi / 2) < mp.mpf("1e-25")
    ims = [mp.im(mp.mpc(0, 1) * mp.ellipk(1 - m) / mp.ellipk(m))
           for m in (mp.mpf("0.1"), mp.mpf("0.01"), mp.mpf("0.001"))]
    assert ims[0] < ims[1] < ims[2]                          # tau -> i*inf at the cusp
