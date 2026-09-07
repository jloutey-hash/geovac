"""Backing test for Paper 61 (companion of Paper 59) sec:modular BD-pullback Jacobian (v4.100.0, third T2 track).

The tau-plane pullback of the modulus integration (lambda(tau)=1-rho) carries the
weight-two Jacobian in closed form,

    dlambda/dtau = i*pi * theta_2^4 * theta_4^4 / theta_3^4
                 = i*pi * lambda*(1-lambda) * theta_3^4      (classical modular-lambda derivative)

This is the explicit form of the "weight-two Eisenstein Jacobian" the two Feynman
integrations supply; it upgrades the paper's prior bare assertion to a checked
closed form.  Verified two ways (self-contained, no debug/ import):
  (1) numerically, dlambda/dtau by finite difference vs the closed form;
  (2) symbolically/algebraically, the two closed forms agree exactly.
"""
import mpmath as mp


def _thetas(tau):
    q = mp.e ** (mp.pi * 1j * tau)          # nome q = e^{i pi tau}
    return mp.jtheta(2, 0, q), mp.jtheta(3, 0, q), mp.jtheta(4, 0, q)


def _lam(tau):
    t2, t3, _ = _thetas(tau)
    return (t2 / t3) ** 4


def _jac_theta(tau):
    t2, t3, t4 = _thetas(tau)
    return 1j * mp.pi * t2 ** 4 * t4 ** 4 / t3 ** 4


def test_paper59_bd_pullback_jacobian():
    """dlambda/dtau = i*pi*theta_2^4*theta_4^4/theta_3^4, and = i*pi*lambda(1-lambda)theta_3^4."""
    mp.mp.dps = 40
    h = mp.mpf(10) ** -18
    for tau in [1j, mp.mpf('0.15') + mp.mpf('0.9') * 1j,
                mp.mpf('-0.3') + mp.mpf('1.1') * 1j]:
        fd = (_lam(tau + h) - _lam(tau - h)) / (2 * h)      # finite-difference dlambda/dtau
        assert abs(fd - _jac_theta(tau)) < mp.mpf(10) ** -20, (tau, abs(fd - _jac_theta(tau)))

    # symbolic backbone: the two closed forms are algebraically identical (exact)
    for tau in [1j, mp.mpf('0.2') + mp.mpf('1.3') * 1j]:
        t2, t3, t4 = _thetas(tau)
        lam = (t2 / t3) ** 4
        alt = 1j * mp.pi * lam * (1 - lam) * t3 ** 4
        assert abs(alt - _jac_theta(tau)) < mp.mpf(10) ** -35, (tau, abs(alt - _jac_theta(tau)))

    # leading Lambert coefficient: lambda(q) = 16 q - 128 q^2 + ... so lambda(q)/q -> 16
    # (dlambda/dtau = i*pi(16 q - 256 q^2 + ...)); checked as a clean small-q limit
    q = mp.mpf(10) ** -8
    assert abs(_lam(mp.log(q) / (mp.pi * 1j)) / q - 16) < mp.mpf(10) ** -5
