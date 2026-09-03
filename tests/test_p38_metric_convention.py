"""Paper 38 sec:ch_triple + Lemma L2 -- the metric is the DUAL-COXETER one,
and its geodesic distance is the rotation angle.

PI direction 2026-09-03 (Option 1): Paper 38 is restated on the dual-Coxeter
sphere throughout, so that its moment and its Lipschitz seminorm live in one
metric and the paper is exactly the rank-1 case of Paper 40.

The derivation pinned here.  On su(2) the dual-Coxeter rule Cas(ad) = h^v = 2
fixes the Ad-invariant inner product to <X,Y> = -2 tr(XY) in the defining
representation.  In that inner product the geodesic distance from the identity
to exp(i theta n.sigma) -- a group element with eigenphases +-theta, i.e. at
unit-S^3 distance theta -- is 2 theta, the SO(3) rotation angle chi.  So:

  * the dual-Coxeter metric on SU(2) is the round 3-sphere of RADIUS 2
    (diameter 2 pi, volume 16 pi^2), not the unit sphere;
  * `central_fejer_su2.gamma_rate`, which integrates against chi, computes the
    dual-Coxeter moment -- it was right all along, and what was inconsistent
    was the setup section it was paired with;
  * the rate constant is 4/pi in that metric and 2/pi on the unit sphere,
    the two differing by the metric scale (a distance moment is homogeneous
    of degree one in the metric).

Also pinned, both halves of the M1 label question.  POSITIVE: the constant IS
a quotient of sphere volumes -- Vol(S^2)/Vol(S^3) = 2/pi exactly at unit radius,
twice that in the dual-Coxeter metric -- so M1's volume content and its
Q[pi, 1/pi] period ring stand.  NEGATIVE: it is not the Hopf fibration's
base-to-total ratio, which is 1/(2 pi) at unit radius and 1/(4 pi) in the
dual-Coxeter metric (the Hopf map is a Riemannian submersion S^3(r) ->
S^2(r/2)); the printed identity 4/pi = Vol(S^2)/pi^2 pairs the dual-Coxeter
base volume with half the unit-metric group volume.  So the name attached the
right kind of object to the wrong sphere.

CONVENTION SCOPE (DELTA #4, 2026-09-03).  The rule Cas(ad) = h^v is the
corpus's declared one, not the field-standard one: Kac's basic form
(theta|theta) = 2 gives Cas(ad) = 2 h^v, radius sqrt(2) and constant
2 sqrt(2)/pi.  4/pi is the value in the declared rule, not a canonical number.
"""
from __future__ import annotations

import mpmath
import numpy as np
import pytest
import sympy as sp

from geovac.central_fejer_su2 import _chi, central_fejer_kernel_su2, gamma_rate


def unit_s3_moment(n: int, prec: int = 30) -> mpmath.mpf:
    """integral_{SU(2)} K_n(g) d_round(e, g) dg with d_round = theta in [0, pi]
    (unit S^3) and Haar class measure (2/pi) sin^2(theta) dtheta; the kernel is
    evaluated at the rotation angle chi = 2 theta."""
    mpmath.mp.dps = prec
    K = sp.lambdify(_chi, central_fejer_kernel_su2(n, _chi), modules="mpmath")
    return mpmath.quad(lambda th: K(2 * th) * th * (2 / mpmath.pi) * mpmath.sin(th) ** 2, [0, mpmath.pi])


def test_rotation_angle_is_twice_the_geodesic_distance():
    th = 0.7
    g = np.cos(th) * np.eye(2) + 1j * np.sin(th) * np.array([[0, 1], [1, 0]])
    phases = np.sort(np.angle(np.linalg.eigvals(g)))
    assert np.allclose(phases, [-th, th])                 # eigenphases +-theta
    chi = sp.Symbol("chi")
    char_half = sp.sin(2 * chi / 2) / sp.sin(chi / 2)      # module's spin-1/2 character
    assert abs(float(char_half.subs(chi, 2 * th)) - np.trace(g).real) < 1e-12


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_module_gamma_is_twice_the_unit_s3_moment(n):
    ratio = gamma_rate(n, prec=30) / unit_s3_moment(n)
    assert abs(ratio - 2) < 1e-12, ratio


def test_gamma_1_is_the_dual_coxeter_haar_mean_distance():
    """K_1 == 1, so gamma_1 is the Haar-mean distance: pi/2 on the unit S^3
    (pi is the diameter, which no mean distance can reach)."""
    mpmath.mp.dps = 30
    assert sp.simplify(central_fejer_kernel_su2(1, _chi)) == 1
    assert abs(gamma_rate(1, prec=30) - mpmath.pi) < 1e-20
    assert abs(unit_s3_moment(1) - mpmath.pi / 2) < 1e-20


def test_rate_constant_by_convention():
    """Doubling estimator a_n = (2n g_{2n} - n g_n)/log 2 at n = 25 in both
    conventions: the unit-S^3 value is exactly half the module's."""
    g25, g50 = gamma_rate(25, prec=30), gamma_rate(50, prec=30)
    a_rot = (2 * 25 * g50 - 25 * g25) / mpmath.log(2)
    u25, u50 = unit_s3_moment(25), unit_s3_moment(50)
    a_unit = (2 * 25 * u50 - 25 * u25) / mpmath.log(2)
    assert abs(a_rot / a_unit - 2) < 1e-10
    # the estimators sit above their limits 4/pi and 2/pi at this n (approach from above)
    assert a_rot > 4 / mpmath.pi and a_unit > 2 / mpmath.pi


# ---------------------------------------------------------------------------
# The dual-Coxeter derivation (added 2026-09-03, PI Option 1)
# ---------------------------------------------------------------------------

SIGMA = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]], complex)]


def _adjoint_matrices():
    """ad(X_a) on su(2) in the basis X_a = -i sigma_a / 2, which satisfies
    [X_a, X_b] = eps_abc X_c."""
    X = [-1j * s / 2 for s in SIGMA]
    ad = np.zeros((3, 3, 3))
    for a in range(3):
        for b in range(3):
            C = X[a] @ X[b] - X[b] @ X[a]
            for c in range(3):
                ad[a][c][b] = np.real(np.trace(C @ np.conj(X[c]).T)
                                      / np.trace(X[c] @ np.conj(X[c]).T))
    return X, ad


def _su3_basis() -> list[np.ndarray]:
    """X_a = -i lambda_a / 2, the anti-Hermitian su(3) generators."""
    l = [np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
         np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]]),
         np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
         np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
         np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]]),
         np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
         np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]),
         np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3)]
    return [-1j * m / 2 for m in l]


def _ad(basis: list[np.ndarray]) -> list[np.ndarray]:
    """ad(X_a) in the given basis, recovered by the (orthogonal) trace form."""
    d = len(basis)
    nrm = [np.trace(Xc @ np.conj(Xc).T) for Xc in basis]
    out = []
    for a in range(d):
        M = np.zeros((d, d), dtype=complex)
        for b in range(d):
            C = basis[a] @ basis[b] - basis[b] @ basis[a]
            for c in range(d):
                M[c, b] = np.trace(C @ np.conj(basis[c]).T) / nrm[c]
        out.append(np.real(M))
    return out


def _dual_coxeter_number(basis: list[np.ndarray]) -> float:
    r"""h^v computed, not asserted: the Killing form of su(N) is
    K(X, Y) = 2 h^v (X|Y) with the basic form (X|Y) = tr_def(XY), so
    h^v = K(X, X) / (2 tr_def(X X)) on any generator."""
    ad = _ad(basis)
    return float(np.trace(ad[0] @ ad[0]).real
                 / (2 * np.trace(basis[0] @ basis[0]).real))


@pytest.mark.parametrize("name, basis_fn, h_vee", [
    ("su(2)", lambda: [-1j * s / 2 for s in SIGMA], 2.0),
    ("su(3)", _su3_basis, 3.0),
])
def test_dual_coxeter_number_is_computed_not_assumed(name, basis_fn, h_vee):
    """The h^v that the rule Cas(ad) = h^v refers to, derived from the Killing
    form rather than looked up (DELTA #4 F1b: it used to be a literal, so the
    'forces lambda = 2' step was 4/2 with the 2 written in by hand)."""
    assert abs(_dual_coxeter_number(basis_fn()) - h_vee) < 1e-10


def test_dual_coxeter_rule_fixes_the_inner_product():
    """Cas(ad) = h^v forces <X,Y> = -2 tr(XY) on su(2).  h^v is computed here,
    not assumed (F1b), and the tr(XY) side is exercised on both candidate
    normalisations rather than asserted."""
    X, ad = _adjoint_matrices()
    S = sum(ad[a] @ ad[a] for a in range(3))
    assert np.allclose(S, -2 * np.eye(3))          # sum_a ad(X_a)^2 = -2 I
    h_vee = _dual_coxeter_number([-1j * s / 2 for s in SIGMA])
    # with <X,Y> = -lam tr(XY) the dual basis is (2/lam) X_a, so |Cas| = 4/lam
    for lam, expected in ((1.0, 4.0), (2.0, 2.0), (4.0, 1.0)):
        assert abs(-(2 / lam) * S[0][0] - expected) < 1e-12
        # ... and the inner product it names really is -lam tr(XY)
        assert abs(-lam * np.trace(X[0] @ X[0]).real - lam / 2) < 1e-12
    lam_dual_coxeter = 4.0 / h_vee                  # |Cas(ad)| = h^v
    assert abs(lam_dual_coxeter - 2.0) < 1e-10
    # the FIELD-STANDARD rule is a different one: (theta|theta) = 2 gives
    # Cas(ad) = 2 h^v, hence lam = 1 and the sphere of radius sqrt(2).
    assert abs(4.0 / (2 * h_vee) - 1.0) < 1e-10


@pytest.mark.parametrize("theta", [0.3, 1.0, np.pi / 2, np.pi - 0.1])
def test_dual_coxeter_distance_is_the_rotation_angle(theta):
    """d(e, exp(i theta n.sigma)) = 2 theta = chi in the dual-Coxeter metric,
    while the unit-S^3 distance is theta."""
    g_gen = 1j * theta * SIGMA[2]                   # in su(2)
    d_dc = np.sqrt(-2.0 * np.trace(g_gen @ g_gen).real)
    assert abs(d_dc - 2 * theta) < 1e-12
    eigphases = np.sort(np.angle(np.linalg.eigvals(
        np.cos(theta) * np.eye(2) + 1j * np.sin(theta) * SIGMA[2])))
    assert np.allclose(eigphases, [-theta, theta])  # unit-S^3 distance is theta


def test_dual_coxeter_sphere_has_radius_two():
    """r = 2, derived from the metric rather than restated (DELTA #4 F8: the
    previous body was an identity on its own literals and passed unchanged when
    retargeted to 'radius three').

    Route: the one-parameter subgroup t -> exp(i t n.sigma) closes at t = 2 pi
    (exp is 2 pi-periodic in t), and its tangent X = i n.sigma has dual-Coxeter
    norm ||X|| = sqrt(-2 tr X^2) = 2.  So that closed geodesic has length
    2 * 2 pi = 4 pi = 2 pi r, giving r = 2; the diameter pi r and the volume
    2 pi^2 r^3 then follow."""
    n_hat = np.array([0.6, 0.0, 0.8])                # any unit vector
    Xgen = 1j * sum(n_hat[k] * SIGMA[k] for k in range(3))
    speed = np.sqrt(-2.0 * np.trace(Xgen @ Xgen).real)
    assert abs(speed - 2.0) < 1e-12
    n_sigma = sum(n_hat[k] * SIGMA[k] for k in range(3))
    period = 2 * np.pi                               # exp(i t n.sigma) is 2pi-periodic
    closed = np.cos(period) * np.eye(2) + 1j * np.sin(period) * n_sigma
    assert np.allclose(closed, np.eye(2), atol=1e-12)
    half = np.cos(period / 2) * np.eye(2) + 1j * np.sin(period / 2) * n_sigma
    assert not np.allclose(half, np.eye(2), atol=1e-6)   # 2pi is the FIRST period
    circumference = speed * period
    radius = circumference / (2 * np.pi)
    assert abs(radius - 2.0) < 1e-12
    assert abs(np.pi * radius - 2 * np.pi) < 1e-12               # diameter
    assert abs(2 * np.pi ** 2 * radius ** 3 - 16 * np.pi ** 2) < 1e-9   # volume


def test_constant_is_the_unit_sphere_volume_quotient():
    """The POSITIVE half (2026-09-03): Vol(S^2)/Vol(S^3) = 2/pi exactly for
    unit-radius spheres, and the dual-Coxeter constant is twice that -- so the
    constant carries genuine volume content (Paper 18's M1 slot, period ring
    Q[pi, 1/pi]).  Only the 'Hopf base' label is wrong."""
    vol_s2 = 4 * np.pi
    vol_s3 = 2 * np.pi ** 2
    assert abs(vol_s2 / vol_s3 - 2 / np.pi) < 1e-12
    assert abs(2 * vol_s2 / vol_s3 - 4 / np.pi) < 1e-12


def test_hopf_base_ratio_is_not_the_rate_constant():
    """Guard on the retired reading: Vol(base)/Vol(group) is 1/(4 pi) in the
    dual-Coxeter metric and 1/(2 pi) on the unit sphere -- neither is 4/pi or
    2/pi.  The printed identity 4/pi = Vol(S^2)/pi^2 mixes normalisations."""
    for r, expected_ratio in ((1.0, 1 / (2 * np.pi)), (2.0, 1 / (4 * np.pi))):
        vol_s3 = 2 * np.pi ** 2 * r ** 3
        # the base radius is DERIVED, not assumed: the Hopf map is a
        # Riemannian submersion with fibre length 2 pi r, so
        # Vol(S^3(r)) = Vol(S^2(rho)) * 2 pi r forces rho = r/2.
        rho = np.sqrt(vol_s3 / (2 * np.pi * r) / (4 * np.pi))
        assert abs(rho - r / 2) < 1e-12
        vol_base = 4 * np.pi * rho ** 2              # Hopf base S^2(r/2)
        assert abs(vol_base / vol_s3 - expected_ratio) < 1e-12
        assert abs(vol_base / vol_s3 - 4 / np.pi) > 1.0
        assert abs(vol_base / vol_s3 - 2 / np.pi) > 0.4
    # the identity as printed: dual-Coxeter base volume over half the UNIT
    # group volume
    assert abs(4 * np.pi / np.pi ** 2 - 4 / np.pi) < 1e-12


# --- kernel sensitivity, and the seminorm normalisation ---------------------
# DELTA #4 F5A/B: `test_module_gamma_is_twice_the_unit_s3_moment` pins a ratio
# that a change of variables makes 2 for ANY kernel, so it cannot see a kernel
# regression.  It is a convention pin (SYMBOLIC), and this is the companion
# that is actually sensitive to which kernel `gamma_rate` integrates.

@pytest.mark.parametrize("n, expected", [(1, 3.14159265358979),
                                         (3, 1.61005996806574),
                                         (5, 1.13021895482038)])
def test_gamma_rate_values_are_kernel_sensitive(n, expected):
    """MEASURED values of the central-Fejer first moment.  Unlike the ratio
    pin, these move if the kernel, its cutoff, or the class measure changes --
    verified by planting a Cesaro kernel, which the ratio pin does not see."""
    assert abs(float(gamma_rate(n)) - expected) < 1e-9


def test_seminorm_normalisation_is_the_metric_scale():
    """Paper 38 eq:seminorm_normalisation, L(f) = (1/2) ||[D_CH, M_f]||.

    Content: ||[D_CH, M_f]|| is the UNIT-metric Lipschitz seminorm (it is the
    sup of the unit-metric gradient), and the dual-Coxeter distance is twice the
    unit one, so the dual-Coxeter Lipschitz seminorm is half of it.  Checked on
    sampled pairs of a concrete function, not asserted."""
    rng = np.random.default_rng(38)
    # f(g) = Re tr(g) on SU(2), i.e. 2 cos(theta) in terms of the unit-S^3
    # distance theta from the identity.
    def theta_of(q):                                  # q a unit quaternion
        return float(np.arccos(np.clip(q[0], -1.0, 1.0)))
    qs = rng.normal(size=(400, 4))
    qs /= np.linalg.norm(qs, axis=1, keepdims=True)
    best_unit, best_dc = 0.0, 0.0
    for i in range(len(qs)):
        for j in range(i + 1, len(qs)):
            # relative rotation q_i^{-1} q_j -> its unit-S^3 distance
            qi, qj = qs[i], qs[j]
            dot = abs(float(np.dot(qi, qj)))
            d_unit = float(np.arccos(np.clip(dot, -1.0, 1.0)))
            if d_unit < 1e-6:
                continue
            df = abs(2 * np.cos(theta_of(qi)) - 2 * np.cos(theta_of(qj)))
            best_unit = max(best_unit, df / d_unit)
            best_dc = max(best_dc, df / (2 * d_unit))
    assert best_unit > 1.0                            # non-degenerate sample
    assert abs(best_dc - best_unit / 2) < 1e-12       # L = (1/2) Lipnorm
