"""Backing tests for the Paper 32 remark rem:config_berry_curvature: the configuration operator
F = sum_i P_i carries a FLAT Z2 line bundle over the nuclear-geometry manifold, whose curvature
is a sum of pi-flux deltas at the conical intersections (Longuet-Higgins / molecular Aharonov-
Bohm structure, for F rather than the electronic Hamiltonian).

The overlap engine (fast prolate-spheroidal Gauss-Laguerre x Gauss-Legendre) is embedded here
self-contained and is the SAME evaluator validated to 2.3e-14 vs mpmath adaptive quadrature
(debug/beh2_ci_exact_landscape.py; permanent record CHANGELOG v5.1.0). The load-bearing logic
tested is the geometric-phase STRUCTURE (Z2 quantization, a pi-source at a known CI + flatness
away, the loop law, and the real-monopole 2x2 model), computed on the validated overlaps.

Claims pinned:
  (1) F is real-symmetric => every plaquette Berry flux is exactly 0 or pi (Z2), machine-exact;
  (2) the crossing band carries pi flux at the central conical intersection (2.445,2.445) and
      0 on a CI-free plaquette (the curvature is delta-supported at the CIs, flat elsewhere);
  (3) loop law: parallel transport gives phase pi around a loop enclosing the CI and 0 around a
      CI-free loop (Z2 holonomy = enclosed-CI count mod 2);
  (4) each CI is the real (equatorial) section of a charge-1/2 Berry monopole: the local 2x2
      model has sigma_y coefficient == 0 (real => equator) and nonzero branching Jacobian;
  (5) U(1) LIFT (slow): breaking time-reversal with a magnetic (Peierls) phase phi makes F
      complex-Hermitian and opens the sigma_y axis; the first Chern number over a sphere
      enclosing the degeneracy is +-1 (opposite on the two crossing bands) -- an integer Berry
      monopole. Reality (time-reversal), NOT planarity, is what protects the Z2.
"""
import numpy as np
import pytest
from scipy.special import eval_genlaguerre, lpmv, factorial
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss

STATES = [(1, 0), (2, 1)]                       # sigma pair 1s, 2p0 (m=0)
PAR = np.diag([1.0, -1.0])                       # 2p0 odd under z->-z
_LAG_X, _LAG_W = laggauss(64)
_LEG_X, _LEG_W = leggauss(96)


def _R_nl(Z, n, l, r):
    Z = float(Z); rho = 2 * Z * r / n
    norm = np.sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * eval_genlaguerre(n - l - 1, 2 * l + 1, rho)


def _ang(l, m, ct):
    m = abs(m)
    tn = np.sqrt((2 * l + 1) / 2.0 * factorial(l - m) / factorial(l + m))
    return tn * lpmv(m, l, ct)


def _overlap(Z1, n1, l1, Z2, n2, l2, R):
    """<phi(0,Z1)|phi(R zhat,Z2)>, m=0, exact per-geometry prolate-spheroidal quadrature."""
    half = R / 2.0
    a = half * (Z1 / n1 + Z2 / n2)
    xi = 1.0 + _LAG_X / a
    XI, ETA = np.meshgrid(xi, _LEG_X, indexing="ij")
    r1 = half * (XI + ETA); r2 = half * (XI - ETA)
    with np.errstate(divide="ignore", invalid="ignore"):
        ct1 = (1 + XI * ETA) / (XI + ETA); ct2 = (XI * ETA - 1) / (XI - ETA)
    val = (_R_nl(Z1, n1, l1, r1) * _ang(l1, 0, ct1)
           * _R_nl(Z2, n2, l2, r2) * _ang(l2, 0, ct2) * half ** 3 * (XI ** 2 - ETA ** 2))
    val = np.where((r1 <= 0) | (r2 <= 0) | ~np.isfinite(val), 0.0, val)
    wl = (_LAG_W * np.exp(_LAG_X) / a)[:, None]
    return float(np.sum(wl * (_LEG_W[None, :] * val)))


def _Sblock(R, Za, Zb):
    return np.array([[_overlap(Za, na, la, Zb, nb, lb, R) for (nb, lb) in STATES] for (na, la) in STATES])


def _projectors(d1, d2):
    S1 = _Sblock(d1, 2, 1)
    S2 = PAR @ _Sblock(d2, 2, 1) @ PAR
    SHH = PAR @ _Sblock(d1 + d2, 1, 1) @ PAR
    I = np.eye(2)
    G = np.block([[I, S1, S2], [S1.T, I, SHH], [S2.T, SHH.T, I]])
    if np.linalg.eigvalsh(G).min() < 1e-9:
        return None
    Xh = np.linalg.cholesky(G).T
    return [Xh[:, 2 * k:2 * k + 2] @ np.linalg.pinv(Xh[:, 2 * k:2 * k + 2]) for k in range(3)]


def _F_bands(d1, d2):
    Ps = _projectors(d1, d2)
    return None if Ps is None else np.linalg.eigh(Ps[0] + Ps[1] + Ps[2])


def _plaq_flux(d1, d2, h, band):
    vs = [_F_bands(*p)[1][:, band] for p in [(d1, d2), (d1 + h, d2), (d1 + h, d2 + h), (d1, d2 + h)]]
    prod = (vs[0] @ vs[1]) * (vs[1] @ vs[2]) * (vs[2] @ vs[3]) * (vs[3] @ vs[0])
    return float(np.angle(prod))


def _transport_sign(cx, cy, r, N=96):
    v0 = vp = None
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        w, V = _F_bands(cx + r * np.cos(th), cy + r * np.sin(th))
        if i == 0:
            k = int(np.argmin(np.diff(w))); v = V[:, k]; v0 = v.copy(); vp = v
        else:
            ov = V.T @ vp; j = int(np.argmax(np.abs(ov))); v = V[:, j] * np.sign(ov[j]); vp = v
    return float(np.sign(vp @ v0))


# ============================== tests ==============================
def test_flux_is_Z2_real_symmetric():
    """(1) Every plaquette flux is exactly 0 or pi (real F => Z2), to machine precision."""
    for (d1, d2, band) in [(2.40, 2.40, 2), (2.00, 2.00, 2), (2.60, 2.30, 3), (2.40, 2.40, 0)]:
        f = _plaq_flux(d1, d2, 0.09, band)
        assert min(abs(f), abs(abs(f) - np.pi)) < 1e-9      # in {0, pi}


def test_pi_source_at_central_ci_flat_away():
    """(2) The crossing band (2) carries pi flux at the central CI (2.445,2.445) and 0 far off."""
    assert abs(abs(_plaq_flux(2.40, 2.40, 0.09, 2)) - np.pi) < 1e-6    # pi at the CI
    assert abs(_plaq_flux(2.00, 2.00, 0.09, 2)) < 1e-6                 # flat (0) CI-free
    assert abs(_plaq_flux(2.40, 2.40, 0.09, 0)) < 1e-6                 # non-crossing band: no source


def test_loop_law_z2_holonomy():
    """(3) Transport phase = pi around a loop enclosing the CI, 0 around a CI-free loop."""
    assert _transport_sign(2.445, 2.445, 0.12) < 0                    # encloses 1 CI -> pi
    assert _transport_sign(2.10, 2.10, 0.08) > 0                      # encloses 0 CI -> 0


def test_ci_is_real_section_of_berry_monopole():
    """(4) Local 2x2 model at the central CI: sigma_y coeff == 0 (real => equator of the Bloch
    sphere); nonzero branching Jacobian (a_x,a_z winds once => pi). Each CI = charge-1/2 monopole
    seen edge-on; the sigma_y (reality-breaking) axis is the U(1)-lift direction."""
    c = 2.445
    w0, V0 = _F_bands(c, c)
    kc = int(np.argmin(np.diff(w0)))
    Q = V0[:, [kc, kc + 1]]
    sx = np.array([[0, 1], [1, 0]], float); sz = np.array([[1, 0], [0, -1]], float)
    sy = np.array([[0, -1j], [1j, 0]])

    def pauli(dd1, dd2):
        Ps = _projectors(c + dd1, c + dd2)
        M = Q.T @ (Ps[0] + Ps[1] + Ps[2]) @ Q
        M = M - 0.5 * np.trace(M) * np.eye(2)
        return (0.5 * np.trace(M @ sx).real, 0.5 * np.trace(M @ sy).real, 0.5 * np.trace(M @ sz).real)

    h = 1e-3
    axp, axm, ayp, aym = pauli(h, 0), pauli(-h, 0), pauli(0, h), pauli(0, -h)
    # sigma_y coefficient identically 0 (real family on the equator)
    assert max(abs(pauli(1e-3 * np.cos(t), 1e-3 * np.sin(t))[1]) for t in np.linspace(0, 2 * np.pi, 16)) < 1e-12
    g = np.array([(axp[2] - axm[2]) / (2 * h), (ayp[2] - aym[2]) / (2 * h)])
    hv = np.array([(axp[0] - axm[0]) / (2 * h), (ayp[0] - aym[0]) / (2 * h)])
    assert abs(g[0] * hv[1] - g[1] * hv[0]) > 1e-5                    # nonzero winding => pi holonomy


# ============================== (5) U(1) lift ==============================
def _projectors_phi(d1, d2, phi):
    """Complex-Hermitian G with a Peierls phase e^{i phi} on the H1<->H2 block (T-breaking)."""
    S1 = _Sblock(d1, 2, 1); S2 = PAR @ _Sblock(d2, 2, 1) @ PAR; SHH = PAR @ _Sblock(d1 + d2, 1, 1) @ PAR
    I = np.eye(2); e = np.exp(1j * phi)
    G = np.block([[I + 0j, S1 + 0j, S2 + 0j],
                  [S1.T + 0j, I + 0j, e * SHH],
                  [S2.T + 0j, np.conj(e) * SHH.T, I + 0j]])
    if np.linalg.eigvalsh(G).min() < 1e-9:
        return None
    Xh = np.linalg.cholesky(G).conj().T
    return [Xh[:, 2 * k:2 * k + 2] @ np.linalg.pinv(Xh[:, 2 * k:2 * k + 2]) for k in range(3)]


def _Fphi(d1, d2, phi):
    Ps = _projectors_phi(d1, d2, phi)
    return Ps[0] + Ps[1] + Ps[2]


def _chern(c, r_d, r_phi, band, N=30):
    ths = np.linspace(1e-3, np.pi - 1e-3, N); pss = np.linspace(0, 2 * np.pi, N, endpoint=False)
    V = np.empty((N, N, 6, 6), complex)
    for i, th in enumerate(ths):
        for j, ps in enumerate(pss):
            _, Vij = np.linalg.eigh(_Fphi(c + r_d * np.sin(th) * np.cos(ps),
                                          c + r_d * np.sin(th) * np.sin(ps), r_phi * np.cos(th)))
            V[i, j] = Vij
    tot = 0.0
    for i in range(N - 1):
        for j in range(N):
            jn = (j + 1) % N
            v = [V[i, j, :, band], V[i + 1, j, :, band], V[i + 1, jn, :, band], V[i, jn, :, band]]
            tot += np.angle(np.vdot(v[0], v[1]) * np.vdot(v[1], v[2]) * np.vdot(v[2], v[3]) * np.vdot(v[3], v[0]))
    return tot / (2 * np.pi)


@pytest.mark.slow
def test_u1_lift_is_integer_monopole():
    """(5) Breaking time-reversal with a magnetic phase phi lifts Z2 -> U(1): F(0) is exactly
    real; phi opens the gap linearly; and the first Chern number over a sphere enclosing the
    degeneracy is +-1 (opposite on the two crossing bands) -- a genuine integer Berry monopole."""
    c = 2.445
    w0, _ = np.linalg.eigh(_Fphi(c, c, 0.0))
    kc = int(np.argmin(np.diff(w0)))
    assert np.max(np.abs(_Fphi(c, c, 0.0).imag)) < 1e-12          # T-symmetric point is exactly real

    def gap(phi):
        w = np.linalg.eigvalsh(_Fphi(c, c, phi)); return float(w[kc + 1] - w[kc])
    assert gap(0.10) > 50 * gap(0.0)                              # phi opens the gap strongly
    assert abs(gap(0.10) / gap(0.05) - 2.0) < 0.1                 # ...linearly (a genuine sigma_y axis)

    kphi = gap(0.05) / 0.05
    r_d = 0.03; r_phi = r_d * 0.07 / kphi
    C_lo, C_hi = _chern(c, r_d, r_phi, kc), _chern(c, r_d, r_phi, kc + 1)
    assert abs(round(C_lo) - C_lo) < 0.1 and abs(round(C_hi) - C_hi) < 0.1   # integer Chern
    assert abs(round(C_lo)) == 1 and round(C_lo) == -round(C_hi)             # charge +-1, opposite
