"""Paper 53, Sec. 'The finite disk-with-cone: measurements' (prop:disk, prop:prop2).

REBUILT 2026-09-29.  Earlier versions of this file asserted two claims that are
now WITHDRAWN in the paper:

  * "the Cesaro weight s>=2 PRESERVES positivity on the finite disk" -- it was
    checked on ONE observable, (rho/R)^2.  A non-negative ring observable is
    mapped to a function with negative values (prop:disk (b)).
  * "the interior approximation rate does NOT decay (Lambda^{+0.07}, plateau
    ~0.3)" -- the plateau came from evaluating the 0/0 Markov ratio at the exact
    grid point rho = R, where every J_0(a_j) is ~1e-16 rounding noise.  On the
    open disk the error DECAYS (prop:disk (c), rem:interior_correction).

What this file asserts now, each against a named wrong answer:

  (1) test_disk_cesaro_not_positivity_preserving -- a non-negative ring f has
      min g < 0 at s = 2 (and s = 4) while the denominator stays positive.
      Rejects: "Cesaro s>=2 is positivity-preserving on the disk".
  (2) test_disk_rho2_positive_and_heat_negative -- on (rho/R)^2 the Cesaro g
      stays positive (the single-observable fact the old claim rested on) and
      the heat weight dips below -0.10 at Lambda = 200.  Rejects: a drift of
      those two measured values.
  (3) test_disk_open_disk_error_decays -- sup over rho < R of |g - f|/|grad f|
      decays (fitted exponent < -0.3; last < half of first) for (rho/R)^2 and
      for an interior bump.  Rejects: "the finite-disk error does not decay"
      (the withdrawn Lambda^{+0.07} obstruction).
  (4) test_disk_plateau_is_the_boundary_point_artifact -- including the exact
      point rho = R reproduces the old plateau (no decay), and moving that one
      point to R - 1e-9 removes it.  Pins the diagnosis of the retraction.

(5) prop:prop2 -- O^2 fills M_N(C) at four truncations (propagation number 2),
    via the live geovac.operator_system machinery.  [MEASURED], not all cutoffs.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.special import jn_zeros, jv

from geovac.operator_system import operator_system_dim, _extract_matrix_basis

R = 1.0


# --------------------------------------------------------------------------
# Finite-disk Markov-Cesaro map, radial (k=0) sector.
# --------------------------------------------------------------------------
def _k0_modes(Jmax: int, Lam: float):
    """k=0 Dirichlet radial modes R_j(rho)=sqrt2/(R|J_1(a_j)|) J_0(a_j rho/R),
    normalized int_0^R R_j^2 rho drho = 1; eigenvalue (a_j/R)^2 <= Lambda."""
    modes = []
    for a in jn_zeros(0, Jmax):
        lam = (a / R) ** 2
        if lam <= Lam:
            modes.append((a, lam, np.sqrt(2.0) / (R * abs(jv(1, a)))))
    return modes


def _ratio_parts(f_fn, modes, weight, rho, nquad: int = 1200):
    """Numerator and denominator of the Markov ratio g = num/den."""
    xg, wg = leggauss(nquad)
    s = 0.5 * R * (xg + 1.0)
    ws = 0.5 * R * wg
    num = np.zeros_like(rho, dtype=float)
    den = np.zeros_like(rho, dtype=float)
    fs = f_fn(s)
    for (a, lam, norm) in modes:
        Rj_s = norm * jv(0, a * s / R)
        cf = np.sum(Rj_s * fs * s * ws)
        c1 = np.sum(Rj_s * s * ws)
        w = weight(lam)
        Rj_rho = norm * jv(0, a * rho / R)
        num += w * cf * Rj_rho
        den += w * c1 * Rj_rho
    return num, den


def _smooth_ratio(f_fn, modes, weight, rho):
    num, den = _ratio_parts(f_fn, modes, weight, rho)
    return num / den


def _heat_weight(Lam):
    return lambda lam: np.exp(-(1.0 / Lam) * lam)


def _cesaro_weight(Lam, s=2):
    return lambda lam: max(0.0, (1.0 - lam / Lam)) ** s


def _ring(rho0: float, width: float = 0.02):
    return lambda r: np.exp(-(((r - rho0) / width) ** 2))


def _most_negative_kernel_radius(modes, weight) -> float:
    """sigma in (0, 0.9] where K(rho ~ 0, sigma) = sum_a w_a R_a(0) R_a(sigma)
    is most negative (the ring is placed there)."""
    sig = np.linspace(1e-3, 0.9, 900)
    K = np.zeros_like(sig)
    for (a, lam, norm) in modes:
        K += weight(lam) * norm * jv(0, 0.0) * norm * jv(0, a * sig / R)
    return float(sig[np.argmin(K)])


_OPEN_RHO = np.linspace(1e-4, R - 1e-6, 3000)


@pytest.mark.parametrize("s,Lams", [(2, (100.0, 400.0, 1600.0)), (4, (400.0, 1600.0))])
def test_disk_cesaro_not_positivity_preserving(s, Lams):
    """prop:disk (b): a non-negative ring is sent to a function with negative values."""
    for Lam in Lams:
        modes = _k0_modes(400, Lam)
        w = _cesaro_weight(Lam, s=s)
        rho0 = _most_negative_kernel_radius(modes, w)
        num, den = _ratio_parts(_ring(rho0), modes, w, _OPEN_RHO[_OPEN_RHO <= 0.9])
        g = num / den
        i = int(np.argmin(g))
        assert den[i] > 0.1, f"s={s} Lam={Lam}: denominator not positive at the min ({den[i]})"
        assert g[i] < -1e-3, (
            f"s={s} Lam={Lam}: ring f>=0 gave min g = {g[i]:.5f}; a positive map "
            "would give g >= 0 -- the paper says the Cesaro map is NOT positive"
        )


def test_disk_rho2_positive_and_heat_negative():
    """prop:disk (b): on (rho/R)^2 the Cesaro g stays positive; heat dips < -0.10."""
    f1 = lambda r: (r / R) ** 2
    ces, heat = [], []
    for Lam in (100.0, 200.0, 400.0, 800.0, 1600.0):
        modes = _k0_modes(80, Lam)
        ces.append(float(_smooth_ratio(f1, modes, _cesaro_weight(Lam, 2), _OPEN_RHO).min()))
        heat.append(float(_smooth_ratio(f1, modes, _heat_weight(Lam), _OPEN_RHO).min()))
    assert min(ces) > 0.0, f"Cesaro on (rho/R)^2: mins {ces}"
    assert heat[1] < -0.10, f"heat Lambda=200 min g = {heat[1]:.4f}, expected < -0.10"


def _open_disk_errors(f_fn, grad, Lams, rho):
    errs = []
    for Lam in Lams:
        modes = _k0_modes(400, Lam)
        g = _smooth_ratio(f_fn, modes, _cesaro_weight(Lam, 2), rho)
        errs.append(float(np.max(np.abs(g - f_fn(rho))) / grad))
    return errs


_BUMP = lambda r: np.where(r < 0.6, np.cos(np.pi * r / 1.2) ** 2, 0.0)

_DECAY_CASES = {
    "rho2": (lambda r: (r / R) ** 2, 2.0, (100.0, 200.0, 400.0, 800.0, 1600.0)),
    "interior_bump": (_BUMP, np.pi / 1.2, (100.0, 400.0, 1600.0, 6400.0)),
}


@pytest.mark.parametrize("name", list(_DECAY_CASES))
def test_disk_open_disk_error_decays(name):
    """prop:disk (c): on the open disk the reconstruction error DECAYS."""
    f_fn, grad, Lams = _DECAY_CASES[name]
    errs = _open_disk_errors(f_fn, grad, Lams, _OPEN_RHO)
    slope = float(np.polyfit(np.log(Lams), np.log(errs), 1)[0])
    assert slope < -0.3, f"{name}: open-disk error must decay; slope {slope:.3f}, errs {errs}"
    assert errs[-1] < 0.5 * errs[0], f"{name}: errs {errs}"


def test_disk_plateau_is_the_boundary_point_artifact():
    """rem:interior_correction: the old plateau is the single point rho = R."""
    f1 = lambda r: (r / R) ** 2
    Lams = (100.0, 400.0, 1600.0)
    with_R = np.linspace(1e-4, R, 2000)          # the old grid, endpoint included
    moved = with_R.copy()
    moved[-1] = R - 1e-9                          # the same grid, one point moved
    e_with, e_moved = [], []
    for Lam in Lams:
        modes = _k0_modes(80, Lam)
        w = _cesaro_weight(Lam, 2)
        e_with.append(float(np.max(np.abs(_smooth_ratio(f1, modes, w, with_R) - f1(with_R))) / 2.0))
        e_moved.append(float(np.max(np.abs(_smooth_ratio(f1, modes, w, moved) - f1(moved))) / 2.0))
    # old grid: no decay (the withdrawn plateau)
    assert min(e_with) > 0.15, f"old-grid plateau not reproduced: {e_with}"
    # one point moved: decays
    assert e_moved[-1] < 0.5 * e_moved[0], f"moved-point errors do not decay: {e_moved}"
    assert e_moved[-1] < 0.1, f"moved-point error at Lambda=1600: {e_moved[-1]}"


# --------------------------------------------------------------------------
# (2) prop:prop2 — assembled disk operator system O_n = P_n C(disk) P_n has
#     propagation number 2 (O² fills the full M_N(C) envelope).
# --------------------------------------------------------------------------
def _radial_modes(m: float, R_: float, N_rho: int, J: int):
    from scipy.linalg import eigh_tridiagonal
    a = R_ / N_rho
    k = np.arange(1, N_rho + 1)
    rho = k * a
    diag = 2.0 / a**2 + (m**2 - 0.25) / rho**2
    off = -np.ones(N_rho - 1) / a**2
    _, v = eigh_tridiagonal(diag, off, select="i", select_range=(0, J - 1))
    return v, rho, a


def _build_disk_operator_system(K: int, J: int, R_: float, N_rho: int, max_deg: int):
    ks = list(range(-K, K))           # 2K azimuthal modes
    n_modes = len(ks) * J
    mode_vecs, rho, a = {}, None, None
    for k in ks:
        v, rho, a = _radial_modes(abs(k + 0.5), R_, N_rho, J)
        mode_vecs[k] = v

    def idx(k, j):
        return (k + K) * J + j

    def mult_matrix(p: int, q: int):
        M = np.zeros((n_modes, n_modes), dtype=complex)
        for k in ks:
            kp = k + q
            if kp not in ks:
                continue
            wgt = (rho ** p) * a
            Ovl = (mode_vecs[kp] * wgt[:, None]).T @ mode_vecs[k]
            for jp in range(J):
                for j in range(J):
                    M[idx(kp, jp), idx(k, j)] = Ovl[jp, j]
        return M

    gens = [np.eye(n_modes, dtype=complex)]
    seen = set()
    for a_ in range(max_deg + 1):
        for b_ in range(max_deg + 1):
            if a_ + b_ == 0 or a_ + b_ > max_deg:
                continue
            q, p = a_ - b_, a_ + b_
            if abs(q) >= 2 * K or (p, q) in seen:
                continue
            seen.add((p, q))
            gens.append(mult_matrix(p, q))
    return gens, n_modes


def _prop_from_generators(gens, max_k=4, tol=1e-9):
    N = gens[0].shape[0]
    target = N * N
    cur = list(gens)
    dims = [operator_system_dim(cur, tol=tol)]
    if dims[0] == target:
        return 1, dims
    for k in range(2, max_k + 1):
        basis = _extract_matrix_basis(cur, tol=tol)
        newg = [A @ B for A in basis for B in gens]
        dk = operator_system_dim(newg, tol=tol)
        dims.append(dk)
        if dk == target:
            return k, dims
        if dk == dims[-2]:
            return -1, dims
        cur = newg
    return -1, dims


@pytest.mark.parametrize("K,J,N", [(2, 2, 8), (2, 3, 12), (3, 2, 12), (3, 3, 18)])
def test_disk_operator_system_prop2(K, J, N):
    """prop:prop2 — O² already spans the full M_N(C) envelope (propagation 2)."""
    gens, n_modes = _build_disk_operator_system(K, J, R_=8.0, N_rho=120, max_deg=8)
    assert n_modes == N, f"(K={K},J={J}) expected N={N}, got {n_modes}"
    prop, dims = _prop_from_generators(gens)
    assert prop == 2, f"(K={K},J={J}) expected prop=2; got {prop} (dims={dims})"
    # the 2-step generating signature: O² == full M_N(C) envelope (N²)
    assert dims[1] == N * N, (
        f"(K={K},J={J}) O² should fill M_{N}(C) envelope {N*N}; got {dims[1]}"
    )


# --------------------------------------------------------------------------
# Guard pass (/qa trunk DELTA #17): the printed prop:disk (b)/(c) numbers,
# the interior Cesaro kernel ratios, and the operator-level non-positivity
# of P M_g P, each pinned to the precision the paper prints.
# --------------------------------------------------------------------------
def _ring_ratio_min(s: int, Lam: float) -> float:
    modes = _k0_modes(400, Lam)
    w = _cesaro_weight(Lam, s=s)
    rho0 = _most_negative_kernel_radius(modes, w)
    num, den = _ratio_parts(_ring(rho0), modes, w, _OPEN_RHO[_OPEN_RHO <= 0.9])
    return float((num / den).min())


# printed in prop:disk (b): min g on the ring
_PRINTED_RING_MIN = {(2, 100.0): -0.016, (2, 400.0): -0.025, (2, 1600.0): -0.045,
                     (4, 400.0): -0.006, (4, 1600.0): -0.011}


@pytest.mark.parametrize("key", list(_PRINTED_RING_MIN))
def test_disk_printed_ring_minima(key):
    """prop:disk (b) printed min g values, to 3 decimals.  Rejects a drift of any
    printed value (e.g. the earlier -0.044 at s=2, Lambda=1600; measured -0.0448)."""
    s, Lam = key
    assert round(_ring_ratio_min(s, Lam), 3) == _PRINTED_RING_MIN[key]


def test_disk_printed_heat_and_rho2_errors():
    """prop:disk (b) heat min ~ -0.13 at Lambda=200; prop:disk (c) errors
    0.237, 0.179, 0.133, 0.097, 0.070 (exponent -0.44) and bump 0.163 -> 0.0032
    (exponent roughly -0.9)."""
    f1 = lambda r: (r / R) ** 2
    heat = float(_smooth_ratio(f1, _k0_modes(80, 200.0), _heat_weight(200.0), _OPEN_RHO).min())
    assert round(heat, 2) == -0.13
    f_fn, grad, Lams = _DECAY_CASES["rho2"]
    errs = _open_disk_errors(f_fn, grad, Lams, _OPEN_RHO)
    assert [round(e, 3) for e in errs] == [0.237, 0.179, 0.133, 0.097, 0.070]
    assert round(float(np.polyfit(np.log(Lams), np.log(errs), 1)[0]), 2) == -0.44
    f_fn, grad, Lams = _DECAY_CASES["interior_bump"]
    errs = _open_disk_errors(f_fn, grad, Lams, _OPEN_RHO)
    assert round(errs[0], 3) == 0.163 and round(errs[-1], 4) == 0.0032
    assert abs(float(np.polyfit(np.log(Lams), np.log(errs), 1)[0]) + 0.9) < 0.1


@pytest.mark.parametrize("s,printed", [(2, -0.03), (4, -0.009)])
def test_disk_interior_cesaro_kernel_ratio(s, printed):
    """prop:disk (b): interior Cesaro kernel K(0, sigma) has min/max ~ -0.03 (s=2),
    -0.009 (s=4) -- the plane Bochner-Riesz values.  Rejects: a non-negative
    kernel, or a drift of the printed ratios."""
    Lam = 1600.0
    modes = _k0_modes(800, Lam)
    w = _cesaro_weight(Lam, s=s)
    sig = np.linspace(1e-4, 0.9, 20000)
    K = sum(w(lam) * nrm * nrm * jv(0, a * sig / R) for (a, lam, nrm) in modes)
    ratio = float(K.min() / K.max())
    assert ratio < 0.0
    d = 2 if s == 2 else 3
    assert round(ratio, d) == printed, f"s={s}: min/max {ratio:.5f}"


def _pmgp_min_eig(s: int, Lam: float, nquad: int = 1200) -> float:
    """Smallest eigenvalue of P M_g P in the k=0 sector, g = Markov-Cesaro
    reconstruction of the non-negative ring."""
    modes = _k0_modes(400, Lam)
    w = _cesaro_weight(Lam, s=s)
    rho0 = _most_negative_kernel_radius(modes, w)
    xg, wg = leggauss(nquad)
    x = 0.5 * R * (xg + 1.0)
    wx = 0.5 * R * wg
    num, den = _ratio_parts(_ring(rho0), modes, w, x)
    g = num / den
    B = np.array([nrm * jv(0, a * x / R) for (a, lam, nrm) in modes])
    M = (B * (g * x * wx)) @ B.T
    return float(np.linalg.eigvalsh(M)[0])


# printed in prop:disk (b): smallest eigenvalue of the compressed operator
_PRINTED_PMGP = {(2, 100.0): -0.0034, (2, 400.0): -0.0022, (2, 1600.0): -0.0022,
                 (4, 1600.0): -0.00044, (4, 400.0): 0.0015}


@pytest.mark.parametrize("key", list(_PRINTED_PMGP))
def test_disk_compressed_operator_not_positive(key):
    """prop:disk (b), operator level: B_Lambda(f) = P M_g P has a negative
    eigenvalue for a non-negative f (except s=4, Lambda=400, printed positive).
    Rejects: 'the failure is pointwise only; the compressed operator stays
    positive'."""
    s, Lam = key
    lam_min = _pmgp_min_eig(s, Lam)
    printed = _PRINTED_PMGP[key]
    d = 5 if abs(printed) < 1e-3 else 4
    assert round(lam_min, d) == printed, f"s={s} Lam={Lam}: {lam_min:.6f}"
    if printed < 0:
        assert lam_min < 0.0
    # quadrature-converged: doubling the nodes moves it by < 1e-9
    assert abs(_pmgp_min_eig(s, Lam, nquad=2400) - lam_min) < 1e-9
