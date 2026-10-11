"""Pairwise three-centre ERIs for general Slater orbitals (s, p), water geometry.

(XY|XZ) = int a(1) b(1) r12^-1 c(2) d(2), a and c on the shared centre X (origin),
b on Y, d on Z.  Expand each two-centre density about X in real harmonics:
rho(r) = sum_LM rho_LM(r) Y_LM(rhat).  The Laplace expansion of 1/r12 couples
only equal (L, M), so EXACTLY

    (XY|XZ) = sum_LM int r^2 dr rho^{XZ}_LM(r) V^{XY}_LM(r),
    V_LM = potential of the channel rho^{XY}_LM Y_LM (radial Poisson per L).

Every rho_LM is a ONE-PAIR object; the third centre enters only through which
direction each pair points.  The truncation in L is the only approximation.

Numerics: rho_LM(r) by angular quadrature about X in a frame whose pole points
at the pair's other centre B, with the substitution w = |r - B| (removes the
cusp at B: integrand smooth in w); FE-DVR radial grid with element boundaries
at every |B| (the channels kink there); Poisson per L in the same DVR.
Referee: McMurchie-Davidson over Gaussian fits (geovac.noci_engine), n_gauss
sequence.  Driver: debug/pairwise_general.py.
"""
import importlib.util
import math
import os

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln, lpmv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('walsh', os.path.join(ROOT, 'tests', 'test_walsh_ah2.py'))
_w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(_w)
fedvr_w = _w.fedvr_w

# water: O at origin, H in the xz-plane, R = 1.809 bohr, angle 104.5 deg
_HALF = math.radians(104.5) / 2
R_OH = 1.809
O = np.zeros(3)
H1 = R_OH * np.array([math.sin(_HALF), 0.0, math.cos(_HALF)])
H2 = R_OH * np.array([-math.sin(_HALF), 0.0, math.cos(_HALF)])

ANG = {'s': None, 'px': 0, 'py': 1, 'pz': 2}


def orbital(center, n, zeta, kind):
    return dict(center=np.asarray(center, float), n=n, zeta=zeta, kind=kind)


def orb_value(o, P):
    """Normalised Slater orbital N r^{n-1} e^{-zeta r} S(rhat) at points P (..., 3)."""
    d = P - o['center']
    r = np.linalg.norm(d, axis=-1)
    n, z = o['n'], o['zeta']
    N = (2 * z) ** (n + 0.5) / math.sqrt(math.factorial(2 * n))
    rad = N * r ** (n - 1) * np.exp(-z * r)
    if o['kind'] == 's':
        return rad / math.sqrt(4 * math.pi)
    comp = d[..., ANG[o['kind']]]
    # p: S = sqrt(3/4pi) x_i / r; r^{n-1} x_i / r = r^{n-2} x_i (n >= 2)
    return N * r ** (n - 2) * comp * np.exp(-z * r) * math.sqrt(3 / (4 * math.pi))


def real_Y_all(L, rhat):
    """Real orthonormal harmonics Y_LM, L <= Lmax, at directions rhat (N, 3).
    Returns (N, (L+1)^2) ordered (l, m) with m = -l..l (m<0 sine)."""
    x, y, z = rhat[:, 0], rhat[:, 1], np.clip(rhat[:, 2], -1, 1)
    phi = np.arctan2(y, x)
    out = []
    for l in range(L + 1):
        for m in range(-l, l + 1):
            am = abs(m)
            norm = math.sqrt((2 * l + 1) / (4 * math.pi) * math.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
            P = norm * lpmv(am, l, z) * (-1) ** am
            if m == 0:
                out.append(P)
            elif m > 0:
                out.append(math.sqrt(2) * P * np.cos(am * phi))
            else:
                out.append(math.sqrt(2) * P * np.sin(am * phi))
    return np.array(out).T


def channels(a, b, r, L, nw=48, nphi=None):
    """rho_LM(r_i) of the density a(r) b(r); a on the origin, b at B."""
    B = b['center']
    D = float(np.linalg.norm(B))
    ez = B / D
    tmp = np.array([1.0, 0.0, 0.0]) if abs(ez[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = tmp - (tmp @ ez) * ez
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(ez, e1)
    nphi = nphi or (2 * L + 8)
    phi = 2 * np.pi * np.arange(nphi) / nphi
    xg, wg = leggauss(nw)
    out = np.zeros((len(r), (L + 1) ** 2))
    for i, ri in enumerate(r):
        lo, hi = abs(ri - D), ri + D
        w = lo + (hi - lo) * (xg + 1) / 2
        ww = wg * (hi - lo) / 2 * w / (ri * D)          # du = w dw / (r D)
        u = np.clip((ri ** 2 + D ** 2 - w ** 2) / (2 * ri * D), -1, 1)
        st = np.sqrt(1 - u ** 2)
        dirs = (st[:, None, None] * (np.cos(phi)[None, :, None] * e1 + np.sin(phi)[None, :, None] * e2)
                + u[:, None, None] * ez).reshape(-1, 3)
        wt = (ww[:, None] * np.full(nphi, 2 * np.pi / nphi)[None, :]).ravel()
        P = ri * dirs
        f = orb_value(a, P) * orb_value(b, P) * wt
        out[i] = f @ real_Y_all(L, dirs)
    return out


def coulomb_channels(rho, r, T, W, rmax, L):
    """V_LM(r_i) from rho_LM(r_i): (2T + L(L+1)/r^2) U = 4 pi r rho, U = r V."""
    V = np.zeros_like(rho)
    col = 0
    for l in range(L + 1):
        A = 2 * T + np.diag(l * (l + 1) / r ** 2)
        Ainv_rhs = np.linalg.solve(A, (np.sqrt(W) * 4 * np.pi * r)[:, None] * rho[:, col:col + 2 * l + 1])
        q = np.sum((W * r ** (l + 2))[:, None] * rho[:, col:col + 2 * l + 1], axis=0)
        Urmax = 4 * np.pi * q / ((2 * l + 1) * rmax ** l)
        U = Ainv_rhs / np.sqrt(W)[:, None] + Urmax[None, :] * (r[:, None] / rmax) ** (l + 1)
        V[:, col:col + 2 * l + 1] = U / r[:, None]
        col += 2 * l + 1
    return V


def grid_for(distances, order=12):
    pts = sorted(set([0, 0.04, 0.12, 0.3, 0.6] + [round(d, 6) for d in distances]))
    top = max(distances)
    pts += [top + x for x in (0.5, 1.2, 2.2, 3.6, 5.5, 8, 11.5, 16, 22, 30, 41, 55)]
    return fedvr_w(sorted(set(pts)), order)


def eri_pairwise(a, b, c, d, L, order=12, nw=48):
    """(ab|cd) with a, c on the origin; returns partial sums by L."""
    dist = [float(np.linalg.norm(b['center'])), float(np.linalg.norm(d['center']))]
    r, T, W, rmax = grid_for(dist, order)
    rho1 = channels(a, b, r, L, nw)
    rho2 = channels(c, d, r, L, nw)
    V1 = coulomb_channels(rho1, r, T, W, rmax, L)
    per = np.sum((W * r ** 2)[:, None] * rho2 * V1, axis=0)
    partial, col = [], 0
    for l in range(L + 1):
        partial.append(float(np.sum(per[:col + 2 * l + 1])))
        col += 2 * l + 1
    return partial


# ---------------------------------------------------------------- Gaussian referee
def eri_gauss(a, b, c, d, n_gauss):
    from geovac import noci_engine as GE
    lmn = {'s': (0, 0, 0), 'px': (1, 0, 0), 'py': (0, 1, 0), 'pz': (0, 0, 1)}
    fns = []
    for o in (a, b, c, d):
        l = 0 if o['kind'] == 's' else 1
        al, dc, _ = GE.fit_sto_shape(l, o['n'], n_gauss=n_gauss)
        fns.append(GE.sto_shape_basis(o['center'], 'x', o['zeta'], {'x': (al, dc)}, lmn[o['kind']]))
    nrm = [math.sqrt(GE.overlap_md(f, f)) for f in fns]
    return GE.eri_md(*fns) / (nrm[0] * nrm[1] * nrm[2] * nrm[3])


# ================================================================ tests
# Converged water values: L = 20, order 12, nw = 48 (debug/pairwise_general.py,
# 2026-10-10), each confirmed by the Gaussian-fit sequence (n_gauss 8/12/16,
# agreement ~1e-9; the O 1s core case needs n_gauss -> 18, approaching from below).
WATER_PINS = {'(O2px H1|O2pz H2)': 0.027937789713, '(H1 O2px|H1 H2)': 0.052798542906}


def _water(name):
    O_, H1_, H2_ = O, H1, H2
    o = lambda kind: orbital(O_, 2, 2.275, kind)
    h1, h2 = orbital(H1_, 1, 1.0, 's'), orbital(H2_, 1, 1.0, 's')
    if name == '(O2px H1|O2pz H2)':
        return o('px'), h1, o('pz'), h2
    sh = lambda orb: dict(orb, center=orb['center'] - H1_)
    return tuple(sh(x) for x in (h1, o('px'), h1, h2))


def test_position_route_matches_momentum_route_on_1s():
    """Two independent constructions of the pairwise sum (position-space
    multipoles here; momentum-space channels in test_pairwise_three_center)."""
    spec = importlib.util.spec_from_file_location('pw', os.path.join(ROOT, 'tests', 'test_pairwise_three_center.py'))
    pw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pw)
    mom = pw.channels(20)[0][-1]['partial_sum']
    a = orbital(pw.X, 1, 1.0, 's')
    pos = eri_pairwise(a, orbital(pw.Y, 1, 1.0, 's'), a, orbital(pw.Z, 1, 1.0, 's'), 12, 8, 32)[-1]
    assert abs(pos - mom) < 5e-9


def test_water_three_centre_values():
    """Shared-O (p-p) and shared-H cases at water's geometry, cheap settings."""
    for name, pin in WATER_PINS.items():
        val = eri_pairwise(*_water(name), 12, 8, 32)[-1]
        assert abs(val - pin) < 1e-8, (name, val)


def test_channel_convergence_is_fast():
    part = eri_pairwise(*_water('(O2px H1|O2pz H2)'), 12, 8, 32)
    assert abs(part[0] - part[-1]) > 1e-3          # the l = 0 term alone is not the answer
    assert abs(part[8] - part[-1]) < 2e-7


@pytest.mark.slow
def test_water_value_matches_gaussian_referee():
    a, b, c, d = _water('(O2px H1|O2pz H2)')
    assert abs(eri_gauss(a, b, c, d, 12) - WATER_PINS['(O2px H1|O2pz H2)']) < 1e-8
