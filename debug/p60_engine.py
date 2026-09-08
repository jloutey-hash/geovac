"""Grid harness for the Paper-60 secular matrix.

Two numerical routes over the SAME physics code (geovac.sturmian_secular):

  route 'uni'   -- the production uniform mesh, r = linspace(1e-7, R_MAX, N)
  route 'grade' -- a power-graded mesh r = R_MAX * t^p, t uniform in (0,1]

The graded route requires the two cumulative integrators to become
non-uniform-safe; those are the ONLY two places the module assumes a constant
dr (``cumulative_trapezoid(y, dx=dr)``).  Everything else already integrates
against ``r`` explicitly, so patching them is sufficient and is verified by
reproducing the uniform route at matched resolution.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

import geovac.sturmian_secular as S

_ORIG_FWD = S._ctrap_fwd
_ORIG_REV = S._ctrap_rev


def _fwd_nonuniform(y: np.ndarray) -> np.ndarray:
    return np.concatenate(([0.0], cumulative_trapezoid(y, x=S.r)))


def _rev_nonuniform(y: np.ndarray) -> np.ndarray:
    F = _fwd_nonuniform(y)
    return F[-1] - F


def set_grid(rmax: float, npts: int, kind: str = "uni", p: float = 2.0) -> None:
    """Install a radial mesh.  r, dr, r2 are patched together (they must be)."""
    if kind == "uni":
        S._ctrap_fwd, S._ctrap_rev = _ORIG_FWD, _ORIG_REV
        S.r = np.linspace(1e-7, float(rmax), int(npts))
    elif kind == "grade":
        S._ctrap_fwd, S._ctrap_rev = _fwd_nonuniform, _rev_nonuniform
        t = np.linspace(1e-4, 1.0, int(npts))
        S.r = float(rmax) * t ** p
    else:
        raise ValueError(kind)
    S.R_MAX = float(rmax)
    S.N_GRID = int(npts)
    S.dr = S.r[1] - S.r[0]
    S.r2 = S.r * S.r
    S._GAUNT_CACHE.clear()
    S.reset_caches()


def norms(cts, Z: float = 2.0) -> dict:
    """Every 1-norm leg of M = diag(Z R_nu) + T'."""
    cfgs = S.build_configs(cts)
    Tp = S.build_Tprime(cfgs)
    Rnu = np.array([c.Rnu for c in cfgs])
    M = Tp.copy()
    M[np.diag_indices_from(M)] += Z * Rnu
    dTp = float(np.abs(np.diag(Tp)).sum())
    tot = float(np.abs(Tp).sum())
    return dict(
        K=len(cfgs),
        M_total=float(np.abs(M).sum()),
        M_diag=float(np.abs(np.diag(M)).sum()),
        M_off=tot - dTp,
        Tp_full=tot,
        Tp_diag=dTp,
        T0=float(Z * Rnu.sum()),
        E=float(-np.sort(np.linalg.eigvalsh(M))[-1] ** 2 / 2),
    )


def family(nmax: int, lmax: int = 3):
    return S.gen_configs(lmax, {l: nmax for l in range(lmax + 1)})
