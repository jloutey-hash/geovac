"""Variational-CI companion route to the Paper-60 isoenergetic secular equation.

This module is the SECOND, independent numerical route behind Paper 60's
``eq:scale_lock``.  Where :mod:`geovac.sturmian_secular` poses the metric-free
isoenergetic problem (``M = diag(Z R_nu) + T'``, never forming an overlap
metric), this module builds a genuine variational CI over the *same* Goscinskian
configuration span, with the global scale ``lambda`` as the variational
parameter -- the variational analogue of the isoenergetic ``p_kappa``.  The two
routes share only the primitive ``sturmian_secular.repulsion_terms``; route A
never forms ``S`` or ``W``, route B never forms ``M``.  That disjointness is
what makes the scale lock a result rather than a tautology, and it is asserted
operationally in ``tests/test_paper60_scale_lock.py``.

Contents
--------
Grid harness (``set_grid`` and the two non-uniform-safe cumulative integrators
it installs).  Two numerical routes over the SAME physics code:

  route ``'uni'``   -- the production uniform mesh, ``r = linspace(1e-7, R_MAX, N)``
  route ``'grade'`` -- a power-graded mesh ``r = R_MAX * t^p``, ``t`` uniform in (0,1]

The graded route requires the two cumulative integrators of
:mod:`geovac.sturmian_secular` to become non-uniform-safe; those are the ONLY
two places that module assumes a constant ``dr``
(``cumulative_trapezoid(y, dx=dr)``).  Everything else already integrates
against ``r`` explicitly, so patching them is sufficient, and it is verified by
reproducing the uniform route at matched resolution.  The harness is
load-bearing rather than cosmetic: the box-convergence discipline it enforces
(``R_MAX = max(80, 5 n_max^2)`` at ``N = 24000`` graded points) is what caught a
retired Paper-60 exponent.

Variational assembly (``build``) and the whitened generalized eigenproblem
(``var_energy`` / ``var_levels``).  The one-body part uses the hydrogenic
eigen-identity on the ket,

    sum_j (-1/2 grad_j^2) |Phi_nu> = -(Q_nu^2 R_nu^2 / 2)|Phi_nu> + Q_nu U |Phi_nu>,

with ``U = sum_j 1/r_j`` and ``Q_nu = pk_ref / R_nu``.  At ``pk_ref = 1`` the
product ``Q_nu R_nu = 1``, so the first term is exactly ``-1/2``.  The bra/ket
asymmetry of the resulting kinetic matrix is returned as a grid-quality check
(it is zero for the exact integrals and measures radial-grid error only).

Grid mutation
-------------
:func:`set_grid` MUTATES module globals of :mod:`geovac.sturmian_secular`
(``r``, ``dr``, ``r2``, ``R_MAX``, ``N_GRID``, ``_ctrap_fwd``, ``_ctrap_rev``)
and clears that module's caches.  Callers that must leave the grid untouched
for other consumers are responsible for snapshotting and restoring those names
(``tests/test_paper60_scale_lock.py`` does this in a module-scoped autouse
fixture).

Provenance: promoted verbatim (physics unchanged) from the transient drivers
``debug/p60_engine.py`` and ``debug/p60_variational_probe.py``, which are now
thin compatibility shims re-exporting from here, so that the Paper-60 guards do
not depend on the prunable ``debug/`` tree (CLAUDE.md SS9 clean-room rule; gate
C22 check D).
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.integrate import cumulative_trapezoid

import geovac.sturmian_secular as S

# An "orbital" is the dict {rid, l, m, P} of geovac.sturmian_secular; a "term"
# is one (coeff, orb_e1, orb_e2) entry of a Config's term list.
Orbital = Dict[str, Any]
Term = Tuple[float, Orbital, Orbital]

# Pristine (uniform-grid) cumulative integrators, captured before any patching.
_ORIG_FWD = S._ctrap_fwd
_ORIG_REV = S._ctrap_rev


# --------------------------------------------------------------------------------------
# Radial-grid harness.
# --------------------------------------------------------------------------------------
def _fwd_nonuniform(y: np.ndarray) -> np.ndarray:
    """Forward cumulative integral against the (possibly non-uniform) mesh ``S.r``."""
    return np.concatenate(([0.0], cumulative_trapezoid(y, x=S.r)))


def _rev_nonuniform(y: np.ndarray) -> np.ndarray:
    """Reverse cumulative integral against the (possibly non-uniform) mesh ``S.r``."""
    F = _fwd_nonuniform(y)
    return F[-1] - F


def set_grid(rmax: float, npts: int, kind: str = "uni", p: float = 2.0) -> None:
    """Install a radial mesh on :mod:`geovac.sturmian_secular`.

    ``r``, ``dr`` and ``r2`` are patched together (they must be), the matching
    cumulative integrators are installed, and both caches are cleared.

    Parameters
    ----------
    rmax : float
        Box radius ``R_MAX``.
    npts : int
        Number of radial points.
    kind : str
        ``'uni'`` for the production uniform mesh, ``'grade'`` for the
        power-graded mesh ``r = rmax * t^p`` with ``t`` uniform on (0, 1].
    p : float
        Grading exponent, used only for ``kind='grade'``.
    """
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


def family(nmax: int, lmax: int = 3) -> List[Tuple[int, int, int]]:
    """The Paper-60 configuration family: every l <= ``lmax`` at the same ``nmax``."""
    return S.gen_configs(lmax, {l: nmax for l in range(lmax + 1)})


def norms(cts: List[Tuple[int, int, int]], Z: float = 2.0) -> Dict[str, float]:
    """Every 1-norm leg of the isoenergetic secular matrix ``M = diag(Z R_nu) + T'``.

    Diagnostic accounting behind Paper 60 ``eq:sublinear`` (which leg of the
    1-norm carries the sublinear growth); returns the config count ``K``, the
    total / diagonal / off-diagonal 1-norms of ``M`` and ``T'``, the nuclear
    diagonal ``T0 = Z sum R_nu``, and the resulting energy.
    """
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


# --------------------------------------------------------------------------------------
# One-body 1/r primitives (the metric W).
# --------------------------------------------------------------------------------------
def radial_1r(Pa: np.ndarray, Pb: np.ndarray) -> float:
    """``int Pa Pb (1/r) r^2 dr = int Pa Pb r dr`` on the current module grid."""
    return float(np.trapezoid(Pa * Pb * S.r, S.r))


def radial_overlap_c(oa: Orbital, ob: Orbital) -> float:
    """Radial overlap of two orbitals, short-circuited on a shared radial id."""
    if oa['rid'] == ob['rid']:
        return 1.0
    return S.radial_overlap(oa['P'], ob['P'])


def u_terms(termsA: List[Term], termsB: List[Term]) -> float:
    """``<Psi_A | 1/r_1 + 1/r_2 | Psi_B>``, unnormalized."""
    tot = 0.0
    for (wa, ua, va) in termsA:
        for (wb, ub, vb) in termsB:
            if ua['l'] != ub['l'] or ua['m'] != ub['m']:
                continue
            if va['l'] != vb['l'] or va['m'] != vb['m']:
                continue
            su = radial_overlap_c(ua, ub)
            sv = radial_overlap_c(va, vb)
            ru = radial_1r(ua['P'], ub['P'])
            rv = radial_1r(va['P'], vb['P'])
            tot += wa * wb * (ru * sv + su * rv)
    return tot


# --------------------------------------------------------------------------------------
# Variational assembly over the Goscinskian span.
# --------------------------------------------------------------------------------------
def build(nmax: int, lmax: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                         np.ndarray, int, float]:
    """Assemble the variational CI over the ``(nmax, lmax)`` Goscinskian span.

    Returns ``(S_overlap, T, W, G, K, bra_ket_asym)``:

      ``S_overlap``  L2 overlap metric between normalized configurations;
      ``T``          symmetrized one-body kinetic matrix, built from the ket-side
                     hydrogenic eigen-identity
                     ``T[:, nu] = -(Q_nu R_nu)^2 / 2 * S[:, nu] + Q_nu W[:, nu]``;
      ``W``          one-body ``sum_j 1/r_j`` metric (diagonal, ``= R_nu``, by
                     Paper 60 ``eq:W_diagonal``);
      ``G``          interelectron repulsion ``<Psi_i | 1/r12 | Psi_j>``;
      ``K``          configuration count;
      ``bra_ket_asym`` relative bra/ket asymmetry of the unsymmetrized ``T``,
                     a radial-grid quality check (exactly zero in exact arithmetic).
    """
    cfgs = S.build_configs(family(nmax, lmax))
    K = len(cfgs)
    Smat = np.zeros((K, K))
    W = np.zeros((K, K))
    G = np.zeros((K, K))
    for i in range(K):
        ci = cfgs[i]
        for j in range(i, K):
            cj = cfgs[j]
            nn = ci.norm * cj.norm
            Smat[i, j] = Smat[j, i] = nn * S.overlap_terms(ci.terms, cj.terms)
            W[i, j] = W[j, i] = nn * u_terms(ci.terms, cj.terms)
            G[i, j] = G[j, i] = nn * S.repulsion_terms(ci.terms, cj.terms)
    Q = np.array([c.Q for c in cfgs])
    Rn = np.array([c.Rnu for c in cfgs])
    # ket-side kinetic:  T[:,nu] = -(Q_nu^2 R_nu^2/2) S[:,nu] + Q_nu W[:,nu]
    Tk = -0.5 * (Q * Rn) ** 2 * Smat + W * Q
    asym = np.abs(Tk - Tk.T).max() / max(np.abs(Tk).max(), 1e-30)
    return Smat, (Tk + Tk.T) / 2, W, G, K, asym


def _whiten(Smat: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Whitening transform ``X`` with ``X^T S X = I``, dropping the null space.

    Eigenvalues of ``S`` at or below ``tol * max(w)`` are discarded. This is a
    numerical safety net, NOT a physical truncation, and on this family it
    provably never fires: 0 of 78 / 105 / 290 directions were dropped at every
    measured case (cond(S) = 33 / 43 / 56). The whitened problem therefore spans
    the IDENTICAL space, which is what lets `var_energy` be compared with the
    locked posing over 'the same span'.

    An earlier version of this docstring said the metric is 'ill-conditioned by
    construction ... so this truncation is required, not cosmetic'. That rested
    on the retired cond(S) 4 -> 3673 divergence [retracted 2026-09-07: p60-l2-metric-diverges] (a radial-box artifact) and was
    wrong twice over -- the metric is ordinary, and the truncation never fires.
    """
    w, V = np.linalg.eigh(Smat)
    keep = w > tol * w.max()
    return V[:, keep] / np.sqrt(w[keep])


def _hamiltonian(T: np.ndarray, W: np.ndarray, G: np.ndarray,
                 Z: float, lam: float) -> np.ndarray:
    """``H(lam) = lam^2 T + lam (-Z W + G)`` -- the scaled variational Hamiltonian."""
    return lam ** 2 * T + lam * (-Z * W + G)


def var_energy(Smat: np.ndarray, T: np.ndarray, W: np.ndarray, G: np.ndarray,
               Z: float, lam: float, tol: float = 1e-10) -> float:
    """Lowest eigenvalue of ``H(lam) C = E S C`` in the ``S`` metric."""
    H = _hamiltonian(T, W, G, Z, lam)
    X = _whiten(Smat, tol)
    return float(np.linalg.eigvalsh(X.T @ H @ X).min())


def var_levels(Smat: np.ndarray, T: np.ndarray, W: np.ndarray, G: np.ndarray,
               Z: float, lam: float, nlev: int, tol: float = 1e-10) -> np.ndarray:
    """The ``nlev`` lowest roots of ``H(lam) C = E S C`` (ascending)."""
    H = _hamiltonian(T, W, G, Z, lam)
    X = _whiten(Smat, tol)
    return np.linalg.eigvalsh(X.T @ H @ X)[:nlev]
