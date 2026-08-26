"""Probe A -- LP core: minimal QSVT degree for x^{-1/2} under a boundedness constraint.

The honest constrained object (see debug/probeA_findings.md):

    find p, deg p <= d,   |p(x)| <= 1  for all x in [-1, 1]   (QSP realizability)
    with  |p(x_i) - c * lam_i^{-1/2}| <= eps * c * lam_i^{-1/2}   at the accuracy set.

We solve, per degree d, the *linear program*

    minimise r   s.t.  |p(x_i) - y_i| <= r * y_i  (accuracy set),  |p(t_g)| <= 1 (grid),

whose optimum r*(d) is the best achievable RELATIVE error at degree d; the minimal
degree at tolerance eps is then  d(eps) = min{ d : r*(d) <= eps }, found by bisection.
r*(d) is monotone non-increasing in d, so bisection is valid.

Two accuracy sets are compared with identical machinery:
  * "aware"   -- only the 2n known eigenvalues (the exploit under test);
  * "generic" -- a dense grid on the whole continuum [lam_min, lam_max].

Two rescaling conventions:
  * "p60"   -- x = lam / lam_max in [1/kappa, 1], |p| <= 1 on [-1,1]
               (exactly the convention Paper 60 sec:resource prices, d ~ kappa ln(kappa/eps));
  * "shift" -- x = affine([lam_min, lam_max] -> [-1,1]) (needs lam_min, i.e. the sigma law).

Chebyshev basis: a polynomial with |p| <= 1 on [-1,1] has |a_j| <= 2, so the LP is
well scaled; we impose -2.05 <= a_j <= 2.05 as harmless variable bounds.

Boundedness is enforced on m = BND_PER_DEG * (d+1) Chebyshev extreme points.  By the
Ehlich-Zeller bound the true supremum of a degree-d polynomial bounded by 1 on m such
points is at most 1/cos(pi*d/(2m)); with BND_PER_DEG = 12 that is 1.0086.  Every reported
solution is re-checked a posteriori on a 10x finer grid.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

BND_PER_DEG = 12
COEF_BOUND = 2.05


def cheb_design(x: np.ndarray, d: int) -> np.ndarray:
    """T_j(x) for j = 0..d, evaluated stably as cos(j arccos x)."""
    x = np.clip(np.asarray(x, dtype=float), -1.0, 1.0)
    th = np.arccos(x)
    j = np.arange(d + 1)
    return np.cos(np.outer(th, j))


def bound_grid(d: int, extra: np.ndarray | None = None) -> np.ndarray:
    m = BND_PER_DEG * (d + 1) + 1
    g = np.cos(np.pi * np.arange(m) / (m - 1))
    if extra is not None:
        g = np.concatenate([g, np.clip(extra, -1.0, 1.0)])
    return np.unique(g)


def min_rel_error(nodes_x: np.ndarray, nodes_y: np.ndarray, d: int,
                  bound: float = 1.0) -> tuple[float, np.ndarray]:
    """min over deg<=d of max_i |p(x_i)-y_i|/y_i  s.t.  |p| <= bound on [-1,1]."""
    A_acc = cheb_design(nodes_x, d)
    g = bound_grid(d, extra=nodes_x)
    A_bnd = cheb_design(g, d)
    nv = d + 2                       # coefficients + r
    n_acc, n_bnd = A_acc.shape[0], A_bnd.shape[0]
    A_ub = np.zeros((2 * n_acc + 2 * n_bnd, nv))
    b_ub = np.zeros(2 * n_acc + 2 * n_bnd)
    # p(x_i) - y_i <= r y_i
    A_ub[:n_acc, :d + 1] = A_acc
    A_ub[:n_acc, -1] = -nodes_y
    b_ub[:n_acc] = nodes_y
    # -(p(x_i) - y_i) <= r y_i
    A_ub[n_acc:2 * n_acc, :d + 1] = -A_acc
    A_ub[n_acc:2 * n_acc, -1] = -nodes_y
    b_ub[n_acc:2 * n_acc] = -nodes_y
    # |p| <= bound
    A_ub[2 * n_acc:2 * n_acc + n_bnd, :d + 1] = A_bnd
    b_ub[2 * n_acc:2 * n_acc + n_bnd] = bound
    A_ub[2 * n_acc + n_bnd:, :d + 1] = -A_bnd
    b_ub[2 * n_acc + n_bnd:] = bound
    c = np.zeros(nv)
    c[-1] = 1.0
    cb = COEF_BOUND * max(1.0, bound)
    bnds = [(-cb, cb)] * (d + 1) + [(0.0, None)]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bnds, method="highs")
    if not res.success:
        return np.inf, np.zeros(d + 1)
    return float(res.x[-1]), np.asarray(res.x[:d + 1])


def true_sup(coef: np.ndarray) -> float:
    """A-posteriori sup |p| on [-1,1] on a 10x finer Chebyshev grid."""
    d = len(coef) - 1
    m = 60 * (d + 1) + 1
    g = np.cos(np.pi * np.arange(m) / (m - 1))
    return float(np.max(np.abs(cheb_design(g, d) @ coef)))


def min_degree(nodes_x: np.ndarray, nodes_y: np.ndarray, eps: float,
               d_lo: int = 1, d_max: int = 600) -> tuple[int, float, np.ndarray]:
    """Smallest d with r*(d) <= eps, by doubling + bisection.  (-1 if > d_max)."""
    d = max(d_lo, 2)
    r, coef = min_rel_error(nodes_x, nodes_y, d)
    if r <= eps:
        hi = d
        lo = 1
    else:
        lo = d
        hi = -1
        while d <= d_max:
            d = min(2 * d, d_max)
            r, coef = min_rel_error(nodes_x, nodes_y, d)
            if r <= eps:
                hi = d
                break
            lo = d
            if d == d_max:
                break
        if hi < 0:
            return -1, r, coef
    best_coef, best_r = coef, r
    while hi - lo > 1:
        mid = (hi + lo) // 2
        r, c = min_rel_error(nodes_x, nodes_y, mid)
        if r <= eps:
            hi, best_coef, best_r = mid, c, r
        else:
            lo = mid
    return hi, best_r, best_coef


# ---------------------------------------------------------------- conventions

def make_problem(lam: np.ndarray, mode: str, convention: str,
                 n_grid: int = 500, headroom: float = 0.5) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (nodes_x in [-1,1], target y, subnormalisation c).

    Normalisation convention: c = headroom * sqrt(lam_min), so the largest required
    value is `headroom`.  The LP's relative-error objective is invariant under
    rescaling of c EXCEPT through the boundedness constraint, so `headroom` is a real
    parameter:  headroom = 1 forces p to touch the QSP ceiling exactly at lam_min
    (a tangency that only converges algebraically), headroom = 1/2 is the standard
    factor-2 subnormalisation head-room of QSVT inversion constructions.  A smaller
    c costs 1/c amplitude-amplification rounds downstream, so it is held FIXED across
    every arm of the comparison and the reported degrees are directly comparable.
    """
    lam = np.sort(np.asarray(lam, dtype=float))
    lo, hi = lam[0], lam[-1]
    if mode == "aware":
        lam_acc = lam
    elif mode == "generic":
        geo = lo * (hi / lo) ** np.linspace(0.0, 1.0, n_grid)
        th = np.pi * np.arange(n_grid) / (n_grid - 1)
        chb = 0.5 * (lo + hi) - 0.5 * (hi - lo) * np.cos(th)
        lam_acc = np.unique(np.concatenate([geo, chb, lam]))
    else:
        raise ValueError(mode)
    if convention == "p60":
        x = lam_acc / hi
    elif convention == "shift":
        x = (2.0 * lam_acc - (lo + hi)) / (hi - lo)
    else:
        raise ValueError(convention)
    c = headroom * np.sqrt(lo)
    y = c / np.sqrt(lam_acc)
    return x, y, c


def map_x(lam: np.ndarray, lam_ref: np.ndarray, convention: str) -> np.ndarray:
    lo, hi = float(np.min(lam_ref)), float(np.max(lam_ref))
    if convention == "p60":
        return np.asarray(lam) / hi
    return (2.0 * np.asarray(lam) - (lo + hi)) / (hi - lo)


# ------------------------------------------------- cutting-plane (fast) solver

def min_rel_error_cp(nodes_x: np.ndarray, nodes_y: np.ndarray, d: int,
                     bound: float = 1.0, n_fine: int = 12, n_start: int = 2,
                     max_rounds: int = 12) -> tuple[float, np.ndarray]:
    """Same LP as min_rel_error, solved by cutting planes on the boundedness set.

    The boundedness constraint is imposed on a growing ACTIVE subset of a fine
    Chebyshev grid (n_fine*(d+1)+1 points); after each solve, every fine-grid point
    where |p| exceeds the bound is added.  On exit the returned p satisfies
    |p| <= bound on the whole fine grid, so (Ehlich-Zeller) sup_{[-1,1]}|p| <=
    bound / cos(pi d / (2 m)) with m = n_fine*(d+1)+1 -- 1.0086 at n_fine=12.
    """
    A_acc = cheb_design(nodes_x, d)
    m = n_fine * (d + 1) + 1
    gf = np.cos(np.pi * np.arange(m) / (m - 1))
    A_fine = cheb_design(gf, d)
    step = max(1, n_fine // n_start)
    active = set(range(0, m, step)) | {0, m - 1}
    active |= set(np.searchsorted(gf[::-1], np.clip(nodes_x, -1, 1)).tolist())
    n_acc = A_acc.shape[0]
    cb = COEF_BOUND * max(1.0, bound)
    bnds = [(-cb, cb)] * (d + 1) + [(0.0, None)]
    cvec = np.zeros(d + 2)
    cvec[-1] = 1.0
    r, coef = np.inf, np.zeros(d + 1)
    for _ in range(max_rounds):
        idx = np.fromiter(sorted(i for i in active if 0 <= i < m), dtype=int)
        A_b = A_fine[idx]
        nb = A_b.shape[0]
        A_ub = np.zeros((2 * n_acc + 2 * nb, d + 2))
        b_ub = np.zeros(2 * n_acc + 2 * nb)
        A_ub[:n_acc, :d + 1] = A_acc;   A_ub[:n_acc, -1] = -nodes_y;  b_ub[:n_acc] = nodes_y
        A_ub[n_acc:2 * n_acc, :d + 1] = -A_acc
        A_ub[n_acc:2 * n_acc, -1] = -nodes_y
        b_ub[n_acc:2 * n_acc] = -nodes_y
        A_ub[2 * n_acc:2 * n_acc + nb, :d + 1] = A_b;   b_ub[2 * n_acc:2 * n_acc + nb] = bound
        A_ub[2 * n_acc + nb:, :d + 1] = -A_b;           b_ub[2 * n_acc + nb:] = bound
        res = linprog(cvec, A_ub=A_ub, b_ub=b_ub, bounds=bnds, method="highs")
        if not res.success:
            return np.inf, np.zeros(d + 1)
        r = float(res.x[-1]); coef = np.asarray(res.x[:d + 1])
        vals = A_fine @ coef
        viol = np.nonzero(np.abs(vals) > bound * (1.0 + 1e-9))[0]
        if viol.size == 0:
            return r, coef
        active |= set(viol.tolist())
    return r, coef


def min_degree_cp(nodes_x: np.ndarray, nodes_y: np.ndarray, eps: float,
                  d_lo: int = 2, d_max: int = 4096,
                  **kw) -> tuple[int, float, np.ndarray]:
    """Smallest d with r*(d) <= eps (cutting-plane LP), doubling + bisection."""
    lo, hi, best = 1, -1, (np.inf, np.zeros(1))
    d = max(2, d_lo)
    while d <= d_max:
        r, coef = min_rel_error_cp(nodes_x, nodes_y, d, **kw)
        if r <= eps:
            hi, best = d, (r, coef)
            break
        lo = d
        if d == d_max:
            return -1, r, coef
        d = min(2 * d, d_max)
    if hi < 0:
        return -1, best[0], best[1]
    while hi - lo > 1:
        mid = (hi + lo) // 2
        r, c = min_rel_error_cp(nodes_x, nodes_y, mid, **kw)
        if r <= eps:
            hi, best = mid, (r, c)
        else:
            lo = mid
    return hi, best[0], best[1]
