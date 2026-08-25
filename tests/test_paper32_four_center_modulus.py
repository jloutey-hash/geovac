"""Backing tests for the Paper 32 remark rem:four_center_modulus: a FOUR-center configuration
operator F = sum_{i=1}^4 P_i carries a genuine CONTINUOUS four-body modulus -- the gauge-invariant
non-abelian loop holonomy of the overlap connection around the 4-cycle -- which the three-center
case (discrete Z2 conical intersections) did not have.  It requires rank >= 2 (four real rank-1
subspaces are fully pairwise) and complexifies to a U(1) holonomy under a magnetic flux.

Claims pinned:
  (1) THRESHOLD: the four-vertex invariant Tr(P1P2P3P4) is a FUNCTION of pairwise+triple data for
      four rank-1 subspaces (residual ~0) but a GENUINE four-body dof for rank>=2 (residual >0.1)
      -- local Jacobian nullspace test;
  (2) the loop holonomy is GAUGE-INVARIANT under per-center O(r) (a real physical invariant);
  (3) molecular realization: over a rhombus shape sweep the holonomy is a smooth, monotone,
      continuous four-body modulus (embedded validated 1s/2p Slater-Koster overlaps);
  (4) U(1) lift: a small magnetic flux on the 4-cycle makes the holonomy phase grow continuously
      from 0 (Z2 at zero flux -> continuous U(1)).
"""
import numpy as np
import pytest
from scipy.special import eval_genlaguerre, lpmv, factorial
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss

# --------------------------- abstract subspace machinery (claims 1,2) ---------------------------
def _proj(Q):
    return Q @ np.linalg.pinv(Q)


def _feat_target(theta, d, r):
    Qs = [theta[k*d*r:(k+1)*d*r].reshape(d, r) for k in range(4)]
    Ps = [_proj(Q) for Q in Qs]
    tr = lambda *ix: float(np.trace(np.linalg.multi_dot([Ps[i] for i in ix])))
    feats = np.array([tr(i, j) for i in range(4) for j in range(i+1, 4)] +
                     [tr(i, j, k) for i in range(4) for j in range(i+1, 4) for k in range(j+1, 4)])
    return feats, tr(0, 1, 2, 3)


def _fourbody_residual(d, r, seed):
    rng = np.random.default_rng(seed)
    theta = rng.standard_normal(4*d*r); n = theta.size; h = 1e-6
    Jf = np.zeros((10, n)); Jt = np.zeros(n)
    for p in range(n):
        tp = theta.copy(); tp[p] += h; tm = theta.copy(); tm[p] -= h
        fp, tpv = _feat_target(tp, d, r); fm, tmv = _feat_target(tm, d, r)
        Jf[:, p] = (fp - fm)/(2*h); Jt[p] = (tpv - tmv)/(2*h)
    proj = Jf.T @ np.linalg.pinv(Jf @ Jf.T) @ Jf @ Jt
    return np.linalg.norm(Jt - proj)/max(np.linalg.norm(Jt), 1e-30)


def test_fourbody_threshold_rank1_pairwise_rank2_genuine():
    """(1) rank-1 four real subspaces: Tr(P1P2P3P4) is pairwise+triple-determined (residual ~0);
    rank-2: a genuine four-body dof (residual > 0.1)."""
    r1 = np.mean([_fourbody_residual(6, 1, s) for s in range(3)])
    r2 = np.mean([_fourbody_residual(6, 2, s) for s in range(3)])
    assert r1 < 1e-6            # rank-1: no four-body dof (fully pairwise)
    assert r2 > 0.1            # rank-2: genuine four-body degree of freedom


def test_loop_holonomy_gauge_invariant():
    """(2) Tr(P1P2P3P4)-type loop holonomy is invariant under per-center O(r) gauge."""
    rng = np.random.default_rng(0)
    Qs = [np.linalg.qr(rng.standard_normal((8, 2)))[0] for _ in range(4)]
    holo = lambda Q: np.trace(Q[0].T@Q[1]@Q[1].T@Q[2]@Q[2].T@Q[3]@Q[3].T@Q[0])
    Qg = [Qs[k] @ np.linalg.qr(rng.standard_normal((2, 2)))[0] for k in range(4)]
    assert abs(holo(Qs) - holo(Qg)) < 1e-12


# --------------------------- embedded overlap engine (claims 3,4) ---------------------------
STATES = [(1, 0), (2, 1)]
_LAG_X, _LAG_W = laggauss(64); _LEG_X, _LEG_W = leggauss(96)


def _R_nl(Z, n, l, r):
    Z = float(Z); rho = 2*Z*r/n
    norm = np.sqrt((2*Z/n)**3 * factorial(n-l-1)/(2*n*factorial(n+l)))
    return norm*np.exp(-rho/2)*rho**l*eval_genlaguerre(n-l-1, 2*l+1, rho)


def _ang(l, m, ct):
    m = abs(m); tn = np.sqrt((2*l+1)/2.0*factorial(l-m)/factorial(l+m))
    return tn*lpmv(m, l, ct)


def _ov(n1, l1, n2, l2, R):
    half = R/2.0; a = half*(1/n1 + 1/n2); xi = 1.0 + _LAG_X/a
    XI, ETA = np.meshgrid(xi, _LEG_X, indexing="ij")
    r1 = half*(XI+ETA); r2 = half*(XI-ETA)
    with np.errstate(divide="ignore", invalid="ignore"):
        ct1 = (1+XI*ETA)/(XI+ETA); ct2 = (XI*ETA-1)/(XI-ETA)
    val = _R_nl(1, n1, l1, r1)*_ang(l1, 0, ct1)*_R_nl(1, n2, l2, r2)*_ang(l2, 0, ct2)*half**3*(XI**2-ETA**2)
    val = np.where((r1 <= 0) | (r2 <= 0) | ~np.isfinite(val), 0.0, val)
    return float(np.sum((_LAG_W*np.exp(_LAG_X)/a)[:, None]*(_LEG_W[None, :]*val)))


def _sk(R, u):
    ss, SP, PS = _ov(1, 0, 1, 0, R), _ov(1, 0, 2, 1, R), _ov(2, 1, 1, 0, R)
    pps, ppp = _ov(2, 1, 2, 1, R), _ov(2, 1, 2, 1, R)  # note ppp needs m=1; recompute:
    ppp = _ov_m1(R)
    Mloc = np.array([[ss, 0.0, SP], [0.0, ppp, 0.0], [PS, 0.0, pps]])
    ux, uz = u; C = np.array([[1.0, 0, 0], [0, uz, ux], [0, -ux, uz]])
    return C @ Mloc @ C.T


def _ov_m1(R):
    half = R/2.0; a = half*(1/2 + 1/2); xi = 1.0 + _LAG_X/a
    XI, ETA = np.meshgrid(xi, _LEG_X, indexing="ij")
    r1 = half*(XI+ETA); r2 = half*(XI-ETA)
    with np.errstate(divide="ignore", invalid="ignore"):
        ct1 = (1+XI*ETA)/(XI+ETA); ct2 = (XI*ETA-1)/(XI-ETA)
    val = _R_nl(1, 2, 1, r1)*_ang(1, 1, ct1)*_R_nl(1, 2, 1, r2)*_ang(1, 1, ct2)*half**3*(XI**2-ETA**2)
    val = np.where((r1 <= 0) | (r2 <= 0) | ~np.isfinite(val), 0.0, val)
    return float(np.sum((_LAG_W*np.exp(_LAG_X)/a)[:, None]*(_LEG_W[None, :]*val)))


def _G4(positions):
    G = np.zeros((12, 12))
    for i in range(4):
        G[3*i:3*i+3, 3*i:3*i+3] = np.eye(3)
        for j in range(i+1, 4):
            d = np.array(positions[j]) - np.array(positions[i]); R = np.hypot(*d)
            S = _sk(R, (d[0]/R, d[1]/R)); G[3*i:3*i+3, 3*j:3*j+3] = S; G[3*j:3*j+3, 3*i:3*i+3] = S.T
    return G


def _holo(G):
    Xh = np.linalg.cholesky(G).conj().T if np.iscomplexobj(G) else np.linalg.cholesky(G).T
    P = [Xh[:, 3*k:3*k+3] @ np.linalg.pinv(Xh[:, 3*k:3*k+3]) for k in range(4)]
    return complex(np.trace(P[0] @ P[1] @ P[2] @ P[3]))


def test_molecular_continuous_modulus_and_u1_flux():
    """(3) smooth monotone continuous modulus over a rhombus sweep; (4) U(1): flux -> continuous arg."""
    R = 2.5
    ts = [0.6, 0.9, 1.2, 1.5]
    hol = [_holo(_G4([(-R, 0), (0, -R*t), (R, 0), (0, R*t)])).real for t in ts]
    diffs = np.diff(hol)
    assert all(d < 0 for d in diffs)                       # strictly decreasing (continuous modulus)
    assert max(abs(np.diff(diffs))) < abs(diffs).max()     # smooth (2nd diff < 1st diff scale)

    sq = [(-R/2, -R/2), (R/2, -R/2), (R/2, R/2), (-R/2, R/2)]
    args = []
    for eps in [0.0, 0.1, 0.2]:
        G = _G4(sq).astype(complex)
        for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 0)]:
            G[3*i, 3*j] *= np.exp(1j*eps); G[3*j, 3*i] *= np.exp(-1j*eps)
        args.append(np.angle(_holo(G)))
    assert abs(args[0]) < 1e-9                              # Z2 (real) at zero flux
    assert args[1] > 1e-3 and args[2] > args[1]            # continuous U(1) phase for eps != 0
