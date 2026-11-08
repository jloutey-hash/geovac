"""Decisive test: is the four-vertex invariant Tr(P1P2P3P4) GENUINELY four-body, i.e. NOT a
function of the pairwise {Tr(PiPj)} and triple {Tr(PiPjPk)} data?  (2026-08-25.)

Local test (rigorous): parametrize four rank-r subspaces by their frames; compute the Jacobian
of the feature vector f = [6 pairwise traces, 4 triple traces] and of the target
t = Tr(P1P2P3P4).  Project grad(t) onto the row space of J_f.  The residual
|| grad t - proj || / || grad t ||  is 0 iff t is locally a function of f (pairwise+triple
determines it).  Residual > 0  =>  a genuine four-body degree of freedom.

Prediction: rank 1 -> residual 0 (four real lines are fully pairwise); rank >= 2 -> residual > 0
(genuine four-body content; the non-abelian loop holonomy of the overlap connection).  This is
the geometry-axis analog of the three-center nested commutator, and the seat of a continuous
four-body modulus (the D~4 cross-ratio), which the three-center Z2 CIs did not have.
"""
from __future__ import annotations
import numpy as np


def _proj(Q):
    """projector onto column span of Q (d x r), standard inner product."""
    return Q @ np.linalg.pinv(Q)


def _features_target(theta, d, r):
    """theta: flat (4*d*r,) frames -> (features[10], target) via projectors."""
    Qs = [theta[k * d * r:(k + 1) * d * r].reshape(d, r) for k in range(4)]
    Ps = [_proj(Q) for Q in Qs]
    tr = lambda *ix: float(np.trace(np.linalg.multi_dot([Ps[i] for i in ix])))
    pairs = [tr(i, j) for i in range(4) for j in range(i + 1, 4)]          # 6
    triples = [tr(i, j, k) for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4)]  # 4
    feats = np.array(pairs + triples)
    target = tr(0, 1, 2, 3)
    return feats, target


def residual(d, r, seed):
    rng = np.random.default_rng(seed)
    theta = rng.standard_normal(4 * d * r)
    n = theta.size
    h = 1e-6
    f0, t0 = _features_target(theta, d, r)
    Jf = np.zeros((f0.size, n)); Jt = np.zeros(n)
    for p in range(n):
        tp = theta.copy(); tp[p] += h
        tm = theta.copy(); tm[p] -= h
        fp, tpv = _features_target(tp, d, r)
        fm, tmv = _features_target(tm, d, r)
        Jf[:, p] = (fp - fm) / (2 * h)
        Jt[p] = (tpv - tmv) / (2 * h)
    # project grad(t) onto row space of Jf: proj = Jf^T (Jf Jf^T)^+ Jf grad t
    P_rows = Jf.T @ np.linalg.pinv(Jf @ Jf.T) @ Jf
    proj = P_rows @ Jt
    res = np.linalg.norm(Jt - proj) / max(np.linalg.norm(Jt), 1e-30)
    return res, np.linalg.norm(Jt)


def holonomy_gauge_invariance(d, r, seed):
    """Tr(S12 S23 S34 S41)-type loop holonomy is gauge-invariant under per-center O(r)."""
    rng = np.random.default_rng(seed)
    Qs = [np.linalg.qr(rng.standard_normal((d, r)))[0] for _ in range(4)]
    S = lambda i, j: Qs[i].T @ Qs[j]
    holo = lambda Q: np.trace(Q[0].T @ Q[1] @ Q[1].T @ Q[2] @ Q[2].T @ Q[3] @ Q[3].T @ Q[0])
    h0 = holo(Qs)
    Qg = [Qs[k] @ np.linalg.qr(rng.standard_normal((r, r)))[0] for k in range(4)]   # per-center gauge
    h1 = holo(Qg)
    return abs(h0 - h1)


if __name__ == "__main__":
    print("Genuine-four-body residual  || grad t - proj_f(grad t) || / || grad t ||")
    print("(0 => Tr(P1P2P3P4) is a function of pairwise+triple; >0 => genuine four-body dof)\n")
    for r in [1, 2, 3]:
        d = 8
        rr = [residual(d, r, s) for s in range(4)]
        res = np.mean([x[0] for x in rr]); gt = np.mean([x[1] for x in rr])
        verdict = "PAIRWISE (no four-body dof)" if res < 1e-6 else "GENUINE FOUR-BODY"
        print("  rank r=%d (dim %d): residual = %.2e   (|grad t|=%.2f)  -> %s"
              % (r, d, res, gt, verdict))
    print("\nLoop holonomy gauge-invariance |h(gauged)-h| (should be ~0):")
    for r in [1, 2, 3]:
        print("  rank r=%d: %.2e" % (r, np.mean([holonomy_gauge_invariance(8, r, s) for s in range(4)])))
