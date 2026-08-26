"""Backing tests for the Paper 32 remark rem:config_operator_is_metric: the configuration
operator F = sum_i P_i is not an abstract object -- it is SIMILAR to the composed basis's
overlap (Gram) matrix G, whenever the intra-center blocks of G are the identity.

Claims pinned:
  (1) IDENTITY: F = X G X^{-1} with X = chol(G)^T, hence spec(F) = spec(G).  Pinned both
      abstractly (random admissible G) and on validated molecular overlaps (BeH2 at the
      conical intersection d* = 2.445 and away from it).
  (2) NON-TAUTOLOGY: the intra-center orthonormality hypothesis is load-bearing -- perturbing
      one intra-center block off the identity breaks the identity by ~1e-2, six orders of
      magnitude above the 1e-13 tolerance the true case satisfies.
  (3) FOUR CENTERS, blockwise: Tr(P1P2P3P4) = Tr(G12 G23 G34 G41) -- the loop holonomy of
      rem:four_center_modulus is the same statement read on the metric's off-diagonal blocks.
  (4) ORTHOGONALIZATION CONSEQUENCE: around a conical intersection the METRIC's eigenvectors
      are double-valued (transport return -1) while a CI-free control returns +1, whereas
      symmetric Lowdin G^{-1/2} -- a function of the matrix, not of its eigenvectors --
      returns bit-exactly on BOTH loops.
  (5) HONEST NEGATIVE: the intersection is a MID-spectrum degeneracy and cond(G) is monotone
      through d*; this is deliberately NOT the Paper 60 metric-conditioning wall.
"""
import numpy as np
import pytest
from scipy.linalg import fractional_matrix_power
from scipy.special import eval_genlaguerre, lpmv, factorial
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss

# --------------------------------------------------------------- overlap engine (validated)
STATES = [(1, 0), (2, 1)]                        # sigma pair 1s, 2p0 (m=0)
PAR = np.diag([1.0, -1.0])                       # 2p0 odd under z->-z
_LAG_X, _LAG_W = laggauss(64)
_LEG_X, _LEG_W = leggauss(96)
DSTAR = 2.445                                    # central conical intersection


def _R_nl(Z, n, l, r):
    Z = float(Z); rho = 2 * Z * r / n
    norm = np.sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * eval_genlaguerre(n - l - 1, 2 * l + 1, rho)


def _ang(l, m, ct):
    m = abs(m)
    tn = np.sqrt((2 * l + 1) / 2.0 * factorial(l - m) / factorial(l + m))
    return tn * lpmv(m, l, ct)


def _overlap(Z1, n1, l1, Z2, n2, l2, R):
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
    return np.array([[_overlap(Za, na, la, Zb, nb, lb, R)
                      for (nb, lb) in STATES] for (na, la) in STATES])


def gram_beh2(d1, d2):
    """Union-basis Gram of linear BeH2 (Be at 0, H1 at +d1, H2 at -d2); intra blocks = I."""
    S1 = _Sblock(d1, 2, 1)
    S2 = PAR @ _Sblock(d2, 2, 1) @ PAR
    SH = PAR @ _Sblock(d1 + d2, 1, 1) @ PAR
    I = np.eye(2)
    return np.block([[I, S1, S2], [S1.T, I, SH], [S2.T, SH.T, I]])


# --------------------------------------------------------------- the construction under test
def F_and_X(G, block_sizes):
    """F = sum_k P_k built by the drivers' own recipe: Cholesky-whiten, then project."""
    X = np.linalg.cholesky(G).T
    F = np.zeros_like(G)
    o = 0
    for b in block_sizes:
        W = X[:, o:o + b]
        F += W @ np.linalg.pinv(W)          # pinv, NOT W W^T: collapses only if W^T W = I
        o += b
    return F, X


def similarity_residuals(G, block_sizes):
    F, X = F_and_X(G, block_sizes)
    dspec = float(np.abs(np.linalg.eigvalsh(F) - np.linalg.eigvalsh(G)).max())
    dsim = float(np.abs(F @ X - X @ G).max())
    return dspec, dsim


# --------------------------------------------------------------- (1) the identity
@pytest.mark.parametrize("d1,d2", [(2.0, 2.0), (2.3, 2.3), (DSTAR, DSTAR),
                                   (2.5, 2.5), (3.0, 3.0), (2.698, 2.266), (6.0, 6.0)])
def test_F_is_similar_to_the_metric_molecular(d1, d2):
    """F ~ G on validated BeH2 overlaps, including AT the conical intersection."""
    G = gram_beh2(d1, d2)
    assert np.linalg.eigvalsh(G).min() > 1e-9, "geometry must not be overcomplete"
    dspec, dsim = similarity_residuals(G, [2, 2, 2])
    assert dspec < 1e-13, f"spec(F) != spec(G): {dspec:.2e}"
    assert dsim < 1e-13, f"F X != X G: {dsim:.2e}"


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_F_is_similar_to_the_metric_abstract(seed):
    """The identity is linear algebra: it holds for ANY admissible G, not just molecular ones."""
    rng = np.random.default_rng(seed)
    n_blocks, b = 3, 2
    A = rng.standard_normal((n_blocks * b, n_blocks * b))
    G = A.T @ A + 3.0 * np.eye(n_blocks * b)
    D = np.diag(1.0 / np.sqrt(np.diag(G)))
    G = D @ G @ D                                  # unit diagonal
    for k in range(n_blocks):                      # force intra blocks to identity
        G[k * b:(k + 1) * b, k * b:(k + 1) * b] = np.eye(b)
    if np.linalg.eigvalsh(G).min() <= 1e-9:
        pytest.skip("random draw not positive definite")
    dspec, dsim = similarity_residuals(G, [b] * n_blocks)
    assert dspec < 1e-12 and dsim < 1e-12


# --------------------------------------------------------------- (2) non-tautology control
def test_intra_orthonormality_is_load_bearing():
    """Break W_k^T W_k = I and the identity must FAIL, by orders of magnitude."""
    G = gram_beh2(2.6, 2.6)
    good, _ = similarity_residuals(G, [2, 2, 2])

    Gb = G.copy()
    Gb[0, 1] = Gb[1, 0] = 0.3                      # center 0's two functions now non-orthogonal
    assert np.linalg.eigvalsh(Gb).min() > 1e-9
    bad, _ = similarity_residuals(Gb, [2, 2, 2])

    assert good < 1e-13, f"true case should hold: {good:.2e}"
    assert bad > 1e-3, f"perturbed case should break: {bad:.2e}"
    assert bad / max(good, 1e-16) > 1e6, "the hypothesis must be discriminating"


# --------------------------------------------------------------- (3) four-center, blockwise
def test_four_center_holonomy_is_a_metric_block_product():
    """Tr(P1P2P3P4) = Tr(G12 G23 G34 G41): the loop holonomy read off the metric blocks."""
    rng = np.random.default_rng(7)
    nb, b = 4, 2
    A = rng.standard_normal((nb * b, nb * b))
    G = A.T @ A + 4.0 * np.eye(nb * b)
    D = np.diag(1.0 / np.sqrt(np.diag(G)))
    G = D @ G @ D
    for k in range(nb):
        G[k * b:(k + 1) * b, k * b:(k + 1) * b] = np.eye(b)
    assert np.linalg.eigvalsh(G).min() > 1e-9

    X = np.linalg.cholesky(G).T
    Ps = [X[:, k * b:(k + 1) * b] @ np.linalg.pinv(X[:, k * b:(k + 1) * b]) for k in range(nb)]
    w_proj = float(np.real(np.trace(Ps[0] @ Ps[1] @ Ps[2] @ Ps[3])))

    def blk(i, j):
        return G[i * b:(i + 1) * b, j * b:(j + 1) * b]
    w_gram = float(np.trace(blk(0, 1) @ blk(1, 2) @ blk(2, 3) @ blk(3, 0)))

    assert abs(w_proj - w_gram) < 1e-12, f"{w_proj} vs {w_gram}"
    assert abs(w_proj) > 1e-3, "guard: a null holonomy would make the test vacuous"


# --------------------------------------------------------------- (4) orthogonalization
def _loop(cx, cy, r, N=120):
    """Transport the near-degenerate METRIC eigenvector once around a loop.

    Returns (eigenvector return overlap, Lowdin RETURN drift).  Single-valuedness is a
    return property: the Lowdin drift compares G^{-1/2} at loop end against loop start.
    """
    v0 = vp = X0 = X_end = None
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        G = gram_beh2(cx + r * np.cos(th), cy + r * np.sin(th))
        w, V = np.linalg.eigh(G)
        X = fractional_matrix_power(G, -0.5).real
        if i == 0:
            k = int(np.argmin(np.diff(w)))
            vp = v0 = V[:, k].copy()
            X0 = X
        else:
            ov = V.T @ vp
            j = int(np.argmax(np.abs(ov)))
            vp = V[:, j] * np.sign(ov[j])
        X_end = X
    return float(vp @ v0), float(np.abs(X_end - X0).max())


def test_metric_eigenvectors_double_valued_lowdin_single_valued():
    sgn_ci, drift_ci = _loop(DSTAR, DSTAR, 0.05)
    sgn_ok, drift_ok = _loop(1.60, 1.60, 0.05)

    assert sgn_ci < -0.99, f"expected -1 around the CI, got {sgn_ci:+.4f}"
    assert sgn_ok > +0.99, f"expected +1 on the CI-free control, got {sgn_ok:+.4f}"
    # Lowdin is a matrix function => single-valued on BOTH loops.
    assert drift_ci < 1e-12 and drift_ok < 1e-12


# --------------------------------------------------------------- (5) honest negative
def test_ci_is_midspectrum_and_cond_is_monotone_through_dstar():
    """Deliberately pins the NEGATIVE: this is not the Paper 60 conditioning wall."""
    w = np.linalg.eigvalsh(gram_beh2(DSTAR, DSTAR))
    gaps = np.diff(w)
    k = int(np.argmin(gaps))
    assert gaps[k] < 1e-4, "the CI must actually be near-degenerate here"
    # crossing sits mid-spectrum, well above the smallest eigenvalue
    assert w[k] > 3 * w[0], f"crossing {w[k]:.3f} should be far above lambda_min {w[0]:.3f}"
    assert 0 < k < len(w) - 1, "crossing must be interior to the spectrum"

    conds = [np.linalg.cond(gram_beh2(d, d)) for d in (2.3, DSTAR, 2.5)]
    assert conds[0] > conds[1] > conds[2], f"cond should be monotone through d*: {conds}"
