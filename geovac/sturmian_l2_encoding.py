"""
Atomic L2-Loewdin Jordan-Wigner LCU 1-Norm (Paper 60 eq:blowup)
=================================================================

Faithful block-encoding 1-norm lambda for a naive L2 (Loewdin-orthonormalized,
Jordan-Wigner-encoded) treatment of a single-center, s-only atomic model, built
in two radial families:

    hydrogenic : R_{n0}, decay a = Z/n            (L2-orthonormal across shells)
    sturmian   : S_{n0}, shared decay k = Z fixed  (the genuine Coulomb-Sturmian;
                 L2-NON-orthogonal across shells)

Both bases are Loewdin-orthonormalized through the identical pipeline, their
one- and two-body integrals transformed into the orthonormal basis, encoded as
a restricted (spin-summed) fermionic Hamiltonian, Jordan-Wigner mapped via
openfermion, and the LCU 1-norm lambda = sum_i |c_i| (excluding the identity
term) is measured as a function of shell count N (Q = 2N spin-orbitals).

Paper 60 eq:blowup claims this naive encoding inflates lambda as:
    hydrogenic : lambda ~ Q^1.19
    sturmian   : lambda ~ Q^3.33  (the shared-scale overlap ill-conditions,
                 and Loewdin whitening spreads that ill-conditioning into the
                 transformed integrals -- the motivation for the isoenergetic
                 reformulation, which removes the metric for atoms).

Anchoring k = Z for the Sturmian family makes the n=1 shell identical to the
hydrogenic n=1 shell, so the two sweeps are on a fair, apples-to-apples footing.

Ported from the diagnostic driver debug/io_ladder_sturmian_lambda.py (untracked,
not regression-protected) to close that gap; see CLAUDE.md CHANGELOG entry for
the promotion sprint.
"""

from typing import Dict, List, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray
from openfermion import FermionOperator, jordan_wigner
from scipy.special import genlaguerre

FloatArray = NDArray[np.float64]

# Default radial grid (fine; smooth integrands, exponential decay) -- matches
# the diagnostic driver so the exact target values reproduce.
_R_MAX = 60.0
_N_PTS = 40000


def radial_grid(r_max: float = _R_MAX, n_pts: int = _N_PTS) -> Tuple[FloatArray, float]:
    """Uniform radial grid on (0, r_max], avoiding the r=0 singularity."""
    r = np.linspace(1e-6, r_max, n_pts)
    return r, float(r[1] - r[0])


def _l2_normalize(f: FloatArray, r: FloatArray) -> FloatArray:
    """L2-normalize a radial function against the r^2 dr measure."""
    return f / np.sqrt(np.trapezoid(f * f * r * r, r))


def hydrogenic_s(n: int, Z: float, r: FloatArray) -> FloatArray:
    """Hydrogenic s radial R_{n0}(r), decay a = Z/n. L2-orthonormal across n."""
    a = Z / n
    L = genlaguerre(n - 1, 1)(2 * a * r)
    f = np.exp(-a * r) * L
    return _l2_normalize(f, r)


def sturmian_s(n: int, k: float, r: FloatArray) -> FloatArray:
    """Genuine shared-scale Coulomb-Sturmian s radial S_{n0}(r), fixed decay k for
    every n. L2-NON-orthogonal across n (only 1/r-weighted-orthogonal)."""
    L = genlaguerre(n - 1, 1)(2 * k * r)
    f = np.exp(-k * r) * L
    return _l2_normalize(f, r)


def build_basis(N: int, family: str, Z: float, r: FloatArray) -> List[FloatArray]:
    """Build the first N s-shell radial functions for `family` in {'hydrogenic',
    'sturmian'}, anchored so shell n=1 is identical between the two families
    (Sturmian scale k is fixed at Z)."""
    if family == "hydrogenic":
        return [hydrogenic_s(n, Z, r) for n in range(1, N + 1)]
    if family == "sturmian":
        return [sturmian_s(n, Z, r) for n in range(1, N + 1)]
    raise ValueError(f"unknown family {family!r}; expected 'hydrogenic' or 'sturmian'")


def overlap_matrix(fs: Sequence[FloatArray], r: FloatArray) -> FloatArray:
    """L2 overlap matrix S_ij = <f_i|f_j> of the given radial functions."""
    N = len(fs)
    S = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            S[i, j] = np.trapezoid(fs[i] * fs[j] * r * r, r)
    return S


def h1_matrix(fs: Sequence[FloatArray], r: FloatArray, Z: float) -> FloatArray:
    """One-body matrix h1_ij = <f_i| -1/2 d^2/dr^2 (l=0 kinetic) - Z/r |f_j> in the
    (non-orthogonal) radial basis fs."""
    N = len(fs)
    H = np.zeros((N, N))
    dfs = [np.gradient(f, r) for f in fs]
    for i in range(N):
        for j in range(N):
            T = 0.5 * np.trapezoid(dfs[i] * dfs[j] * r * r, r)
            V = -Z * np.trapezoid(fs[i] * fs[j] * r, r)
            H[i, j] = T + V
    return H


def _u_of(gdens: FloatArray, r: FloatArray, dr: float) -> FloatArray:
    """Classical (k=0 multipole) electrostatic potential of a density gdens(r):
    U(r1) = int_0^r1 gdens(r2) r2^2 dr2 / r1 + int_r1^inf gdens(r2) r2 dr2."""
    cum_in = np.cumsum(gdens * r * r) * dr
    tail = np.cumsum((gdens * r)[::-1])[::-1] * dr
    return cum_in / r + tail


def eri_tensor(fs: Sequence[FloatArray], r: FloatArray, dr: float) -> NDArray[np.float64]:
    """Chemist-notation s-only (k=0 multipole) two-electron repulsion integrals
    g[i,j,k,l] = (ij|kl) in the (non-orthogonal) radial basis fs."""
    N = len(fs)
    g = np.zeros((N, N, N, N))
    u_cache: Dict[Tuple[int, int], FloatArray] = {}
    for k_ in range(N):
        for l_ in range(N):
            u_cache[(k_, l_)] = _u_of(fs[k_] * fs[l_], r, dr)
    for i in range(N):
        for j in range(N):
            fij = fs[i] * fs[j]
            for k_ in range(N):
                for l_ in range(N):
                    g[i, j, k_, l_] = np.trapezoid(fij * u_cache[(k_, l_)] * r * r, r)
    return g


def lowdin_transform(S: FloatArray) -> FloatArray:
    """Symmetric (Loewdin) orthonormalization transform X = S^{-1/2} such that
    X^T S X = I."""
    w, V = np.linalg.eigh(S)
    return V @ np.diag(w ** -0.5) @ V.T


def jw_lcu_lambda(h1_orth: FloatArray, eri_orth: NDArray[np.float64]) -> float:
    """Faithful LCU 1-norm lambda = sum_i |c_i| (excluding the identity term) of
    the Jordan-Wigner-encoded restricted (spin-summed) molecular Hamiltonian built
    from h1_orth/eri_orth in an ORTHONORMAL spatial-orbital basis."""
    N = h1_orth.shape[0]
    H = FermionOperator()
    for p in range(N):
        for q in range(N):
            c = h1_orth[p, q]
            if abs(c) > 1e-12:
                for s in (0, 1):
                    H += FermionOperator(((2 * p + s, 1), (2 * q + s, 0)), c)
    for p in range(N):
        for q in range(N):
            for rr in range(N):
                for ss in range(N):
                    c = 0.5 * eri_orth[p, q, rr, ss]
                    if abs(c) > 1e-12:
                        for s1 in (0, 1):
                            for s2 in (0, 1):
                                H += FermionOperator(
                                    (
                                        (2 * p + s1, 1),
                                        (2 * rr + s2, 1),
                                        (2 * ss + s2, 0),
                                        (2 * q + s1, 0),
                                    ),
                                    c,
                                )
    jw = jordan_wigner(H)
    return float(sum(abs(c) for key, c in jw.terms.items() if key != ()))


def atomic_l2_lambda(
    N: int,
    family: str,
    Z: float = 2.0,
    r_max: float = _R_MAX,
    n_pts: int = _N_PTS,
) -> float:
    """The Paper 60 eq:blowup atomic L2-Loewdin Jordan-Wigner LCU 1-norm lambda
    for an N-shell s-only single-center model (He by default, Z=2), in the given
    radial family.

    family='hydrogenic': decay a=Z/n, L2-orthonormal across shells (Q^1.19 growth).
    family='sturmian':   shared decay k=Z for all shells, the genuine Coulomb-Sturmian
                          (L2-NON-orthogonal across shells; Q^3.33 growth after Loewdin).

    Returns lambda = sum |Pauli coefficient| over the non-identity Jordan-Wigner
    terms of the Loewdin-orthonormalized, spin-summed Hamiltonian.
    """
    r, dr = radial_grid(r_max, n_pts)
    fs = build_basis(N, family, Z, r)
    S = overlap_matrix(fs, r)
    X = lowdin_transform(S)
    h1 = h1_matrix(fs, r, Z)
    g4 = eri_tensor(fs, r, dr)
    h1o = X.T @ h1 @ X
    go = np.einsum("ip,jq,kr,ls,ijkl->pqrs", X, X, X, X, g4, optimize=True)
    return jw_lcu_lambda(h1o, go)


def validate_n1(Z: float = 2.0, r_max: float = _R_MAX, n_pts: int = _N_PTS) -> Tuple[float, float]:
    """Single-shell (N=1) sanity values: F0(1s,1s) (exact 5Z/8) and h1(1s)
    (exact -Z^2/2). Returns (F0, h1_11)."""
    r, dr = radial_grid(r_max, n_pts)
    f1 = [hydrogenic_s(1, Z, r)]
    g = eri_tensor(f1, r, dr)
    h = h1_matrix(f1, r, Z)
    return float(g[0, 0, 0, 0]), float(h[0, 0])


def lambda_sweep(
    N_values: Sequence[int],
    family: str,
    Z: float = 2.0,
    r_max: float = _R_MAX,
    n_pts: int = _N_PTS,
) -> List[Dict[str, float]]:
    """Sweep atomic_l2_lambda over a sequence of shell counts N, returning a list
    of {'N', 'Q', 'lam_excl'} records for `family`."""
    out: List[Dict[str, float]] = []
    for N in N_values:
        lam = atomic_l2_lambda(N, family, Z, r_max, n_pts)
        out.append({"N": float(N), "Q": float(2 * N), "lam_excl": lam})
    return out


def fit_lambda_exponent(
    N_values: Sequence[int],
    family: str,
    Z: float = 2.0,
    r_max: float = _R_MAX,
    n_pts: int = _N_PTS,
) -> float:
    """Fit lambda ~ Q^p (Q = 2N spin-orbitals) over the given shell counts by a
    log-log linear fit, excluding the degenerate N=1 point. Returns the fitted
    exponent p."""
    sweep = lambda_sweep([N for N in N_values if N > 1], family, Z, r_max, n_pts)
    xs = np.log([rec["Q"] for rec in sweep])
    ys = np.log([rec["lam_excl"] for rec in sweep])
    return float(np.polyfit(xs, ys, 1)[0])
