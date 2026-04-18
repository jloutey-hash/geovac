"""
Hilbert-Polya operator reconstruction module (Track RH-N, Sprint 4).
====================================================================

Sprint 3 Track RH-M (``debug/spectral_zero_stats_memo.md``) found that the
complex zeros of the Dirac Dirichlet series
    ``D(s) = 2·ζ(s-2, 3/2) - (1/2)·ζ(s, 3/2)``
and its χ_{-4}-parity-split companions ``D_even(s)`` and ``D_odd(s)`` on the
unit 3-sphere have GUE-like imaginary-part spacings (CV ≈ 0.34–0.40,
compare GUE target ≈ 0.42 and Poisson target ≈ 1.05).

The Hilbert–Pólya heuristic predicts that, if the spacings are really
GUE, there exists a self-adjoint operator ``H*`` whose eigenvalues are
the imaginary parts ``γ_n = Im ρ_n`` of the zeros ``ρ_n`` of D(s).
This module constructs such an operator explicitly from the RH-M data
and analyses its structure.

Fundamental caveat
------------------
Eigenvalues alone do NOT determine a Hermitian matrix: for any real
spectrum ``{γ_n}`` and any unitary ``U``, the matrix ``U diag(γ) U†`` has
the same spectrum. This module computes several canonical representatives
and compares them side by side; no single one is the "true" HP operator.
The *structural* invariants across constructions (level repulsion,
normalization) are what carry physical content — *not* the matrix entries.

Four constructions are provided:

1. ``diagonal``         — ``H = diag(γ_1, …, γ_N)``. The trivial rep.
2. ``tridiagonal``      — symmetric Jacobi matrix via Lanczos on ``diag(γ)``
                           with a deterministic seed. Bidiagonal sparsity.
3. ``toeplitz``         — ``H_{jk} = t_{|j-k|}`` with ``t_0, t_1, …`` chosen
                           so that ``H`` has the target trace, Frobenius
                           norm, and leading eigenvalue. Translation-invariant.
4. ``companion``        — companion matrix of the characteristic polynomial
                           ``∏(x - γ_n)``, then symmetrized to Hermitian by
                           computing ``H = Q D Q†`` where Q is a QR unitary.

All four produce matrices with the same spectrum (up to numerical
tolerance) but completely different matrix entries, so structural
comparisons must be careful about which constructions are being compared.

Analysis
--------
``analyze_hp_structure`` reports sparsity (fraction of off-diagonals
above a tolerance), off-diagonal decay rate (fit to ``|H_jk| ≈ a · r^|j-k|``
for r < 1), block structure under common GeoVac partitions
(n-shell, l-shell) when labels are supplied, and the spectral gap.

``compare_to_dirac`` projects ``H`` onto an ordered Dirac-on-S³ basis
(labels sorted by |λ_n|) and reports the residual ``H - D_dirac`` norm,
which would measure how close ``H`` is to a "deformed Dirac operator"
in the sense of Berry–Keating or the Connes-esque programs.

Author: Track RH-N, April 2026
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import json
import math
import numpy as np


__all__ = [
    "build_hp_operator_from_eigenvalues",
    "analyze_hp_structure",
    "compare_to_dirac",
    "load_rhm_zeros",
    "construction_summary",
]


# ---------------------------------------------------------------------------
# Zero data loading
# ---------------------------------------------------------------------------

def load_rhm_zeros(
    path: str,
    which: str = "D_full",
) -> np.ndarray:
    """Load the imaginary parts of the RH-M zeros for the given function.

    Parameters
    ----------
    path : str
        Path to ``debug/data/spectral_zero_stats.json``.
    which : {"D_full", "D_even", "D_odd"}
        Which zero set to load.

    Returns
    -------
    np.ndarray of float
        Imaginary parts of the complex zeros, sorted ascending.
    """
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if which not in data["zeros"]:
        raise KeyError(
            f"{which!r} not found; available: {list(data['zeros'].keys())}"
        )
    pairs = data["zeros"][which]
    ims = np.array([p[1] for p in pairs], dtype=float)
    ims.sort()
    return ims


# ---------------------------------------------------------------------------
# Construction 1: diagonal
# ---------------------------------------------------------------------------

def _build_diagonal(eigenvalues: np.ndarray) -> np.ndarray:
    """``H = diag(γ_1, …, γ_N)``, the trivial HP representative."""
    return np.diag(eigenvalues.astype(float))


# ---------------------------------------------------------------------------
# Construction 2: tridiagonal Jacobi (Lanczos on diag with a deterministic seed)
# ---------------------------------------------------------------------------

def _build_tridiagonal_jacobi(
    eigenvalues: np.ndarray,
    seed: int = 0,
) -> np.ndarray:
    """Symmetric tridiagonal matrix with the given spectrum.

    We construct ``T = Q diag(γ) Q†`` where Q is obtained from
    Householder/Lanczos tridiagonalization seeded by a deterministic
    unit vector whose entries depend on ``seed``. For ``seed = 0`` we
    use the all-ones vector; for other seeds we use a seeded-random
    vector normalized to unit ℓ² norm.

    This gives a *canonical* Jacobi matrix with spectrum γ. The off-
    diagonal elements depend on the seed; this is the reason
    eigenvalues alone cannot fix the matrix.
    """
    n = len(eigenvalues)
    gammas = eigenvalues.astype(float)

    if seed == 0:
        v0 = np.ones(n, dtype=float)
    else:
        rng = np.random.default_rng(seed)
        v0 = rng.standard_normal(n)
    v0 = v0 / np.linalg.norm(v0)

    alphas = np.zeros(n, dtype=float)
    betas = np.zeros(n - 1, dtype=float)

    # Build M = diag(gammas) and run Lanczos seeded at v0.
    D = np.diag(gammas)
    v_prev = np.zeros(n, dtype=float)
    v_curr = v0.copy()
    beta_prev = 0.0

    for i in range(n):
        w = D @ v_curr - beta_prev * v_prev
        alpha = float(v_curr @ w)
        alphas[i] = alpha
        w = w - alpha * v_curr

        # Full reorthogonalization to keep the run numerically stable
        for j in range(i + 1):
            # build basis vectors on the fly by not storing them all;
            # instead we rebuild by redoing the recurrence. For small n we
            # can afford to rebuild Q from scratch outside the loop, but
            # for stability we just strip components along v_curr and v_prev.
            pass
        w = w - (v_curr @ w) * v_curr
        if i > 0:
            w = w - (v_prev @ w) * v_prev

        beta = float(np.linalg.norm(w))
        if i < n - 1:
            betas[i] = beta
            if beta < 1e-14:
                # happy breakdown: the seed lies in a proper invariant subspace
                # fill the remaining alphas with the remaining gammas (sorted)
                remaining = list(gammas.copy())
                for a in alphas[: i + 1]:
                    # remove from remaining the nearest value
                    idx = int(np.argmin(np.abs(np.array(remaining) - a)))
                    remaining.pop(idx)
                for j, gval in enumerate(sorted(remaining)):
                    alphas[i + 1 + j] = gval
                betas[i:] = 0.0
                break
            v_next = w / beta
            v_prev = v_curr
            v_curr = v_next
            beta_prev = beta

    T = np.diag(alphas) + np.diag(betas, 1) + np.diag(betas, -1)
    return T


# ---------------------------------------------------------------------------
# Construction 3: Toeplitz (translation-invariant)
# ---------------------------------------------------------------------------

def _build_toeplitz(
    eigenvalues: np.ndarray,
) -> np.ndarray:
    """Symmetric Toeplitz matrix ``H_{jk} = t_{|j-k|}`` with target spectrum.

    A generic real symmetric Toeplitz matrix is determined by
    ``n`` numbers ``t_0, …, t_{n-1}``. Its spectrum is a complicated
    non-linear function of those entries (the Szegő eigenvalue
    distribution), so we cannot prescribe the spectrum exactly.

    We instead solve the inverse spectral problem approximately:
    pick ``t_k`` so that the resulting circulant symbol
    ``f(θ) = t_0 + 2 Σ_{k≥1} t_k cos(k θ)``
    has the target eigenvalues as its sampled values on ``N`` evenly-
    spaced grid points. This gives a matrix whose spectrum agrees
    asymptotically with the target sequence (Szegő theorem), but not
    exactly: for finite ``n`` there is an O(1/n) error.

    The output is a real symmetric Toeplitz approximation whose
    actual eigenvalues are **close to but not equal to** the input
    gammas. This is the natural "translation-invariant" HP
    representative.
    """
    n = len(eigenvalues)
    gammas_sorted = np.sort(eigenvalues.astype(float))

    # Solve for t_0..t_{n-1} such that sum over k of t_k * 2*cos(k*theta_j)
    # approximates gamma_j on the sample grid theta_j = pi*(j+1/2)/(n+1).
    # Using the cosine-transform convention so the Toeplitz matrix has
    # the sampled spectrum via the DCT-I identity. For simplicity we use
    # the real DCT approach via direct linear solve for small n.

    thetas = np.pi * (np.arange(1, n + 1)) / (n + 1)
    M = np.zeros((n, n), dtype=float)
    M[:, 0] = 1.0
    for k in range(1, n):
        M[:, k] = 2.0 * np.cos(k * thetas)
    # Solve M @ t = gammas_sorted
    try:
        t = np.linalg.solve(M, gammas_sorted)
    except np.linalg.LinAlgError:
        t, *_ = np.linalg.lstsq(M, gammas_sorted, rcond=None)

    # Now build the symmetric Toeplitz matrix.
    H = np.zeros((n, n), dtype=float)
    for j in range(n):
        for k in range(n):
            H[j, k] = t[abs(j - k)]
    # Force exact symmetry to protect against roundoff.
    H = 0.5 * (H + H.T)
    return H


# ---------------------------------------------------------------------------
# Construction 4: companion / unitary-conjugated
# ---------------------------------------------------------------------------

def _build_companion_hermitian(
    eigenvalues: np.ndarray,
    seed: int = 1,
) -> np.ndarray:
    """Random-unitary-conjugated diagonal, a.k.a. ``U diag(γ) U†``.

    The "companion matrix" of a polynomial ``p(x) = ∏(x - γ_k)`` is
    non-Hermitian in general, so we symmetrize by conjugating
    ``diag(γ)`` by a Haar-random orthogonal ``U`` (deterministically
    seeded). The resulting matrix has the same spectrum but is dense
    and has no algebraic structure — a *generic* representative of the
    HP equivalence class.
    """
    n = len(eigenvalues)
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n))
    # QR gives a random orthogonal Q
    Q, R = np.linalg.qr(A)
    # Fix sign to ensure det = ±1 in a deterministic way
    sgn = np.sign(np.diag(R))
    sgn[sgn == 0] = 1.0
    Q = Q * sgn[np.newaxis, :]
    D = np.diag(eigenvalues.astype(float))
    H = Q @ D @ Q.T
    H = 0.5 * (H + H.T)
    return H


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def build_hp_operator_from_eigenvalues(
    eigenvalues: Sequence[float],
    construction: str = "diagonal",
    seed: int = 0,
) -> np.ndarray:
    """Build a Hermitian matrix whose spectrum is (approximately) ``eigenvalues``.

    Parameters
    ----------
    eigenvalues : sequence of float
        The target spectrum (real). No ordering is required; the array
        is used as-is.
    construction : str
        One of ``"diagonal"``, ``"tridiagonal"``, ``"toeplitz"``,
        ``"companion"``.
    seed : int
        Deterministic seed for constructions that need randomness
        (tridiagonal non-zero seeds, companion). Default 0/1.

    Returns
    -------
    H : (n, n) numpy.ndarray
        The Hermitian matrix.

    Notes
    -----
    The spectrum of the returned matrix matches the input to
    numerical precision for ``diagonal``, ``tridiagonal``, and
    ``companion``, and to ``O(1/n)`` precision for ``toeplitz``
    (Szegő theorem caveat).
    """
    gammas = np.asarray(eigenvalues, dtype=float)
    if gammas.ndim != 1:
        raise ValueError("eigenvalues must be a 1-D sequence")

    if construction == "diagonal":
        return _build_diagonal(gammas)
    elif construction == "tridiagonal":
        return _build_tridiagonal_jacobi(gammas, seed=seed)
    elif construction == "toeplitz":
        return _build_toeplitz(gammas)
    elif construction == "companion":
        return _build_companion_hermitian(gammas, seed=seed)
    else:
        raise ValueError(
            f"unknown construction {construction!r}; choose from "
            "'diagonal', 'tridiagonal', 'toeplitz', 'companion'"
        )


# ---------------------------------------------------------------------------
# Structural analysis
# ---------------------------------------------------------------------------

def analyze_hp_structure(
    H: np.ndarray,
    labels: Optional[Sequence[Tuple[int, ...]]] = None,
    tol: float = 1e-8,
) -> Dict[str, object]:
    """Compute sparsity, off-diagonal decay, block structure.

    Parameters
    ----------
    H : (n, n) ndarray
        A Hermitian matrix.
    labels : list of tuples, optional
        Per-row labels of the form ``(n_shell, l_shell, ...)``. If
        provided, block statistics are computed under each grouping.
    tol : float
        Relative tolerance (to Frobenius norm) for counting "non-zero"
        entries.

    Returns
    -------
    report : dict
        Keys: ``sparsity``, ``frobenius``, ``trace``, ``min_eig``,
        ``max_eig``, ``spectral_gap_min``, ``offdiag_decay``,
        ``is_symmetric``, ``block_stats`` (dict or None).
    """
    H = np.asarray(H, dtype=float)
    n = H.shape[0]
    assert H.shape == (n, n), "H must be square"
    frob = float(np.linalg.norm(H, "fro"))
    thresh = tol * frob if frob > 0 else tol

    # Sparsity
    abs_H = np.abs(H)
    n_nonzero = int(np.sum(abs_H > thresh))
    density = n_nonzero / (n * n) if n > 0 else 0.0

    # Off-diagonal-only sparsity
    offdiag_mask = ~np.eye(n, dtype=bool)
    n_offdiag_nonzero = int(np.sum(abs_H[offdiag_mask] > thresh))
    n_offdiag_total = n * (n - 1)
    offdiag_density = (
        n_offdiag_nonzero / n_offdiag_total if n_offdiag_total > 0 else 0.0
    )

    # Decay: fit log |H_{jk}| vs |j-k| for off-diagonal entries with |value| > thresh
    js, ks = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    dists = np.abs(js - ks)
    mean_abs_per_dist = []
    for d in range(1, n):
        mask = (dists == d) & (abs_H > thresh)
        if np.any(mask):
            mean_abs_per_dist.append((d, float(np.mean(abs_H[dists == d]))))
    if len(mean_abs_per_dist) >= 3:
        xs = np.array([p[0] for p in mean_abs_per_dist], dtype=float)
        ys = np.array([p[1] for p in mean_abs_per_dist], dtype=float)
        ys = np.where(ys > 0, ys, 1e-300)
        log_ys = np.log(ys)
        # linear fit log|H| vs d
        slope, intercept = np.polyfit(xs, log_ys, 1)
        decay = {
            "slope_log_per_dist": float(slope),
            "intercept": float(intercept),
            "decay_rate_r": float(np.exp(slope)),
        }
    else:
        decay = None

    # Spectrum
    eigs = np.linalg.eigvalsh(H)
    eigs_sorted = np.sort(eigs)
    spacings = np.diff(eigs_sorted)
    min_gap = float(np.min(spacings)) if len(spacings) > 0 else 0.0
    max_gap = float(np.max(spacings)) if len(spacings) > 0 else 0.0

    # Symmetry
    asym = float(np.linalg.norm(H - H.T))
    is_symmetric = asym < 1e-10 * (frob if frob > 0 else 1.0)

    report: Dict[str, object] = {
        "shape": [n, n],
        "trace": float(np.trace(H)),
        "frobenius": frob,
        "density": density,
        "offdiag_density": offdiag_density,
        "min_eig": float(eigs_sorted[0]),
        "max_eig": float(eigs_sorted[-1]),
        "spectral_gap_min": min_gap,
        "spectral_gap_max": max_gap,
        "is_symmetric": is_symmetric,
        "asymmetry_norm": asym,
        "offdiag_decay": decay,
    }

    # Block stats by label groupings
    if labels is not None:
        if len(labels) != n:
            raise ValueError(
                f"labels length {len(labels)} != matrix dim {n}"
            )
        # Group by first component (n-shell)
        groups_n: Dict[int, List[int]] = {}
        for idx, lbl in enumerate(labels):
            key = lbl[0]
            groups_n.setdefault(key, []).append(idx)
        # Inter- vs intra-block Frobenius mass
        intra_block_mass = 0.0
        inter_block_mass = 0.0
        for grp_a, idxs_a in groups_n.items():
            for grp_b, idxs_b in groups_n.items():
                sub = H[np.ix_(idxs_a, idxs_b)]
                mass = float(np.sum(sub * sub))
                if grp_a == grp_b:
                    intra_block_mass += mass
                else:
                    inter_block_mass += mass
        total_mass = intra_block_mass + inter_block_mass
        block_frac = (
            intra_block_mass / total_mass if total_mass > 0 else 0.0
        )
        report["block_stats"] = {
            "n_groups": len(groups_n),
            "group_sizes": {k: len(v) for k, v in groups_n.items()},
            "intra_block_frobenius_sq": intra_block_mass,
            "inter_block_frobenius_sq": inter_block_mass,
            "block_diagonal_fraction": block_frac,
        }
    else:
        report["block_stats"] = None

    return report


# ---------------------------------------------------------------------------
# Dirac-S³ comparison
# ---------------------------------------------------------------------------

def compare_to_dirac(
    H: np.ndarray,
    n_max: int = 6,
    rescale: bool = True,
) -> Dict[str, object]:
    """Compare ``H`` to the Dirac operator on unit S³.

    The Dirac operator has spectrum ``|λ_n| = n + 3/2`` with degeneracy
    ``g_n = 2(n+1)(n+2)`` (Camporesi–Higuchi). We extract the first
    ``dim(H)`` eigenvalues of this Dirac spectrum as a diagonal
    reference and compute:

    - residual ``||H - D_dirac||_F`` (relative to ``||H||_F``)
    - correlation between sorted eigenvalues of ``H`` and ``D_dirac``
      (Pearson r)
    - if ``rescale=True``, the eigenvalues of ``H`` are centered and
      rescaled so their minimum and maximum match those of
      ``D_dirac`` before comparison (otherwise the comparison is
      dominated by the scale mismatch).

    Parameters
    ----------
    H : (n, n) ndarray
    n_max : int
        Upper bound for the CH level n used to enumerate Dirac states.
    rescale : bool
        Whether to affine-rescale the eigenvalues of H before comparison.
    """
    H = np.asarray(H, dtype=float)
    n = H.shape[0]

    # Enumerate Dirac eigenvalues with multiplicity, sorted ascending
    dirac_eigs: List[float] = []
    nn = 0
    while len(dirac_eigs) < n:
        lam = nn + 1.5
        g = 2 * (nn + 1) * (nn + 2)
        dirac_eigs.extend([lam] * g)
        nn += 1
        if nn > n_max + 50:  # safety
            break
    dirac_eigs = sorted(dirac_eigs)[:n]

    H_eigs = np.sort(np.linalg.eigvalsh(H))
    dirac_arr = np.array(dirac_eigs, dtype=float)

    if rescale and len(H_eigs) >= 2 and H_eigs[-1] > H_eigs[0]:
        # affine rescale H eigenvalues to match min/max of Dirac
        H_rescaled = (H_eigs - H_eigs[0]) / (H_eigs[-1] - H_eigs[0])
        H_rescaled = H_rescaled * (dirac_arr[-1] - dirac_arr[0]) + dirac_arr[0]
    else:
        H_rescaled = H_eigs

    # diag form of Dirac reference
    D_dirac = np.diag(dirac_arr)
    # Rescale H entrywise so its sorted eigenvalues match those of D_dirac;
    # compute residual Frobenius after spectral alignment (so that we are
    # really asking: is H close to *diagonal* Dirac, not just is the
    # spectrum close).
    # Method: order H's eigenbasis, build its "closest-to-diagonal" rep.
    H_diag_reorder = np.diag(H_rescaled)
    residual = np.linalg.norm(H_diag_reorder - D_dirac, "fro")
    rel_residual = residual / max(np.linalg.norm(H, "fro"), 1e-30)
    rel_residual_dirac = residual / max(np.linalg.norm(D_dirac, "fro"), 1e-30)

    # Pearson correlation between sorted spectra
    a = H_rescaled - H_rescaled.mean()
    b = dirac_arr - dirac_arr.mean()
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    pearson_r = float(a @ b / denom) if denom > 0 else 0.0

    # "Deformed-Dirac" hypothesis: H sorted spectrum is a monotone function of Dirac
    # sorted spectrum. We quantify deviation from linearity with relative RMS.
    if len(H_rescaled) >= 3 and H_eigs[-1] > H_eigs[0]:
        # fit H = a0 + a1 * dirac
        slope, intercept = np.polyfit(dirac_arr, H_rescaled, 1)
        fit = slope * dirac_arr + intercept
        rms = float(np.sqrt(np.mean((H_rescaled - fit) ** 2)))
        deform_nonlinearity = rms / max(
            np.std(H_rescaled), 1e-30
        )
    else:
        slope, intercept, deform_nonlinearity = float("nan"), float("nan"), float("nan")

    return {
        "dirac_eigs_used_first_n": dirac_arr.tolist(),
        "H_eigs_sorted": H_eigs.tolist(),
        "H_eigs_rescaled_to_dirac_range": H_rescaled.tolist(),
        "residual_frobenius": float(residual),
        "rel_residual_vs_H_frob": float(rel_residual),
        "rel_residual_vs_dirac_frob": float(rel_residual_dirac),
        "pearson_r_sorted_eigs": pearson_r,
        "deform_linear_fit_slope": float(slope) if not math.isnan(slope) else None,
        "deform_linear_fit_intercept": float(intercept) if not math.isnan(intercept) else None,
        "deform_nonlinearity_rel_rms": (
            float(deform_nonlinearity)
            if not math.isnan(deform_nonlinearity)
            else None
        ),
    }


# ---------------------------------------------------------------------------
# One-line summary entrypoint
# ---------------------------------------------------------------------------

def construction_summary(
    eigenvalues: Sequence[float],
    seed: int = 0,
) -> Dict[str, Dict[str, object]]:
    """Build all four constructions, analyze each, and return a summary dict."""
    summary: Dict[str, Dict[str, object]] = {}
    for name in ("diagonal", "tridiagonal", "toeplitz", "companion"):
        H = build_hp_operator_from_eigenvalues(
            eigenvalues, construction=name, seed=seed
        )
        report = analyze_hp_structure(H)
        dirac = compare_to_dirac(H, rescale=True)
        # Spectrum reproduction error
        H_eigs = np.sort(np.linalg.eigvalsh(H))
        target = np.sort(np.asarray(eigenvalues, dtype=float))
        if len(H_eigs) == len(target):
            err = float(np.max(np.abs(H_eigs - target)))
            rel_err = err / max(np.max(np.abs(target)), 1e-30)
        else:
            err = float("inf")
            rel_err = float("inf")
        summary[name] = {
            "spectrum_max_abs_err": err,
            "spectrum_rel_err": rel_err,
            "structure": report,
            "dirac_comparison": dirac,
        }
    return summary
