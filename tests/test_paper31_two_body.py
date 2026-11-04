"""Paper 31 §sec:two_body — two-body verification of the A/D partition.

Pins the published numbers in Paper 31's "Two-body verification of the
partition" section (`papers/group3_foundations/paper_31_universal_coulomb_partition.tex`,
\\label{sec:two_body}). These originate as honest-negative diagnostics
(CLAUDE.md §3 dead-end rows: "Gauged tensor-product spectral action for
two-body Coulomb weights", 2026-06-03; "Resolvent (D^2)^-1 for two-body
Coulomb interaction", 2026-06-01). The paper presents them as a
*verification* of the partition, so this file gives them a backing test
(closes the C1 NO-TEST gap surfaced by `/qa group3`).

Two constructions are pinned:

  1. Gauged spectral action on T_{S^3} (x) T_{S^3} (double-sum inner
     fluctuation A = sum_i sum_j (M_i^dag (x) M_j^dag)[D_total, M_i (x) M_j]):
       - connected (non-factorizable) fraction of {D, A}: 77% (n_max=2),
         32% (n_max=3)
       - support (Paper 54 thm:angular_structure): EXACTLY the direct
         elements V[a,b,a,b] -- all N^2 of them nonzero (25 at n_max=2,
         196 at n_max=3), nothing else; so m-conservation, Delta l = 0 and
         the pure k=0 monopole hold automatically.  The Coulomb tensor has
         82 / 3604 further nonzero elements (exchange and k>0 channels
         included) that are absent from the spectral support, and no
         spectral element lies outside the Coulomb support.
         (Replaces a Frobenius-weighted "Gaunt / k=0" audit that passed for
         ANY direct-diagonal tensor -- planting `if True:` in its Gaunt
         branch did not fire; /qa group3 2026-11-03.)
       - radial Pearson vs exact Coulomb Slater: 0.58 (n_max=2),
         0.41 (n_max=3), decreasing with n_max

  Paper 54 thm:free_factorization + tab:heat are pinned here too:
  {D, gamma} = 0 exactly, the free heat-trace ratio is 1 at
  t in {0.1, 0.5, 1, 2}, 4 of 14 / 29 of 55 multipliers have a nonzero
  [D, M], and the gauged heat-trace ratios of tab:heat (t = 0.5) match to
  the printed six decimals.

  2. Resolvent construction V^res_abcd = sum_NLM M[a,c] conj(M[d,b]) w(N):
       - Dirac resolvent w=1/(N+1/2)^2: Pearson 0.846 (n_max=2),
         0.796 (n_max=3) (0.81/0.75 retired by the 2026-08-29 ERI fix)
       - Laplacian resolvent w=1/(N^2-1): essentially uncorrelated
         (the N=1 zero mode of Delta_{S^3} kills the dominant (1s,1s) term)
       - 100% m-conservation for every weighting

The point of the partition is the *negative*: the angular (Gaunt) part is
shared with Coulomb exactly (100% selection-rule match) while the radial
part is structurally different (Pearson well below 1, decreasing with
n_max). The test pins both the positive angular result and the negative
radial result.

The construction is replicated inline from production `geovac/` modules
(operator_system, casimir_ci, hypergeometric_slater) so the test depends
only on the shipped package, not on the `debug/` diagnostic drivers
`paper54_recompute_both_constructions.py` / `resolvent_two_body_diagnostic.py`
(which it bit-reproduces).
"""
import functools

import numpy as np
import pytest

from geovac.operator_system import (
    HyperLabel,
    build_multiplier_matrix,
    allowed_multiplier_labels,
)
from geovac.casimir_ci import two_electron_integral


# ---------------------------------------------------------------------------
# Shared construction helpers (replicated from the verified debug/ drivers,
# using only production geovac/ modules)
# ---------------------------------------------------------------------------


def _scalar_basis(n_max):
    """Fock-projected scalar basis |n, l, m> up to n_max."""
    return [
        HyperLabel(n=n, l=l, m=m)
        for n in range(1, n_max + 1)
        for l in range(n)
        for m in range(-l, l + 1)
    ]


def _coulomb_tensor(basis):
    """Exact Coulomb two-electron integral tensor in the scalar basis."""
    N = len(basis)
    C = np.zeros((N, N, N, N))
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    val = two_electron_integral(
                        basis[a].n, basis[a].l, basis[a].m,
                        basis[b].n, basis[b].l, basis[b].m,
                        basis[c].n, basis[c].l, basis[c].m,
                        basis[d].n, basis[d].l, basis[d].m,
                    )
                    if abs(val) > 1e-15:
                        C[a, b, c, d] = val
    return C


def _radial_pearson(V, C):
    """Pearson correlation of a candidate interaction tensor vs Coulomb,
    over the matrix elements where BOTH are nonzero.

    This is the spectral-action driver's `radial_vs_slater` convention
    (`debug/paper54_recompute_both_constructions.py`)."""
    s, c = V.ravel(), C.ravel()
    both = (np.abs(s) > 1e-12) & (np.abs(c) > 1e-12)
    if int(both.sum()) < 3:
        return None
    ss, cc = s[both], c[both]
    if np.std(ss) < 1e-15 or np.std(cc) < 1e-15:
        return 0.0
    return float(np.corrcoef(ss, cc)[0, 1])


def _resolvent_pearson(V, C):
    """Pearson over ALL Coulomb-nonzero elements (the candidate value may be
    zero there).

    This is the resolvent driver's `compare_interactions` convention
    (`debug/resolvent_two_body_diagnostic.py`) — it differs from the
    spectral-action both-nonzero convention because the resolvent's coverage
    of the Coulomb support (~50-59%) is itself part of the diagnostic."""
    s, c = V.ravel(), C.ravel()
    mask = np.abs(c) > 1e-14
    if int(mask.sum()) < 3:
        return None
    ss, cc = s[mask], c[mask]
    if np.std(ss) < 1e-15 or np.std(cc) < 1e-15:
        return 0.0
    return float(np.corrcoef(ss, cc)[0, 1])


# ---------------------------------------------------------------------------
# 1. Gauged spectral action (double-sum inner fluctuation)
# ---------------------------------------------------------------------------


def _single_particle_scalar(n_max):
    """Chirality-doubled single-particle scalar triple: D = +-(n+1/2),
    gamma swaps the two chirality copies."""
    basis = _scalar_basis(n_max)
    N = len(basis)
    dim = 2 * N
    D = np.zeros((dim, dim))
    for i, b in enumerate(basis):
        ev = b.n + 0.5
        D[i, i] = +ev
        D[N + i, N + i] = -ev
    gamma = np.zeros((dim, dim))
    gamma[:N, N:] = np.eye(N)
    gamma[N:, :N] = np.eye(N)
    return basis, N, dim, D, gamma


def _embedded_multipliers(basis, dim, N, n_max):
    """3-Y multipliers block-embedded in the chirality-doubled space."""
    out = []
    for (Nq, L, M) in allowed_multiplier_labels(n_max):
        M_sc = build_multiplier_matrix(Nq, L, M, basis)
        if np.linalg.norm(M_sc) < 1e-15:
            continue
        M_full = np.zeros((dim, dim), dtype=complex)
        M_full[:N, :N] = M_sc
        M_full[N:, N:] = M_sc
        out.append(M_full)
    return out


def _commutator_count(n_max):
    """(number of multipliers with nonzero [D, M], number of multipliers)."""
    basis, N, dim, D1, _ = _single_particle_scalar(n_max)
    mults = _embedded_multipliers(basis, dim, N, n_max)
    nz = sum(1 for M in mults if np.linalg.norm(D1 @ M - M @ D1) > 1e-14)
    return nz, len(mults)


@functools.lru_cache(maxsize=None)
def _gauged_build(n_max):
    """Double-sum gauge field. Returns (D_total, A, D1, basis, N, dim).

    A = sum_i sum_j (M_i^dag (x) M_j^dag) [D_total, M_i (x) M_j]
    D_total = D_1 (x) I + gamma_1 (x) D_2
    Cached: the n_max = 3 build is ~45 s and is shared by several tests.
    """
    basis, N, dim, D1, g1 = _single_particle_scalar(n_max)
    d = dim
    D_total = np.kron(D1, np.eye(d)) + np.kron(g1, D1)
    mults = _embedded_multipliers(basis, dim, N, n_max)

    A = np.zeros((d * d, d * d), dtype=complex)
    for Mi in mults:
        ci = D1 @ Mi - Mi @ D1
        if np.linalg.norm(ci) < 1e-14:
            continue
        for Mj in mults:
            cj = D1 @ Mj - Mj @ D1
            if np.linalg.norm(cj) < 1e-14:
                continue
            fluct = np.kron(ci, Mj) + np.kron(g1 @ Mi, cj)
            A += np.kron(Mi.conj().T, Mj.conj().T) @ fluct
    A = ((A + A.conj().T) / 2).real
    return D_total, A, D1, basis, N, d


def _heat_ratio(D_total, A, D1, eps, t):
    """Tr e^{-t (D_total + eps A)^2} / (Tr e^{-t D_1^2})^2 (tab:heat)."""
    z1 = float(np.sum(np.exp(-t * np.linalg.eigvalsh(D1) ** 2)))
    lam = np.linalg.eigvalsh(D_total + eps * A)
    return float(np.sum(np.exp(-t * lam ** 2))) / (z1 * z1)


def _spectral_action_connected(n_max):
    """Return (connected_fraction, connected_pp_tensor, basis) of
    {D_total, A} for the double-sum gauge field."""
    D_total, A, _, basis, N, d = _gauged_build(n_max)

    DA = D_total @ A + A @ D_total

    # connected (non-factorizable) part via partial-trace subtraction
    I = np.eye(d)
    t = DA.reshape(d, d, d, d)
    Tr2 = np.trace(t, axis1=1, axis2=3)
    Tr1 = np.trace(t, axis1=0, axis2=2)
    fac = (np.kron(Tr2 / d, I) + np.kron(I, Tr1 / d)
           - (np.trace(DA) / (d * d)) * np.eye(d * d))
    conn = DA - fac
    frac = float(np.linalg.norm(conn) / np.linalg.norm(DA))

    # project the connected part onto the (a,b,c,d) scalar tensor
    # ordering: full index = scalar_idx * dim + scalar_idx (top-left chirality
    # block carries the scalar content)
    V = np.zeros((N, N, N, N))
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for dd in range(N):
                    V[a, b, c, dd] = conn[a * d + b, c * d + dd].real
    return frac, V, basis


SUPPORT_TOL = 1e-12  # smallest direct |V| is ~1e-2; largest off-direct ~1e-16


def _support_audit(V, C, basis):
    """Support of V against the direct set {[a,b,a,b]} and against Coulomb C.

    Returns a dict of integer counts.  This is the content of Paper 54
    thm:angular_structure: the spectral support IS the direct set, and the
    Coulomb-only elements (82 / 3604) are absent from it.
    """
    N = len(basis)
    i = np.arange(N)
    direct = np.zeros(V.shape, dtype=bool)
    direct[i[:, None], i[None, :], i[:, None], i[None, :]] = True
    spec = np.abs(V) > SUPPORT_TOL
    coul = np.abs(C) > SUPPORT_TOL
    coul_only = coul & ~spec
    l = np.array([b.l for b in basis])
    # Coulomb-only elements in a channel the spectral V cannot reach:
    # exchange [a,b,b,a] (a != b) and a multipole k > 0 (l_a != l_c)
    exch = np.zeros(V.shape, dtype=bool)
    exch[i[:, None], i[None, :], i[None, :], i[:, None]] = True
    exch &= ~direct
    lac_differs = (l[:, None, None, None] != l[None, None, :, None])
    return {
        "n_spec": int(spec.sum()),
        "n_direct": int(direct.sum()),
        "n_direct_nonzero": int((spec & direct).sum()),
        "n_off_direct": int((spec & ~direct).sum()),
        "n_coul": int(coul.sum()),
        "n_coul_only": int(coul_only.sum()),
        "n_spec_only": int((spec & ~coul).sum()),
        "n_coul_only_exchange": int((coul_only & exch).sum()),
        "n_coul_only_k_pos": int((coul_only & lac_differs).sum()),
    }


def _assert_direct_support(V, C, basis, n_paper, n_coul_only_paper):
    s = _support_audit(V, C, basis)
    N = len(basis)
    # support is EXACTLY the direct set: every direct element nonzero ...
    assert s["n_direct"] == N * N == n_paper, s
    assert s["n_direct_nonzero"] == n_paper, s
    # ... and nothing else (a non-direct element, even m-conserving with
    # Delta l = 0, is the wrong answer the retired audit accepted)
    assert s["n_off_direct"] == 0, s
    assert s["n_spec"] == n_paper, s
    # strict subset of Coulomb: 0 spectral-only, 82 / 3604 Coulomb-only
    assert s["n_spec_only"] == 0, s
    assert s["n_coul_only"] == n_coul_only_paper, s
    assert s["n_coul"] == n_paper + n_coul_only_paper, s
    # the Coulomb-only set contains exchange and k > 0 channels (the remark
    # after thm:angular_structure)
    assert s["n_coul_only_exchange"] > 0 and s["n_coul_only_k_pos"] > 0, s


def test_spectral_action_two_body_nmax2():
    """Paper 31 §sec:two_body, gauged spectral action at n_max=2.

    connected fraction ~ 77%; support exactly the 25 direct elements
    (Paper 54 thm:angular_structure), 82 Coulomb-only elements absent;
    radial Pearson ~ 0.58.
    """
    frac, V, basis = _spectral_action_connected(2)
    C = _coulomb_tensor(basis)
    pearson = _radial_pearson(V, C)

    # connected (non-factorizable) fraction: paper says 77%
    assert frac == pytest.approx(0.767, abs=0.01), f"conn frac {frac}"
    # angular structure: direct-only support (implies m-conservation,
    # Delta l = 0 and the pure k = 0 monopole)
    _assert_direct_support(V, C, basis, n_paper=25, n_coul_only_paper=82)
    # radial: structurally different from Coulomb (Pearson well below 1)
    assert pearson == pytest.approx(0.58, abs=0.02), f"Pearson {pearson}"
    assert pearson < 0.9


@pytest.mark.slow
def test_spectral_action_two_body_nmax3_decreasing():
    """Paper 31 §sec:two_body at n_max=3: connected fraction ~ 32%,
    radial Pearson ~ 0.41 and DECREASING vs n_max=2 (~0.58).

    Slow: the double-sum inner fluctuation is an O(n_multipliers^2) sum of
    16x16 (x) 16x16 Kronecker products plus a 256x256 eigen-free assembly.
    """
    frac2, V2, basis2 = _spectral_action_connected(2)
    pearson2 = _radial_pearson(V2, _coulomb_tensor(basis2))

    frac3, V3, basis3 = _spectral_action_connected(3)
    C3 = _coulomb_tensor(basis3)
    pearson3 = _radial_pearson(V3, C3)

    assert frac3 == pytest.approx(0.324, abs=0.01), f"conn frac {frac3}"
    _assert_direct_support(V3, C3, basis3, n_paper=196, n_coul_only_paper=3604)
    assert pearson3 == pytest.approx(0.41, abs=0.02), f"Pearson {pearson3}"
    # the load-bearing negative: radial correlation DECREASES with n_max
    assert pearson3 < pearson2


# ---------------------------------------------------------------------------
# 1b. Paper 54 thm:free_factorization, the commutator counts, tab:heat
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [2, 3])
def test_free_factorization_anticommutator_and_heat_ratio(n_max):
    """Paper 54 thm:free_factorization: {D, gamma} = 0 exactly for the
    chirality-doubled CH Dirac and the swap grading; hence D_total^2 has no
    cross term and the free heat-trace ratio is 1 at t in {0.1, 0.5, 1, 2}
    ("Numerical verification" paragraph).

    Wrong answer rejected: a non-chiral D (both copies +(n+1/2)), for which
    the swap grading commutes with D, the cross term 2 D gamma (x) D
    survives, and the heat trace does not factorize.
    """
    basis, N, dim, D, gamma = _single_particle_scalar(n_max)
    anti = D @ gamma + gamma @ D
    assert np.max(np.abs(anti)) == 0.0
    # gamma is a grading: involutive and not the identity
    assert np.array_equal(gamma @ gamma, np.eye(dim))
    assert not np.array_equal(gamma, np.eye(dim))

    D_total = np.kron(D, np.eye(dim)) + np.kron(gamma, D)
    D2_sum = np.kron(D @ D, np.eye(dim)) + np.kron(np.eye(dim), D @ D)
    assert np.max(np.abs(D_total @ D_total - D2_sum)) < 1e-12

    zero = np.zeros_like(D_total)
    for t in (0.1, 0.5, 1.0, 2.0):
        r = _heat_ratio(D_total, zero, D, 0.0, t)
        assert r == pytest.approx(1.0, abs=1e-12), (t, r)


@pytest.mark.parametrize("n_max,n_nonzero,n_total", [(2, 4, 14), (3, 29, 55)])
def test_commutator_counts(n_max, n_nonzero, n_total):
    """Paper 54 sec:gauged: at n_max = 2 four of 14 multiplier generators
    have [D, M] != 0; at n_max = 3, 29 of 55.

    Wrong answer rejected: a D that is not constant on n-shells (e.g. graded
    by l), which changes which multipliers commute with it.
    """
    nz, total = _commutator_count(n_max)
    assert (nz, total) == (n_nonzero, n_total)


# tab:heat, t = 0.5, printed to six decimals
TAB_HEAT = {
    2: {0.0: 1.000000, 0.1: 0.995465, 0.5: 0.977815, 1.0: 0.956832},
    3: {0.0: 1.000000, 0.1: 0.838438, 0.5: 0.410165, 1.0: 0.185405},
}


def _check_heat_table(n_max):
    D_total, A, D1, _, _, _ = _gauged_build(n_max)
    got = {eps: _heat_ratio(D_total, A, D1, eps, 0.5) for eps in TAB_HEAT[n_max]}
    for eps, ref in TAB_HEAT[n_max].items():
        # half a unit in the sixth decimal (display rounding)
        assert got[eps] == pytest.approx(ref, abs=5.01e-7), (n_max, eps, got[eps])
    # "departs from the factorized value monotonically with eps"
    vals = [got[e] for e in sorted(got)]
    assert all(x > y for x, y in zip(vals, vals[1:])), vals


def test_heat_table_nmax2():
    """Paper 54 tab:heat, n_max = 2 column, for the DOUBLE-sum gauge field.

    Wrong answer rejected: the single-sum (diagonal NLM = N'L'M') gauge
    field, which gives 0.998862 / 0.994340 / 0.988757.
    """
    _check_heat_table(2)


@pytest.mark.slow
def test_heat_table_nmax3():
    """Paper 54 tab:heat, n_max = 3 column (slow: ~45 s double-sum build,
    shared with test_spectral_action_two_body_nmax3_decreasing)."""
    _check_heat_table(3)


# ---------------------------------------------------------------------------
# 2. Resolvent construction
# ---------------------------------------------------------------------------


def _resolvent_interaction(n_max, weight_fn):
    """V^res_abcd = sum_NLM M[a,c] conj(M[d,b]) w(N) (real part)."""
    basis = _scalar_basis(n_max)
    N_dim = len(basis)
    V = np.zeros((N_dim, N_dim, N_dim, N_dim), dtype=np.complex128)
    for (Nq, L, M) in allowed_multiplier_labels(n_max):
        w = weight_fn(Nq)
        if abs(w) < 1e-30:
            continue
        M_mat = build_multiplier_matrix(Nq, L, M, basis)
        if np.linalg.norm(M_mat) < 1e-15:
            continue
        M_dag = M_mat.conj().T
        # V[a,b,c,d] += M[a,c] * M_dag[b,d] * w
        # = einsum over (a,c) outer (b,d)
        V += w * np.einsum('ac,bd->abcd', M_mat, M_dag)
    return V.real, basis


def _m_conservation(V, basis):
    N = len(basis)
    ok = tot = 0
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    if abs(V[a, b, c, d]) > 1e-14:
                        tot += 1
                        if basis[a].m + basis[b].m == basis[c].m + basis[d].m:
                            ok += 1
    return ok / tot if tot else 0.0


def test_resolvent_two_body_dirac_pearson():
    """Paper 31 §sec:two_body resolvent: Dirac weight w=1/(N+1/2)^2 gives
    Pearson 0.846 (n_max=2), 0.796 (n_max=3), decreasing; 100% m-conservation.

    Values updated 2026-08-29 (ERI evaluator fix: 0.81 -> 0.846, 0.75 -> 0.796).
    The NEGATIVE RESULT IS UNCHANGED -- the correlation still DECREASES with
    n_max, so the resolvent still fails to reproduce the Coulomb radial
    weights, and the mechanism (the Fock conformal factor is a function, not a
    constant) stands.  m-conservation remains exactly 1.0."""
    V2, basis2 = _resolvent_interaction(2, lambda Nq: 1.0 / (Nq + 0.5) ** 2)
    p2 = _resolvent_pearson(V2, _coulomb_tensor(basis2))
    assert p2 == pytest.approx(0.846, abs=0.02), f"Dirac Pearson n=2 {p2}"
    assert _m_conservation(V2, basis2) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.slow
def test_resolvent_two_body_dirac_pearson_nmax3_decreasing():
    """n_max=3 Dirac-resolvent Pearson ~ 0.796, decreasing vs n_max=2."""
    V2, basis2 = _resolvent_interaction(2, lambda Nq: 1.0 / (Nq + 0.5) ** 2)
    p2 = _resolvent_pearson(V2, _coulomb_tensor(basis2))
    V3, basis3 = _resolvent_interaction(3, lambda Nq: 1.0 / (Nq + 0.5) ** 2)
    p3 = _resolvent_pearson(V3, _coulomb_tensor(basis3))
    assert p3 == pytest.approx(0.796, abs=0.02), f"Dirac Pearson n=3 {p3}"
    assert p3 < p2  # decreasing with n_max


def test_resolvent_laplacian_zero_mode_decorrelation():
    """Paper 31 §sec:two_body: the Laplacian resolvent w=1/(N^2-1) is
    essentially uncorrelated with Coulomb because the N=1 zero mode of
    Delta_{S^3} excludes the dominant (1s,1s) contribution; the Dirac
    resolvent (which regularizes the zero mode) correlates far better."""
    V_lap, basis = _resolvent_interaction(
        2, lambda Nq: 1.0 / (Nq ** 2 - 1) if Nq >= 2 else 0.0)
    V_dir, _ = _resolvent_interaction(2, lambda Nq: 1.0 / (Nq + 0.5) ** 2)
    C = _coulomb_tensor(basis)
    p_lap = _resolvent_pearson(V_lap, C)
    p_dir = _resolvent_pearson(V_dir, C)
    # Laplacian resolvent essentially uncorrelated: |r_lap| << r_dir.
    # (Against the production casimir_ci Coulomb reference |r_lap| ~ 0.047,
    # consistent with the paper's "r < 0.05"; the debug/ driver's own local
    # Coulomb gives -0.131. Both confirm the decorrelation; pin the robust
    # structural fact.)
    assert abs(p_lap) < 0.15, f"Laplacian |Pearson| {p_lap}"
    assert p_dir > 0.7 > abs(p_lap)
    # the dominant (1s,1s,1s,1s) Coulomb element is killed by the zero mode
    assert abs(V_lap[0, 0, 0, 0]) < 1e-14
    assert abs(C[0, 0, 0, 0]) > 0.1
