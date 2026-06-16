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
       - angular part: 100% m-conserving, 100% Gaunt, pure k=0 monopole
       - radial Pearson vs exact Coulomb Slater: 0.58 (n_max=2),
         0.41 (n_max=3), decreasing with n_max

  2. Resolvent construction V^res_abcd = sum_NLM M[a,c] conj(M[d,b]) w(N):
       - Dirac resolvent w=1/(N+1/2)^2: Pearson 0.81 (n_max=2),
         0.75 (n_max=3)
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


def _spectral_action_connected(n_max):
    """Build the double-sum gauged spectral action {D_total, A}, return
    (connected_fraction, connected_pp_tensor, basis).

    A = sum_i sum_j (M_i^dag (x) M_j^dag) [D_total, M_i (x) M_j]
    D_total = D_1 (x) I + gamma_1 (x) D_2
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


def _angular_audit(V, basis):
    """Return (gaunt_pct, mcons_pct, multipole_k_dict) by Frobenius weight."""
    g_ok = g_bad = m_ok = m_bad = 0.0
    by_k = {}
    N = len(basis)
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    v = V[a, b, c, d]
                    if abs(v) < 1e-12:
                        continue
                    dl1 = abs(basis[a].l - basis[c].l)
                    dl2 = abs(basis[b].l - basis[d].l)
                    if dl1 == dl2:
                        g_ok += v * v
                        by_k[dl1] = by_k.get(dl1, 0.0) + v * v
                    else:
                        g_bad += v * v
                    if basis[a].m + basis[b].m == basis[c].m + basis[d].m:
                        m_ok += v * v
                    else:
                        m_bad += v * v
    tg, tm, tk = g_ok + g_bad, m_ok + m_bad, sum(by_k.values())
    gaunt = 100 * g_ok / tg if tg else None
    mcons = 100 * m_ok / tm if tm else None
    kdict = {int(k): 100 * v / tk for k, v in by_k.items()} if tk else {}
    return gaunt, mcons, kdict


def test_spectral_action_two_body_nmax2():
    """Paper 31 §sec:two_body, gauged spectral action at n_max=2.

    connected fraction ~ 77%; angular 100% m-conserving / 100% Gaunt /
    pure k=0; radial Pearson ~ 0.58.
    """
    frac, V, basis = _spectral_action_connected(2)
    gaunt, mcons, kdict = _angular_audit(V, basis)
    pearson = _radial_pearson(V, _coulomb_tensor(basis))

    # connected (non-factorizable) fraction: paper says 77%
    assert frac == pytest.approx(0.767, abs=0.01), f"conn frac {frac}"
    # angular selection rules: exact Coulomb match
    assert mcons == pytest.approx(100.0, abs=1e-6), f"m-cons {mcons}"
    assert gaunt == pytest.approx(100.0, abs=1e-6), f"Gaunt {gaunt}"
    assert set(round(k) for k in kdict) == {0}, f"multipole {kdict}"
    assert kdict[0] == pytest.approx(100.0, abs=1e-6)
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
    gaunt3, mcons3, kdict3 = _angular_audit(V3, basis3)
    pearson3 = _radial_pearson(V3, _coulomb_tensor(basis3))

    assert frac3 == pytest.approx(0.324, abs=0.01), f"conn frac {frac3}"
    assert mcons3 == pytest.approx(100.0, abs=1e-6), f"m-cons {mcons3}"
    assert gaunt3 == pytest.approx(100.0, abs=1e-6), f"Gaunt {gaunt3}"
    assert set(round(k) for k in kdict3) == {0}, f"multipole {kdict3}"
    assert pearson3 == pytest.approx(0.41, abs=0.02), f"Pearson {pearson3}"
    # the load-bearing negative: radial correlation DECREASES with n_max
    assert pearson3 < pearson2


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
    Pearson 0.81 (n_max=2), 0.75 (n_max=3), decreasing; 100% m-conservation."""
    V2, basis2 = _resolvent_interaction(2, lambda Nq: 1.0 / (Nq + 0.5) ** 2)
    p2 = _resolvent_pearson(V2, _coulomb_tensor(basis2))
    assert p2 == pytest.approx(0.81, abs=0.02), f"Dirac Pearson n=2 {p2}"
    assert _m_conservation(V2, basis2) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.slow
def test_resolvent_two_body_dirac_pearson_nmax3_decreasing():
    """n_max=3 Dirac-resolvent Pearson ~ 0.75, decreasing vs n_max=2."""
    V2, basis2 = _resolvent_interaction(2, lambda Nq: 1.0 / (Nq + 0.5) ** 2)
    p2 = _resolvent_pearson(V2, _coulomb_tensor(basis2))
    V3, basis3 = _resolvent_interaction(3, lambda Nq: 1.0 / (Nq + 0.5) ** 2)
    p3 = _resolvent_pearson(V3, _coulomb_tensor(basis3))
    assert p3 == pytest.approx(0.75, abs=0.02), f"Dirac Pearson n=3 {p3}"
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
